pub mod generate_database;
pub mod flower;
use crate::parse_stars::star;
use crate::coordinates::flower::FlowerPattern;
use rustfft::{FftPlanner, num_complex::Complex};
use nalgebra::{Vector3, Matrix3};

/// gets the vector of captured stars and a dft database which is pre-calculated at the beginning
/// of the observations, and tries to compute the rotation matrix of the camera
/// the return matrix R is such that 
/// R*x_camera = x_real, where x_camera is the radius-vector of a point in the camera coordinate
/// system and x_real is the position of the same point in the geocentric coordinate system
pub fn get_rotation_matrix(captured_stars: &Vec<star::Star>, dft_database: &generate_database::DFT_database, flower_patterns: &Vec<FlowerPattern>) -> Result<nalgebra::Matrix3<f64>, String> {
    // get size of the catalogue
    let n = dft_database.r_dft_coefficients.len();
    let k = dft_database.r_dft_coefficients[0].len();
    let fov = dft_database.fov;

    let observed_pattern = generate_flower_pattern_from_observation(captured_stars, k as u16, fov)?;

    // generate DFT coefficients of r(i) and delta(i) functions, i in [1..k]
    let mut fft_planner = FftPlanner::new();
    let fft = fft_planner.plan_fft_forward(k as usize);
    let mut r_values: Vec<Complex<f64>> = observed_pattern.r.iter().map(|&value_f| Complex::new(value_f, 0.0)).collect();
    let mut delta_values: Vec<Complex<f64>> = observed_pattern.delta.iter().map(|&value_f| Complex::new(value_f, 0.0)).collect();
    fft.process(&mut r_values);
    fft.process(&mut delta_values);
    
    // pass those to the match_catalogue_star_to_central function
    // those are the Rs as in the wiki article
    let mut R_rs: Vec<Vec<Complex<f64>>> = vec![];
    let mut R_deltas: Vec<Vec<Complex<f64>>> = vec![];

    for star_index in 0..n {
        // the R vectors of complex numbers as in the wikipedia page for phase correlation
        let mut R_r: Vec<Complex<f64>> = vec![];
        let mut R_delta: Vec<Complex<f64>> = vec![];
        for i in 0..k {
            let (g_a_r, g_b_r) = (dft_database.r_dft_coefficients[star_index][i], r_values[i]);
            R_r.push(g_a_r*g_b_r.conj()/(g_a_r*g_b_r.conj()).norm());

            let (g_a_delta, g_b_delta) = (dft_database.delta_dft_coefficients[star_index][i], delta_values[i]);
            R_delta.push(g_a_delta*g_b_delta.conj()/(g_a_delta*g_b_delta.conj()).norm());
        }
        R_rs.push(R_r);
        R_deltas.push(R_delta);
    } 

    let (central_star_true_index, tau) = match match_catalogue_star_to_central(&mut R_rs, &mut R_deltas) {
        Ok((index, t)) => (index, t), 
        Err(s) => {
            // TODO handle error of no matched star better
            return Err(s);
        }
    };

    // angle of petel minus angle of corresponding petel in the catalogue flower pattern, it is
    // always positive, so in range [0, 2pi)
    let catalogue_flower_pattern = flower_patterns.get(central_star_true_index as usize).unwrap();
    let average_angle_offset = get_angle_of_relative_rotation(&observed_pattern, &catalogue_flower_pattern, tau);

    let R = construct_rot_matrix_from_data(&observed_pattern, &catalogue_flower_pattern, average_angle_offset);

    Ok(R)
}
/// generates a flower pattern from the observed stars, to be used for matching against the DFT
/// database 
/// TODO create better implementation
/// currently it selects the (k + 1) brightest stars, picks the one closest to the center of the
/// image as the central one, and the other k as petels
fn generate_flower_pattern_from_observation(captured_stars: &Vec<star::Star>, k: u16, fov: f64) -> Result<FlowerPattern, String> {
    if k + 1 > captured_stars.len() as u16 {
        return Err(format!("there are not enough stars in the image to generate a flower pattern of size k={}", k));
    }
    let mut stars_to_sort: Vec<star::Star> = (*captured_stars).clone();
    // sort them by distance to the center of the image, smallest distance first
    let center_of_image = Vector3::new(0.0, 1.0, 0.0);
    stars_to_sort.sort_by(|&a, &b| {
        b.coords.dot(&center_of_image).partial_cmp(&(a.coords.dot(&center_of_image))).unwrap()
    });
    // index is 1 and not 0 because of the implementation of generate, made for index according to
    // the HYG db
    Ok(FlowerPattern::generate(1, k, fov, &stars_to_sort)?)
}
/// returns the index of the catalogue star (starting at 0, not at 1 as the catalogue itself) which best matches the observed flower pattern, and the
/// offset tau, for which r'(i) = r((i - tau) % k)
/// R_rs - a Vec of the same length of the star catalogue, each entry is the R vector corresponding
/// to the functions r(i) and r'(i) for every star from the catalogue and the central star
/// R_deltas - same with R_rs but for delta(i) and delta'(i)
/// if successful, returns Ok(index of the best matching star in the catalogue)
/// if not, it could not match sufficiently well any star on the catalogue to the central star,
/// possibly indicating a false central star, or many false stars as petels
/// TODO make a better choice function
fn match_catalogue_star_to_central(R_rs: &mut Vec<Vec<Complex<f64>>>, R_deltas: &mut Vec<Vec<Complex<f64>>>) -> Result<(u16, u16), String> {
    let n = R_rs.len();
    let mut planner = FftPlanner::new();
    let k = R_rs[0].len();
    let inverse_fft = planner.plan_fft_inverse(k);

    struct BestMatch {
        score: f64, 
        central_star_index: u16, 
        tau: u16
    }
    
    fn get_score_tau(r: &Vec<Complex<f64>>) -> (f64, u16) {
        let (max_r_index, max_r_value): (usize, f64) = r.iter().enumerate()
            .map(|(i, c)| (i, c.norm()))
            .max_by(|&(_, f1), &(_, f2)| f1.partial_cmp(&f2).unwrap()).unwrap();

        let k = r.len();
        let residual_quadratic = r.iter().map(|&c| c.norm().powi(2)).sum::<f64>() - max_r_value.powi(2);
        let score: f64 = 1.0 - (residual_quadratic/(k - 1) as f64).powf(0.5)/max_r_value;
        (score, max_r_index as u16) 
    }

    let mut matches: Vec<BestMatch> = vec![];
    for index in 0..n {
        inverse_fft.process(&mut R_rs[index]);
        let (score, tau) = get_score_tau(&R_rs[index]);
        matches.push(BestMatch {
            score,
            central_star_index: index as u16, 
            tau
        });
    }
    // find best score 
    let best_match = matches.iter().max_by(|&m1, &m2| {
        m1.score.partial_cmp(&(m2.score)).unwrap()
    }).unwrap();

    println!("score of match: {}", best_match.score);

    Ok((best_match.central_star_index, (k as u16 - best_match.tau)%k as u16))
}
/// determines how much the camera is rotated counter-clockwise relative to a sky coordinate system (see README)
/// centered at the determined central star
/// observed_pattern - the observed flower pattern with star coordinates in the camera coordinate
/// system
/// catalogue_flower_pattern - the FP of the star which was determined to be central, coordinates
/// are in the geocentric CS
/// tau - the offset of the indexes of r'(i) and r(i), such that r'(i) = r((i - tau)%k)
fn get_angle_of_relative_rotation(observed_pattern: &FlowerPattern, catalogue_flower_pattern: &FlowerPattern, tau: u16) -> f64 {
    let k = observed_pattern.r.len();
    let mut average_angle_offset = 0.0;
    for petel_index in 0..k as i32 {
        let corresponding_index_of_observed_star = (petel_index - tau as i32 + k as i32) % k as i32;
        let observed_angle = observed_pattern.angle_of_petel(petel_index as u16).unwrap();
        let catalogue_angle = catalogue_flower_pattern.angle_of_petel(corresponding_index_of_observed_star as u16).unwrap();
        let delta_angle = observed_angle - catalogue_angle;
        if delta_angle < 0.0 {
            average_angle_offset += delta_angle + 2.0*std::f64::consts::PI;
        } else {
            average_angle_offset += delta_angle;
        }
    }
    average_angle_offset / k as f64
}
/// after we have determined the best match for the central star and the rotation of the camera
/// relative to the sky coordinate system centered at that star, this fn performs the matrix magic
/// to compute the final rotation matrix of the camera, such that r_geocentric = R*r_camera
fn construct_rot_matrix_from_data(observed_pattern: &FlowerPattern, catalogue_flower_pattern: &FlowerPattern, average_angle_offset: f64) -> Matrix3<f64> {
    // matrix magic begins
    // build T_rc
    let y_rc = catalogue_flower_pattern.central_star.coords; 
    let mut x_rc = y_rc.cross(&Vector3::new(0.0, 0.0, 1.0));
    x_rc = x_rc/x_rc.norm();
    let z_rc= x_rc.cross(&y_rc);

    let T_rc = Matrix3::from_columns(&[x_rc, y_rc, z_rc]); 

    // build Rot_y(-average_angle_offset)
    let x_roty = Vector3::new(average_angle_offset.cos(), 0.0, -(average_angle_offset.sin()));
    let y_roty = Vector3::new(0.0, 1.0, 0.0);
    let z_roty = x_roty.cross(&y_roty);

    let Rot_y = Matrix3::from_columns(&[x_roty, y_roty, z_roty]); 

    // build T_rc' inverse
    let y_rc_prime = observed_pattern.central_star.coords; 
    let mut x_rc_prime = y_rc_prime.cross(&Vector3::new(0.0, 0.0, 1.0));
    x_rc_prime = x_rc_prime/x_rc_prime.norm();
    let z_rc_prime= x_rc_prime.cross(&y_rc_prime);

    let T_rc_prime_inverse = Matrix3::from_columns(&[x_rc_prime, y_rc_prime, z_rc_prime]).try_inverse().unwrap(); 

    T_rc*Rot_y*T_rc_prime_inverse
}
fn print_flower_pattern(fp: &FlowerPattern) {
    let n = fp.outer_stars.len(); 
    let (ra, dec) = fp.central_star.get_ra_dec();
    println!("central star ra, dec: {:.2} h {:.2} deg, magnitude {:.2}", 
        ra*57.3/15.0, 
        dec*57.3, 
        fp.central_star.brightness);
    for i in 0..n {
        let (ra, dec) = fp.outer_stars[i].get_ra_dec();
        println!("ra, dec: {:.2} h {:.2} deg, magnitude {:.2}, distance {:.2}, angle {:.2}", 
            ra*57.3/15.0, 
            dec*57.3, 
            fp.outer_stars[i].brightness,
            fp.r[i]*57.295, 
            flower::angle_of_outer_petel(&fp.central_star, &fp.outer_stars[i])*57.295);
        println!("angle between star {} and {}: {}", i, (i + 1)%n, fp.delta[i]*57.295);
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    extern crate rand;
    use rand::Rng;
    #[test]
    fn test_file_read() {
        let all_stars = generate_database::read_hyd_database(4.0, 0.0017).unwrap();
        for i in 0..5 {
            println!("{:?}", all_stars[i]);
        }
        let (ra, dec) = all_stars[0].get_ra_dec();
        println!("sirius ra {}, dec {}, brightness {}", ra*57.295, dec*57.295, all_stars[0].brightness);
    }
    #[test]
    fn test_sirius_pattern() {
        let all_stars = generate_database::read_hyd_database(4.0, 0.0017).unwrap();
        let sirius_pattern = FlowerPattern::generate(1, 5, 0.35, &all_stars).unwrap();
        let n = sirius_pattern.outer_stars.len(); 
        for i in 0..n {
            let (ra, dec) = sirius_pattern.outer_stars[i].get_ra_dec();
            println!("ra, dec: {} {}, magnitude {}, distance {}, angle {}", 
                ra*57.3, 
                dec*57.3, 
                sirius_pattern.outer_stars[i].brightness,
                sirius_pattern.r[i]*57.295, 
                flower::angle_of_outer_petel(&all_stars[0], &sirius_pattern.outer_stars[i])*57.295);
            println!("angle between star {} and {}: {}", i, (i + 1)%n, sirius_pattern.delta[i]*57.295);
        }
    }
    #[test]
    fn test_canopus_pattern() {
        let all_stars = generate_database::read_hyd_database(4.0, 0.0017).unwrap();
        let canopus_pattern = FlowerPattern::generate(2, 5, 0.35, &all_stars).unwrap();
        let n = canopus_pattern.outer_stars.len(); 
        for i in 0..n {
            let (ra, dec) = canopus_pattern.outer_stars[i].get_ra_dec();
            println!("ra, dec: {} {}, magnitude {}, distance {}, angle {}", 
                ra*57.3, 
                dec*57.3, 
                canopus_pattern.outer_stars[i].brightness,
                canopus_pattern.r[i]*57.295, 
                flower::angle_of_outer_petel(&all_stars[1], &canopus_pattern.outer_stars[i])*57.295);
            println!("angle between star {} and {}: {}", i, (i + 1)%n, canopus_pattern.delta[i]*57.295);
        }
    }
    #[test]
    fn test_polaris_pattern() {
        let all_stars = generate_database::read_hyd_database(4.0, 0.0017).unwrap();
        let canopus_pattern = FlowerPattern::generate(46, 5, 0.35, &all_stars).unwrap();
        let (ra, dec) = canopus_pattern.central_star.get_ra_dec();
        println!("ra, dec of polaris: {:.2}h, {:.2}deg", ra*57.29578/15.0, dec*57.29577);
        let n = canopus_pattern.outer_stars.len(); 
        for i in 0..n {
            let (ra, dec) = canopus_pattern.outer_stars[i].get_ra_dec();
            println!("ra, dec: {:.2} {:.2}, magnitude {:.2}, distance {:.2}, angle {:.2}", 
                ra*57.3, 
                dec*57.3, 
                canopus_pattern.outer_stars[i].brightness,
                canopus_pattern.r[i]*57.295, 
                flower::angle_of_outer_petel(&all_stars[45], &canopus_pattern.outer_stars[i])*57.295);
            println!("angle between star {} and {}: {:.2}", i, (i + 1)%n, canopus_pattern.delta[i]*57.295);
        }
    }
    #[test]
    fn test_db_generation() {
        if let Err(s) = generate_database::generate_db(4.0, 5, 0.52) {
            println!("{}", s);
            panic!();
        }
    }
    fn generate_random_orthogonal_matrix() -> Matrix3<f64> {
        let mut rng = rand::thread_rng();
        let mut x = Vector3::new(rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0));
        x = x / x.norm();
        let mut y = Vector3::new(0.0, rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0)).cross(&x);
        y = y / y.norm();
        let z = x.cross(&y);
        Matrix3::from_columns(&[x, y, z])
    }
    #[test]
    fn test_main_fun_of_module() {
        // camera is facing in the direction R_y
        let R = generate_random_orthogonal_matrix();

        // let R = Matrix3::<f64>::identity();
        // generate list of captured stars
        let R_inv = R.try_inverse().unwrap();
        let fov = 0.52;
        let k = 10;
        let (DFT_db, flower_patterns) = generate_database::generate_db(4.0, k, fov).unwrap();
        let all_stars: Vec<star::Star> = flower_patterns.iter().map(|fp| fp.central_star).collect();

        let R_y = R.column(1);
        let captured_stars: Vec<star::Star> = all_stars
            .iter()
            .filter(|&s| s.coords.dot(&R_y) > (2.0*fov).cos())
            .map(|&s| star::Star {
                index: 6969, 
                coords: R_inv*s.coords,
                brightness: s.brightness
            })
            .collect();
        println!("captured {} stars in fov", captured_stars.len());

        let experiment_R = get_rotation_matrix(&captured_stars, &DFT_db, &flower_patterns).unwrap();
        println!("exp: {:?}", experiment_R);
        println!("real: {:?}", R);
    }
    #[test]
    fn long_test() {
        let fov = 0.52;
        let k = 10;
        let (DFT_db, flower_patterns) = generate_database::generate_db(4.0, k, fov).unwrap();
        let all_stars: Vec<star::Star> = flower_patterns.iter().map(|fp| fp.central_star).collect();
        for _ in 0..1000 {
            let R = generate_random_orthogonal_matrix();
            let R_inv = R.try_inverse().unwrap();

            let R_y = R.column(1);
            let captured_stars: Vec<star::Star> = all_stars
                .iter()
                .filter(|&s| s.coords.dot(&R_y) > (1.5*fov).cos())
                .map(|&s| star::Star {
                    index: 6969, 
                    coords: R_inv*s.coords,
                    brightness: s.brightness
                })
                .collect();

            let experiment_R = get_rotation_matrix(&captured_stars, &DFT_db, &flower_patterns).unwrap();
        }
    }
    #[test]
    fn min_distance_between_stars() {
        let all_stars = generate_database::read_hyd_database(4.0, 0.0017).unwrap(); 
        let mut max_cos_angle: f64 = -1.0;
        for i in 0..all_stars.len() {
            for j in 0..all_stars.len() {
                let angle = star::cos_angle_between_stars(&all_stars[i], &all_stars[j]);
                if max_cos_angle < angle && i != j {
                    max_cos_angle = angle;
                }
            }
        }
        println!("max_cos_angle: {}", max_cos_angle);
    }
}
