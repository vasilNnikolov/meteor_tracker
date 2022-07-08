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
    let k = dft_database.r_dft_coefficients[0].len();
    let fov = dft_database.fov;
    let observed_pattern = generate_flower_pattern_from_observation(captured_stars, k as u16, fov)?;

    let mut fft_planner = FftPlanner::new();
    let fft = fft_planner.plan_fft_forward(k as usize);
    let mut r_values: Vec<Complex<f64>> = observed_pattern.r.iter().map(|&value_f| Complex::new(value_f, 0.0)).collect();
    let mut delta_values: Vec<Complex<f64>> = observed_pattern.delta.iter().map(|&value_f| Complex::new(value_f, 0.0)).collect();

    // generate DFT coefficients of r(i) and delta(i) functions, i in [1..k]
    fft.process(&mut r_values);
    fft.process(&mut delta_values);
    
    // pass those to the match_catalogue_star_to_central function
    let mut R_rs: Vec<Vec<Complex<f64>>> = vec![];
    let mut R_deltas: Vec<Vec<Complex<f64>>> = vec![];

    // get size of the catalogue
    let n = dft_database.r_dft_coefficients.len();
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

    // angle of petel minus angle of corresponding petel in the catalogue flower pattern
    let mut average_angle_offset = 0.0;
    let catalogue_flower_pattern = flower_patterns.get(central_star_true_index as usize).unwrap(); 
    for petel_index in 0..k as i32 {
        let corresponding_index_of_observed_star = ((petel_index - tau as i32) % k as i32 + k as i32) % k as i32;
        average_angle_offset += observed_pattern.angle_of_petel(petel_index as u16).unwrap() - catalogue_flower_pattern.angle_of_petel(corresponding_index_of_observed_star as u16).unwrap();
    }
    average_angle_offset /= k as f64;

    // in geocentric coordinate system
    let central_star_coords = flower_patterns.get(central_star_true_index as usize).unwrap().central_star.coords; 
    // matrix magic begins, please work
    // build T_rc
    let y_rc = central_star_coords; 
    let mut x_rc = y_rc.cross(&Vector3::new(0.0, 0.0, 1.0));
    x_rc = x_rc/x_rc.norm();
    let z_rc= x_rc.cross(&y_rc);

    let T_rc = Matrix3::from_columns(&[x_rc, y_rc, z_rc]); 

    // build Rot_y(-average_angle_offset)
    let x_roty = Vector3::new(average_angle_offset.cos(), 0.0, -(average_angle_offset.sin()));
    let y_roty = Vector3::new(0.0, 1.0, 0.0);
    let z_roty = x_roty.cross(&y_roty);

    let Rot_y = Matrix3::from_columns(&[x_roty, y_roty, z_roty/z_roty.norm()]); 

    // build T_rc' inverse
    let y_rc_prime = observed_pattern.central_star.coords; 
    let mut x_rc_prime = y_rc_prime.cross(&Vector3::new(0.0, 0.0, 1.0));
    x_rc_prime = x_rc_prime/x_rc_prime.norm();
    let z_rc_prime= x_rc_prime.cross(&y_rc_prime);

    let T_rc_prime_inverse = Matrix3::from_columns(&[x_rc_prime, y_rc_prime, z_rc_prime]).try_inverse().unwrap(); 

    let R = T_rc*Rot_y*T_rc_prime_inverse;

    Ok(R)
}

/// generates a flower pattern from the observed stars, to be used for matching against the DFT
/// database 
/// TODO create better implementation
/// currently it selects the (k + 1) brightest stars, picks the one closest to the center of the
/// image as the central one, and the other k as petels
fn generate_flower_pattern_from_observation(captured_stars: &Vec<star::Star>, k: u16, fov: f64) -> Result<FlowerPattern, String> {
    let mut k_plus_one_brightest: Vec<star::Star> = captured_stars.clone();
    // sort by lowest magnitude first
    k_plus_one_brightest.sort_by(|&a, &b| {
        a.brightness.partial_cmp(&b.brightness).unwrap()
    });
    if k + 1 > k_plus_one_brightest.len() as u16 {
        return Err(format!("there are not enough stars in the image to generate a flower pattern of size k={}", k));
    }
    k_plus_one_brightest = k_plus_one_brightest.into_iter().take((k + 1) as usize).collect();
    // sort them by distance to the center of the image, smallest distance first
    k_plus_one_brightest.sort_by(|&a, &b| {
        let a_dist_squared = a.coords.x.powi(2) + a.coords.z.powi(2);
        let b_dist_squared = b.coords.x.powi(2) + b.coords.z.powi(2);
        a_dist_squared.partial_cmp(&b_dist_squared).unwrap()
    });
    println!("k_plus_one_brightest: {}", k_plus_one_brightest.len());
    // index is 1 and not 0 because of the implementation of generate, made for index according to
    // the HYG db
    Ok(FlowerPattern::generate(1, k, fov, &k_plus_one_brightest)?)
}
/// returns the index of the catalogue star which best matches the observed flower pattern, and the
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

    // the inverse fft should have a value equal to K at some i, and close to 0 everywhere else
    // currently we only look at distance data, not angle data
    // TODO do the same process for the delta data, and compare the results
    let mut tau_r: u16 = 0;
    let mut best_highest_r: f64 = 0.0;
    let mut best_matching_star_index_r: u16 = 0;

    for index in 0..n {
        inverse_fft.process(&mut R_rs[index]);
        let (max_r_index, max_r_value) = R_rs[index].iter().enumerate()
            .map(|(i, &c)| (i, c.norm()))
            .max_by(|&(_, f1), &(_, f2)| f1.partial_cmp(&f2).unwrap_or(std::cmp::Ordering::Equal)).unwrap();

        if max_r_value > best_highest_r {
            best_highest_r = max_r_value;
            tau_r = max_r_index as u16;
            best_matching_star_index_r = index as u16;
        }
    }

    Ok((best_matching_star_index_r, (k as u16 - tau_r)%k as u16))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_file_read() {
        let all_stars = generate_database::read_hyd_database(4.0).unwrap();
        for i in 0..5 {
            println!("{:?}", all_stars[i]);
        }
        let (ra, dec) = all_stars[0].get_ra_dec();
        println!("sirius ra {}, dec {}, brightness {}", ra*57.295, dec*57.295, all_stars[0].brightness);
    }
    #[test]
    fn test_sirius_pattern() {
        let all_stars = generate_database::read_hyd_database(4.0).unwrap();
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
        let all_stars = generate_database::read_hyd_database(4.0).unwrap();
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
    fn test_db_generation() {
        if let Err(s) = generate_database::generate_db(4.0, 5, 0.52) {
            println!("{}", s);
            panic!();
        }
    }
    #[test]
    fn test_main_fun_of_module(){
        // camera is facing in the direction R_y
        let R_x = Vector3::new(0.9, 0.0, 0.43588);
        let mut R_y = R_x.cross(&Vector3::new(0.0, 1.0, 0.0));
        R_y = R_y / R_y.norm();
        let R_z = R_x.cross(&R_y);
        let R = Matrix3::from_columns(&[R_x, R_y, R_z]);

        // generate list of captured stars
        let R_inv = R.try_inverse().unwrap();
        let fov = 0.52;
        let k = 5;
        let (DFT_db, flower_patterns) = generate_database::generate_db(4.0, k, fov).unwrap();
        let all_stars: Vec<star::Star> = flower_patterns.iter().map(|fp| fp.central_star).collect();

        let mut captured_stars: Vec<star::Star> = all_stars
            .iter()
            .filter(|&s| s.coords.dot(&R_y) > fov.cos())
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
}
