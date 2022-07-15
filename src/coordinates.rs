mod generate_database;
mod flower;
#[cfg(test)]
mod tests;
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
    Ok(FlowerPattern::generate(0, k, fov, &stars_to_sort)?)
}
/// returns the index of the catalogue star, starting at 0, which best matches the observed flower pattern, and the
/// offset tau, for which r'(i) = r((i - tau) % k)
/// R_rs - a Vec of the same length of the star catalogue, each entry is the R vector corresponding
/// to the functions r(i) and r'(i) for every star from the catalogue and the central star
/// R_deltas - same with R_rs but for delta(i) and delta'(i)
/// if successful, returns Ok(index of the best matching star in the catalogue)
/// if not, it could not match sufficiently well any star on the catalogue to the central star,
/// possibly indicating a false central star, or many false stars as petels
/// TODO make a better choice function
fn match_catalogue_star_to_central(R_rs: &mut Vec<Vec<Complex<f64>>>, R_deltas: &mut Vec<Vec<Complex<f64>>>) -> Result<(u16, u16), String> {
    let p: u16 = 20;
    let best_matches_r = top_match_candidates(R_rs, p)?;
    let best_matches_delta = top_match_candidates(R_deltas, p)?;
    
    // a vector of (central_star_index, tau) elements which are in both the best matches for R and
    // delta
    let mut matching_star: Vec<(u16, u16)> = vec![];
    for i in 0..p as usize{
        for j in 0..p as usize{
            if best_matches_r[i] == best_matches_delta[j] {
                matching_star.push(best_matches_r[i]);
            } 
        }
    }
    if matching_star.len() == 0 {
        return Err("there are no stars which match both the r'(i) and delta'(i) of the observed flower pattern".to_string());
    }
    // this is the match that is the best according to distance from central star measurements
    // (r'(i)). the matching_star vector is sorted by the score of matches from r'(i)
    // most of the times the matching_star vec will have one element anyways(i hope)
    Ok(matching_star[0])
}
/// this function looks at an R vector (as in wiki article) and returns the top p candidates for
/// matches, as per some score function inside this one
/// returns - a result of Vector of tuples (star_index, tau), the vector has length p. The vector
/// is sorted by highest score first
/// if it does not find enough good matches or deems the best matches too bad it returns an error
/// with the corresponding message 
/// R needs to be mutable to allow for inplace inverse fft, using which it determines the best
/// matching stars
fn top_match_candidates(R: &mut Vec<Vec<Complex<f64>>>, p: u16) -> Result<Vec<(u16, u16)>, String> {
    let n = R.len();
    if n < p as usize {
        return Err(format!("The passed value for R vector has less elements than the required {} top matches ({} < {})", 
                p, n, p));
    }
    let mut planner = FftPlanner::new();
    let k = R[0].len();
    let inverse_fft = planner.plan_fft_inverse(k);

    #[derive(Copy, Clone)]
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
        inverse_fft.process(&mut R[index]);
        let (score, tau) = get_score_tau(&R[index]);
        matches.push(BestMatch {
            score,
            central_star_index: index as u16, 
            tau
        });
    }
    // TODO make finding the the top p best matches happen in the upper loop, for better efficiency
    let mut top_p_matches: Vec<BestMatch> = vec![];
    for _ in 0..p {
        let (best_match_index, &best_match) = matches.iter().enumerate().max_by(|&(_, m1), &(_, m2)| {
            m1.score.partial_cmp(&(m2.score)).unwrap()
        }).unwrap();
        top_p_matches.push(best_match);
        // so the next iteration does not pick up the same maximum in score
        matches.remove(best_match_index);
    }

    let final_result: Vec<(u16, u16)> = top_p_matches.iter().map(|&bm| {
        (bm.central_star_index, (k as u16 - bm.tau)%k as u16)
    }).collect(); 

    Ok(final_result)
}
/// determines how much the camera is rotated counter-clockwise(so return value is always in the range [0, 2pi)) relative to a sky coordinate system (see README)
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
