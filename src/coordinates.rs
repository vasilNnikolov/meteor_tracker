pub mod generate_database;
pub mod flower;
use crate::parse_stars::star;
use rustfft::{FftPlanner, num_complex::Complex};


/// gets the vector of captured stars and a dft database which is pre-calculated at the beginning
/// of the observations, and tries to compute the rotation matrix of the camera
/// the return matrix R is such that 
/// R*x_camera = x_real, where x_camera is the radius-vector of a point in the camera coordinate
/// system and x_real is the position of the same point in the geocentric coordinate system
pub fn get_rotation_matrix(captured_stars: &Vec<star::Star>, dft_database: generate_database::DFT_database) -> Result<nalgebra::Matrix3<f64>, String> {
    let k = dft_database.r_dft_coefficients[0].len();
    let fov = dft_database.fov;
    let pattern = flower::FlowerPattern::generate(6969, k as u16, fov, captured_stars)?;

    let mut fft_planner = FftPlanner::new();
    let fft = fft_planner.plan_fft_forward(k as usize);
    let mut r_values: Vec<Complex<f64>> = pattern.r.iter().map(|&value_f| Complex::new(value_f, 0.0)).collect();
    let mut delta_values: Vec<Complex<f64>> = pattern.delta.iter().map(|&value_f| Complex::new(value_f, 0.0)).collect();

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
    let central_star_true_index: u16;
    match match_catalogue_star_to_central(R_rs, R_deltas) {
        Ok(index) => central_star_true_index = index, 
        Err(s) => {
            // TODO handle error of no matched star better
            return Err(s);
        }
    }

    // once we know the central star index, we need the rotation of the camera relative to the Sky
    // Coordinate system (see README) by checking the actual angles of the observed flower pattern
    // and comparing them with the expected angles for the central star. The observed angles should
    // all be shifted some amount
    // TODO implement ^
    
    Err(String::from("get_rotation_matrix not unimplemented!"))
}

/// returns the index of the catalogue star which best matches the observed flower pattern
/// R_rs - a Vec of the same length of the star catalogue, each entry is the R vector corresponding
/// to the functions r(i) and r'(i) for every star from the catalogue and the central star
/// R_deltas - same with R_rs but for delta(i) and delta'(i)
/// if successful, returns Ok(index of the best matching star in the catalogue)
/// if not, it could not match sufficiently well any star on the catalogue to the central star,
/// possibly indicating a false central star, or some many false stars as petels
fn match_catalogue_star_to_central(R_rs: Vec<Vec<Complex<f64>>>, R_deltas: Vec<Vec<Complex<f64>>>) -> Result<u16, String> {
    Err(String::from("match_catalogue_star_to_central not unimplemented!"))
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
        let sirius_pattern = flower::FlowerPattern::generate(1, 5, 0.35, &all_stars).unwrap();
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
        let canopus_pattern = flower::FlowerPattern::generate(2, 5, 0.35, &all_stars).unwrap();
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
}
