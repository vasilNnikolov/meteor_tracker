use crate::parse_stars::star::Star;
use crate::coordinates::flower;
use std::fs::File;
use std::io::{BufReader, prelude::*};
use rustfft::{FftPlanner, num_complex::Complex};

// extern crate serde_json;
// extern crate serde;
extern crate serde_derive;

#[derive(Debug, serde_derive::Serialize, serde_derive::Deserialize)]
pub struct DFT_coeffs {
    r_coeffs: Vec<Complex<f64>>, 
    delta_coeffs: Vec<Complex<f64>>
}

/// tries to read HYG db, if successful returns a vector of stars
/// note: the indexes of the stars in the Star struct start at 1 
/// m_lim - the limiting magnitude of the camera, and tha highest magnitude in the resulting vector
pub fn read_hyd_database(m_lim: f64) -> Result<Vec<Star>, String> {
    let file_path = "data/hyg_data.csv";
    let file = match File::open(file_path) {
        Ok(f) => f, 
        Err(_) => return Err(String::from("could not open hyg database file"))
    };

    let buff_reader = BufReader::new(file);
    let mut lines = buff_reader.lines().map(|l| l.unwrap());
    lines.next(); // the heading 

    let mut stars: Vec<Star> = vec![];
    for (index, line) in lines.enumerate() {
        let line_contents = line.split(',').collect::<Vec<&str>>();
        if line_contents.len() != 4 as usize {
            return Err("The hyg star db has a line with more or less than 3 entries".to_string());
        } 

        let mut line_contents_float: [f64; 4] = [0.0; 4];

        for i in 0..4 {
            let x = line_contents[i].parse::<f64>();
            match x {
                Err(_) => return Err(format!("Bad data on row {} of HYG database", index)), 
                Ok(f) => line_contents_float[i] = f
            }
        }
        let [index, ra, dec, brightness] = line_contents_float;

        if brightness > m_lim {
            break;
        }
        stars.push(Star::new(
                ra*15.0*std::f64::consts::PI/180.0, // ra in the database is in hours
                dec*std::f64::consts::PI/180.0, 
                brightness, 
                index as u16));
    }
    Ok(stars)
}

/// Generates a file of DFT coefficients for flower patterns. 
/// m_lim - the limiting magnitude of the camera
/// k - the number of outer stars in the flower pattern, approx 10, maybe more
/// fov - the radius of the field of view of the camera, in radians
/// It creates a file with N entries, N being the number of stars brighter than m_lim, and for each one stores the DFT coefficients
/// of r(i) - distance between the central star and star number i of the petels, and delta(i) - angle between the petel i and (i + 1)
pub fn generate_db(m_lim: f64, k: u16, fov: f64) -> Result<(), String>{
    let all_stars = read_hyd_database(m_lim)?; 

    let mut fft_planner = FftPlanner::new();
    let fft = fft_planner.plan_fft_forward(k as usize);

    // create file to write to 
    let db_location = format!("databases/k{}fov_deg{}.json", k, (fov*180.0/std::f64::consts::PI) as u16);
    let db_file = match File::create(db_location) {
        Ok(f) => f, 
        Err(_) => { return Err(String::from("could not create file to write the dft coefficients to")); }
    };
    
    for star in all_stars.iter() {
        let pattern = flower::FlowerPattern::generate(star.index, k, fov, &all_stars)?;
        // these will hold the DFT coefficients after the transformation
        let mut r_values: Vec<Complex<f64>> = pattern.r.iter().map(|&value_f| Complex::new(value_f, 0.0)).collect();
        let mut delta_values: Vec<Complex<f64>> = pattern.delta.iter().map(|&value_f| Complex::new(value_f, 0.0)).collect();
        // generate DFT coefficients of r(i) and delta(i) functions, i in [1..k]
        fft.process(&mut r_values);
        fft.process(&mut delta_values);
        // write DFT coefficients of flower patterns to .csv file
        let coeff_struct = DFT_coeffs {
            r_coeffs: r_values, 
            delta_coeffs: delta_values
        };
        
        
    }

    Ok(())
}
