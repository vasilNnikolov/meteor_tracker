use crate::parse_stars::star::Star;
use std::fs::File;
use std::io::{BufReader, prelude::*};

fn read_hyd_database() -> Result<Vec<Star>, String> {
    let file_path = "data/hyg_data.csv";
    let file = match File::open(file_path) {
        Ok(f) => f, 
        Err(_) => return Err(String::from("could not open hyg database file"))
    };

    let buff_reader = BufReader::new(file);
    let mut iter = buff_reader.lines().map(|l| l.unwrap()).enumerate();
    iter.next(); // the heading 

    let mut stars: Vec<Star> = vec![];
    for (index, line) in iter{
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

        stars.push(Star::new(
                ra*std::f64::consts::PI/180.0, 
                dec*std::f64::consts::PI/180.0, 
                brightness, 
                index as i32));
    }
    Ok(stars)
}

/// Generates a file of DFT coefficients for flower patterns. 
/// m_lim - the limiting magnitude of the camera
/// k - the number of outer stars in the flower pattern, approx 10, maybe more
/// It creates a file with N entries, N being the number of stars brighter than m_lim, and for each one stores the DFT coefficients
/// of r(i) - distance between the central star and star number i of the petels, and delta(i) - angle between the petel i and (i + 1)
pub fn generate_db(m_lim: f64, k: u16) -> Result<(), String>{
    let all_stars = read_hyd_database()?; 
    let mut brightest_n_stars: Vec<Star> = vec![];

    for s in all_stars {
        if s.brightness <= m_lim {
            brightest_n_stars.push(s);
        } else { break; }
    }
    // TODO for each star generate the flower pattern, generate their DFTs and write them to
    // database

    Ok(())
}
