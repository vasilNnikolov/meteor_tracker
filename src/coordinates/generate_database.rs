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
    let mut stars: Vec<Star> = vec![];
    let mut iter = buff_reader.lines().map(|l| l.unwrap()).enumerate();
    iter.next(); // the heading 
    for (index, line) in iter{
        let line_contents = line.split(',').collect::<Vec<&str>>();
        if line_contents.len() != 3 as usize {
            return Err("The hyg star db has a line with more or less than 3 entries".to_string());
        } 

        let ra = match line_contents[0].parse::<f64>() {
            Ok(ra) => ra, 
            Err(_) => return Err(format!("rigth asc on row {} is not float", index))
        };
        let dec = match line_contents[1].parse::<f64>() {
            Ok(ra) => ra, 
            Err(_) => return Err(format!("dec on row {} is not float", index))
        };
        let brightness = match line_contents[2].parse::<f64>() {
            Ok(ra) => ra, 
            Err(_) => return Err(format!("brightness on row {} is not float", index))
        };

        stars.push(Star::new(ra*std::f64::consts::PI/180.0, dec*std::f64::consts::PI/180.0, brightness));
    }
    Ok(stars)
}

/// takes the limiting magnitude and the number of stars to be included in the flower pattern, and 
/// generates a database of DFT coefficients for the functions r(i) - the distance from the central
/// star to petel number i, and 
pub fn generate_db(m_lim: f64, k: u16) -> Result<(), String>{
    let stars = read_hyd_database()?; 
    for i in 0..10 {
        println!("{:?}", stars.get(i).unwrap())
    }
    Ok(())
}
