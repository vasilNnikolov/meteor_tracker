pub mod generate_database;
pub mod flower;

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

}
