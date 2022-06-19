use crate::parse_stars::star::Star;  
use crate::parse_stars::star;
pub struct FlowerPattern {
    /// the index of the center star, according to the hyg database
    /// indeces start at 1, and bigger index means less bright star
    pub center_star_index: u16, 
    /// the list of top K outer stars by brightness
    pub outer_stars: Vec<Star>, 
    /// the list of the distances between the central star and each of the outer stars, in radians
    pub r: Vec<f64>,
    /// list of angles between lines connecting outer stars (delta[i] = angle between outer stars i and i+1)
    /// delta[k-1] = angle between star k-1 and star 0
    pub delta: Vec<f64>, 
    /// the fov radius in radians, with which the flower pattern was generated
    pub fov: f64
}

/// a helper function which computes the angle between the x axis of the camera (see Sky coordinate
/// system) in README, and a given star, in radians
/// It is used to order the stars around the central star when generating the flower pattern
fn angle_of_outer_petel(central_star: &Star, outer_star: &Star) -> f64 {
    // TODO implement delta funcion
    1.1
}

impl FlowerPattern {
    /// generates the flower pattern for a star
    /// index - the index of the central star according to the HYG db
    /// k - number of outer stars, used to limit low brightness stars in the petels
    /// fov - RADIUS of the field of view of the camera, fov = min(w, h)/2 where w and h are the
    /// sizes of the camera field in radians
    /// stars - the list of all stars as generated from the HYG db
    pub fn generate(index: u16, k: u16, fov: f64, stars: Vec<Star>) -> Result<FlowerPattern, String> {
        let mut outer_stars: Vec<Star> = vec![];
        let mut radius: Vec<f64> = vec![];
        let mut delta: Vec<f64> = vec![];

        let central_star: &Star = match stars.get(index as usize - 1) {
            Some(s) => s, 
            None => return Err(String::from("The index of the central star is larger than the largest one in the database or is 0"))
        };

        // dot product of the coordinates should be larger than cos(FOV) in order for the star to
        // be in the FOV
        let mut stars_in_fov: Vec<&Star> = stars.iter().filter(|&s| {
            star::cos_angle_between_stars(central_star, s) > fov.cos()
        }).collect();
        // take first k brightest stars
        if stars_in_fov.len() < k as usize {
            return Err(String::from(format!("There are {} stars in the fov, more than k={}", stars_in_fov.len(), k)));
        }
        let mut stars_in_fov: Vec<&Star> = stars_in_fov.into_iter().take(k as usize).collect();
        
        // order the stars by their angle relative to the x axis of the sky camera, in ascending
        // order
        stars_in_fov.sort_by(|&a, &b| {
            angle_of_outer_petel(central_star, a).partial_cmp(&angle_of_outer_petel(central_star, b)).unwrap() 
        });

        Ok(FlowerPattern {
            center_star_index: index, 
            outer_stars, 
            r: radius, 
            delta, 
            fov
        })
    }
}

