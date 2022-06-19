use crate::parse_stars::star::Star;  
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

impl FlowerPattern {
    /// generates the flower pattern for a star
    /// index - the index of the central star according to the HYG db
    /// fov - RADIUS of the field of view of the camera, fov = min(w, h)/2 where w and h are the
    /// sizes of the camera field in radians
    /// stars - the list of all stars as generated from the HYG db
    pub fn generate(index: u16, fov: f64, stars: Vec<Star>) -> FlowerPattern {
        let mut outer_stars: Vec<Star>;
        let mut radius: Vec<f64>;
        let mut delta: Vec<f64>;
        
        FlowerPattern {
            center_star_index: index, 
            outer_stars, 
            r: radius, 
            delta, 
            fov
        }
    } 
}

/// generates the flower pattern for the 
pub fn generate_flower_pattern(stars: &Vec<Star>, index: u16) {
    
}
