use crate::parse_stars::star::Star;  
pub struct Flower {
    /// the index of the center star, according to the hyg database
    /// indeces start at 1, and bigger index means less bright star
    pub center_star_index: u16, 
    /// the list of top K outer stars by brightness
    pub outer_stars: Vec<Star>, 
    /// the list of the distances between the central star and each of the outer stars, in radians
    pub r: Vec<f64>,
    /// list of angles between lines connecting outer stars (delta[i] = angle between outer stars i and i+1)
    /// delta[k-1] = angle between star k-1 and star 0
    pub delta: Vec<f64>
}



/// generates the flower pattern for the 
pub fn generate_flower_pattern(stars: &Vec<Star>, index: u16) {
    
}
