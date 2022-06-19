
#[derive(Debug, Copy, Clone)]
pub struct Star {
    /// the index of the star, serving as its name. Stars are to be refered to by index
    pub index: u16,
    /// the x, y, z coordinates of a star in the earth reference frame
    pub coords: [f64; 3], 
    /// the visual brigthness of the star in magnitudes
    pub brightness: f64,
}

/// returns the cosine of the angle between two stars
pub fn cos_angle_between_stars(s1: &Star, s2: &Star) -> f64 {
    // dot product
    s1.coords.iter().zip(s2.coords.iter()).map(|(&a, &b)| a*b).sum()
}

pub fn angle_between_stars(s1: &Star, s2: &Star) -> f64 {
    cos_angle_between_stars(s1, s2).acos()
}

impl Star {
    pub fn new(ra: f64, dec: f64, brightness: f64, index: u16) -> Star {
        Star {
            index,
            coords: [dec.cos()*ra.cos(), dec.cos()*ra.sin(), dec.sin()],
            brightness
        }     
    }
}


