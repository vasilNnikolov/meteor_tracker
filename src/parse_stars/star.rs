
#[derive(Debug)]
pub struct Star {
    /// the x, y, z coordinates of a star in the earth reference frame
    pub coords: [f64; 3], 
    /// the visual brigthness of the star in magnitudes
    pub brightness: f64,
}

impl Star {
    pub fn new(ra: f64, dec: f64, brightness: f64) -> Star {
        Star {
            coords: [dec.cos()*ra.cos(), dec.cos()*ra.sin(), dec.sin()],
            brightness
        }     
    }
}


