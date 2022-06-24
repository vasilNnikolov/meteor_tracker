use crate::take_picture;
pub mod star;

/// takes a Picture struct, and tries to parse it into a vector of Star structs, in the coordinate
/// system of the camera
/// if it deems the picture bad (too many clouds, lights etc.), it returns a string with a
/// corresponding error message
pub fn parse_image_to_stars(image: take_picture::Picture) -> Result<Vec<star::Star>, String> {
    Err(String::from("parse_image_to_stars not unimplemented"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_ra_dec() {
        let s = star::Star::new(-0.23, -1.0, 3.3, 1);
        let (ra, dec) = s.get_ra_dec();
        assert_eq!(ra, -0.23);
        assert_eq!(dec, -1.0);
    }
}
