use crate::parse_stars::star;
use crate::parse_stars::star::Star;

use nalgebra::{Matrix3, Vector3};
pub struct FlowerPattern {
    /// the central star struct
    pub central_star: Star,
    /// the list of top K outer stars by brightness
    pub outer_stars: Vec<Star>,
    /// the list of the distances between the central star and each of the outer stars, in radians
    pub r: Vec<f64>,
    /// list of angles between lines connecting outer stars (delta[i] = angle between outer stars i and i+1)
    /// delta[k-1] = angle between star k-1 and star 0
    pub delta: Vec<f64>,
    /// the fov radius in radians, with which the flower pattern was generated
    pub fov: f64,
}

/// a helper function which computes the angle between the x axis of the camera (see Sky coordinate
/// system) in README, and a given star, in radians
/// It is used to order the stars around the central star when generating the flower pattern
pub fn angle_of_outer_petel(central_star: &Star, outer_star: &Star) -> f64 {
    let y_prime = central_star.coords;
    let mut x_prime = y_prime.cross(&Vector3::new(0.0, 0.0, 1.0));
    x_prime = x_prime / x_prime.norm();
    let z_prime = x_prime.cross(&y_prime);

    let r = Matrix3::from_columns(&[x_prime, y_prime, z_prime]);
    // TODO work out the inverse by hand and hard-code it, for slightly better code speed
    let r_inv = r.try_inverse().unwrap();

    let petel_new_coords = r_inv * outer_star.coords;
    let petel_angle = petel_new_coords[(2, 0)].atan2(petel_new_coords[(0, 0)]);

    // return angle in range [0, 2pi)
    if petel_angle < 0.0 {
        return petel_angle + 2.0 * std::f64::consts::PI;
    } else {
        return petel_angle;
    }
}

/// the angle between two adjacent outer stars
fn delta_angle(central_star: &Star, petel_i: &Star, petel_i_plus_1: &Star) -> f64 {
    let delta_i = angle_of_outer_petel(central_star, petel_i_plus_1)
        - angle_of_outer_petel(central_star, petel_i);
    if delta_i < 0.0 {
        return delta_i + 2.0 * std::f64::consts::PI;
    } else {
        return delta_i;
    }
}

impl FlowerPattern {
    /// generates the flower pattern for a star
    /// index - the index of the central star according to the stars vector
    /// k - number of outer stars, used to limit low brightness stars in the petels
    /// fov - RADIUS of the field of view of the camera in radians, fov = min(w, h)/2 where w and h are the
    /// sizes of the camera field in radians
    /// stars - the list of all stars as generated from the HYG db, or the list of all observed
    /// stars if the function is used on observational data
    pub fn generate(
        index: u16,
        k: u16,
        fov: f64,
        stars: &Vec<Star>,
    ) -> Result<FlowerPattern, String> {
        let mut outer_stars: Vec<Star> = vec![];
        let mut radius: Vec<f64> = vec![];
        let mut delta: Vec<f64> = vec![];

        let central_star: &Star =
            match stars.get(index as usize) {
                Some(s) => s,
                None => return Err(String::from(
                    "The index of the central star is larger than the largest one in the database",
                )),
            };

        // dot product of the coordinates should be larger than cos(FOV) in order for the star to
        // be in the FOV
        let mut stars_in_fov: Vec<&Star> = (*stars)
            .iter()
            .filter(|&s| {
                let cos_angle = central_star.coords.dot(&s.coords);
                // the second portion is so it does not include the central star itself
                cos_angle > fov.cos() && cos_angle < 0.999999
            })
            .collect();

        // order the stars, lowest magnitude first, then highest distance from central star
        stars_in_fov.sort_by(|&s1, &s2| {
            match s1.brightness.partial_cmp(&(s2.brightness)).unwrap() {
                std::cmp::Ordering::Equal => s1
                    .coords
                    .dot(&central_star.coords)
                    .partial_cmp(&(s2.coords.dot(&central_star.coords)))
                    .unwrap(),
                other => other,
            }
        });

        // take first k brightest stars
        if stars_in_fov.len() < k as usize {
            return Err(String::from(format!(
                "There are {} stars in the fov, less than k={}",
                stars_in_fov.len(),
                k
            )));
        }
        stars_in_fov = stars_in_fov.into_iter().take(k as usize).collect();

        // order the stars by their angle relative to the x axis of the sky camera, in ascending
        // order
        stars_in_fov.sort_by(|&a, &b| {
            angle_of_outer_petel(central_star, a)
                .partial_cmp(&angle_of_outer_petel(central_star, b))
                .unwrap()
        });

        for i in 0..k as usize {
            let petel = stars_in_fov[i];
            let next_petel = stars_in_fov[(i + 1) % k as usize];
            outer_stars.push(*petel);
            radius.push(star::angle_between_stars(petel, central_star));
            delta.push(delta_angle(central_star, petel, next_petel));
        }

        Ok(FlowerPattern {
            central_star: *central_star,
            outer_stars,
            r: radius,
            delta,
            fov,
        })
    }

    /// returns the angle between the petel of index petel_index in the vector of stars, and the x
    /// axis of the Sky coordinate system
    /// petel_index - between 0 and k
    /// angle returned is between 0 and 2 pi
    pub fn angle_of_petel(&self, petel_index: u16) -> Result<f64, String> {
        if petel_index > self.outer_stars.len() as u16 {
            return Err(String::from(
                "petel_index is larger than k, the number of petels in the flower pattern",
            ));
        }
        let y_prime = self.central_star.coords;
        let mut x_prime = y_prime.cross(&Vector3::new(0.0, 0.0, 1.0));
        x_prime = x_prime / x_prime.norm();

        let z_prime = x_prime.cross(&y_prime);

        let r = Matrix3::from_columns(&[x_prime, y_prime, z_prime]);
        // todo work out the inverse by hand and hard-code it, for slightly better code speed
        let r_inv = r.try_inverse().unwrap();

        let petel_new_coords = r_inv * self.outer_stars[petel_index as usize].coords;
        let petel_angle = petel_new_coords[(2, 0)].atan2(petel_new_coords[(0, 0)]);

        // return angle in range [0, 2pi)
        if petel_angle < 0.0 {
            return Ok(petel_angle + 2.0 * std::f64::consts::PI);
        } else {
            return Ok(petel_angle);
        }
    }
}
