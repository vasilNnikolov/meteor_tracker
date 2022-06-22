pub mod star;

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
