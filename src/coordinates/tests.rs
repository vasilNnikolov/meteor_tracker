use super::*;
extern crate rand;
use rand::Rng;
fn print_flower_pattern(fp: &FlowerPattern) {
    let n = fp.outer_stars.len();
    let (ra, dec) = fp.central_star.get_ra_dec();
    // println!("central star ra, dec: {:.2} h {:.2} deg, magnitude {:.2}",
    //     ra*57.3/15.0,
    //     dec*57.3,
    //     fp.central_star.brightness);
    for i in 0..n {
        let (ra, dec) = fp.outer_stars[i].get_ra_dec();
        // println!("ra, dec: {:.2} h {:.2} deg, magnitude {:.2}, distance {:.2}, angle {:.2}",
        //     ra*57.3/15.0,
        //     dec*57.3,
        //     fp.outer_stars[i].brightness,
        //     fp.r[i]*57.295,
        //     flower::angle_of_outer_petel(&fp.central_star, &fp.outer_stars[i])*57.295);
        // println!("angle between star {} and {}: {}", i, (i + 1)%n, fp.delta[i]*57.295);
    }
}
#[test]
fn test_file_read() {
    let all_stars = generate_database::read_hyd_database(4.0, 0.0017).unwrap();
    for i in 0..5 {
        println!("{:?}", all_stars[i]);
    }
    let (ra, dec) = all_stars[0].get_ra_dec();
    println!(
        "sirius ra {}, dec {}, brightness {}",
        ra * 57.295,
        dec * 57.295,
        all_stars[0].brightness
    );
}
#[test]
fn test_sirius_pattern() {
    let all_stars = generate_database::read_hyd_database(4.0, 0.0017).unwrap();
    let sirius_pattern = FlowerPattern::generate(1, 5, 0.35, &all_stars).unwrap();
    let n = sirius_pattern.outer_stars.len();
    for i in 0..n {
        let (ra, dec) = sirius_pattern.outer_stars[i].get_ra_dec();
        println!(
            "ra, dec: {} {}, magnitude {}, distance {}, angle {}",
            ra * 57.3,
            dec * 57.3,
            sirius_pattern.outer_stars[i].brightness,
            sirius_pattern.r[i] * 57.295,
            flower::angle_of_outer_petel(&all_stars[0], &sirius_pattern.outer_stars[i]) * 57.295
        );
        println!(
            "angle between star {} and {}: {}",
            i,
            (i + 1) % n,
            sirius_pattern.delta[i] * 57.295
        );
    }
}
#[test]
fn test_canopus_pattern() {
    let all_stars = generate_database::read_hyd_database(4.0, 0.0017).unwrap();
    let canopus_pattern = FlowerPattern::generate(2, 5, 0.35, &all_stars).unwrap();
    let n = canopus_pattern.outer_stars.len();
    for i in 0..n {
        let (ra, dec) = canopus_pattern.outer_stars[i].get_ra_dec();
        println!(
            "ra, dec: {} {}, magnitude {}, distance {}, angle {}",
            ra * 57.3,
            dec * 57.3,
            canopus_pattern.outer_stars[i].brightness,
            canopus_pattern.r[i] * 57.295,
            flower::angle_of_outer_petel(&all_stars[1], &canopus_pattern.outer_stars[i]) * 57.295
        );
        println!(
            "angle between star {} and {}: {}",
            i,
            (i + 1) % n,
            canopus_pattern.delta[i] * 57.295
        );
    }
}
#[test]
fn test_polaris_pattern() {
    let all_stars = generate_database::read_hyd_database(4.0, 0.0017).unwrap();
    let canopus_pattern = FlowerPattern::generate(46, 5, 0.35, &all_stars).unwrap();
    let (ra, dec) = canopus_pattern.central_star.get_ra_dec();
    println!(
        "ra, dec of polaris: {:.2}h, {:.2}deg",
        ra * 57.29578 / 15.0,
        dec * 57.29577
    );
    let n = canopus_pattern.outer_stars.len();
    for i in 0..n {
        let (ra, dec) = canopus_pattern.outer_stars[i].get_ra_dec();
        println!(
            "ra, dec: {:.2} {:.2}, magnitude {:.2}, distance {:.2}, angle {:.2}",
            ra * 57.3,
            dec * 57.3,
            canopus_pattern.outer_stars[i].brightness,
            canopus_pattern.r[i] * 57.295,
            flower::angle_of_outer_petel(&all_stars[45], &canopus_pattern.outer_stars[i]) * 57.295
        );
        println!(
            "angle between star {} and {}: {:.2}",
            i,
            (i + 1) % n,
            canopus_pattern.delta[i] * 57.295
        );
    }
}
#[test]
fn test_db_generation() {
    if let Err(s) = generate_database::generate_db(4.0, 5, 0.52) {
        println!("{}", s);
        panic!();
    }
}
fn generate_random_orthogonal_matrix() -> Matrix3<f64> {
    let mut rng = rand::thread_rng();
    let mut x = Vector3::new(
        rng.gen_range(-1.0..1.0),
        rng.gen_range(-1.0..1.0),
        rng.gen_range(-1.0..1.0),
    );
    x = x / x.norm();
    let mut y = Vector3::new(0.0, rng.gen_range(-1.0..1.0), rng.gen_range(-1.0..1.0)).cross(&x);
    y = y / y.norm();
    let z = x.cross(&y);
    Matrix3::from_columns(&[x, y, z])
}
#[test]
fn test_main_fun_of_module() {
    // camera is facing in the direction R_y
    let R = generate_random_orthogonal_matrix();

    // generate list of captured stars
    let R_inv = R.try_inverse().unwrap();
    let fov = 0.52;
    let k = 10;
    let (DFT_db, flower_patterns) = generate_database::generate_db(4.0, k, fov).unwrap();
    let all_stars: Vec<star::Star> = flower_patterns.iter().map(|fp| fp.central_star).collect();

    let R_y = R.column(1);
    let captured_stars: Vec<star::Star> = all_stars
        .iter()
        .filter(|&s| s.coords.dot(&R_y) > (2.0 * fov).cos())
        .map(|&s| star::Star {
            index: 6969,
            coords: R_inv * s.coords,
            brightness: s.brightness,
        })
        .collect();
    println!("captured {} stars in fov", captured_stars.len());

    let experiment_R = get_rotation_matrix(&captured_stars, &DFT_db, &flower_patterns).unwrap();
    println!("exp: {}", experiment_R);
    println!("real: {}", R);
}
#[test]
fn normal_long_test() {
    let fov = 0.52;
    let k = 10;
    let (DFT_db, flower_patterns) = generate_database::generate_db(4.0, k, fov).unwrap();
    let all_stars: Vec<star::Star> = flower_patterns.iter().map(|fp| fp.central_star).collect();
    for i in 0..100 {
        let R = generate_random_orthogonal_matrix();
        let R_inv = R.try_inverse().unwrap();

        let R_y = R.column(1);
        let captured_stars: Vec<star::Star> = all_stars
            .iter()
            .filter(|&s| s.coords.dot(&R_y) > (1.5 * fov).cos())
            .map(|&s| star::Star {
                index: 6969,
                coords: R_inv * s.coords,
                brightness: s.brightness,
            })
            .collect();

        let experiment_R = get_rotation_matrix(&captured_stars, &DFT_db, &flower_patterns).unwrap();
        if i % 20 == 0 {
            println!("real: {}", R);
            println!("expe: {}", experiment_R);
        }
    }
}
// returns the average distance between the corresponding coordinate axex of the coordinate systems
fn get_score_of_guess(actual_matrix: Matrix3<f64>, experimental_matrix: Matrix3<f64>) -> f64 {
    let mut final_distance = 0.0;
    for col_index in 0..3 {
        final_distance += actual_matrix
            .column(col_index)
            .dot(&experimental_matrix.column(col_index));
    }
    final_distance / 3.0
}
#[test]
fn randomized_long_test() {
    let mut rng = rand::thread_rng();
    let rand_coefficient = 0.000;
    fn get_random_observed_coords(
        acutal_coords: Vector3<f64>,
        R_inv: Matrix3<f64>,
        rng: &mut impl Rng,
        rand_coefficient: f64,
    ) -> Vector3<f64> {
        let theoretical_coords = R_inv * acutal_coords;
        let dx = rng.gen_range(-1.0..1.0) * rand_coefficient;
        let dz = rng.gen_range(-1.0..1.0) * rand_coefficient;
        let result = theoretical_coords + Vector3::new(dx, 0.0, dz);
        result / result.norm()
    }
    let fov = 0.52;
    let k = 10;
    let (DFT_db, flower_patterns) = generate_database::generate_db(4.0, k, fov).unwrap();
    let all_stars: Vec<star::Star> = flower_patterns.iter().map(|fp| fp.central_star).collect();
    let mut good_tests = 0;
    let N = 500;
    for _ in 0..N {
        let R = generate_random_orthogonal_matrix();
        let R_inv = R.try_inverse().unwrap();

        let R_y = R.column(1);
        let captured_stars: Vec<star::Star> = all_stars
            .iter()
            .filter(|&s| s.coords.dot(&R_y) > (1.5 * fov).cos())
            .map(|&s| star::Star {
                index: 6969,
                coords: get_random_observed_coords(s.coords, R_inv, &mut rng, rand_coefficient),
                brightness: s.brightness,
            })
            .collect();

        let experiment_R = get_rotation_matrix(&captured_stars, &DFT_db, &flower_patterns).unwrap();
        let score = get_score_of_guess(R, experiment_R);
        if score > 0.90 {
            println!("Ok");
            good_tests += 1;
        } else {
            println!("score of guess: {:.5}", score);
            let observed_fp =
                super::generate_flower_pattern_from_observation(&captured_stars, k, fov).unwrap();
            let min_angle_between_stars = observed_fp
                .delta
                .iter()
                .min_by(|&&f1, &&f2| f1.partial_cmp(&f2).unwrap())
                .unwrap();
            println!(
                "min angle between stars: {:.2}",
                min_angle_between_stars * 57.2955
            );
            print_flower_pattern(&observed_fp);
        }
    }
    println!("coverage: {:.2}", 100.0 * good_tests as f64 / N as f64);
}
#[test]
fn min_distance_between_stars() {
    let all_stars = generate_database::read_hyd_database(4.0, 0.0017).unwrap();
    let mut max_cos_angle: f64 = -1.0;
    for i in 0..all_stars.len() {
        for j in 0..all_stars.len() {
            let angle = star::cos_angle_between_stars(&all_stars[i], &all_stars[j]);
            if max_cos_angle < angle && i != j {
                max_cos_angle = angle;
            }
        }
    }
    println!("max_cos_angle: {}", max_cos_angle);
}
