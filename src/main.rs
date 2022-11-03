mod coordinates;
mod parse_stars;
mod take_picture;

use rustfft::{num_complex::Complex, FftPlanner};
fn main() {
    let mut buff: Vec<Complex<f64>> = vec![];
    let n: i32 = 10;
    for i in 0..n {
        buff.push(Complex {
            re: i as f64,
            im: 0.0,
        });
    }

    let mut buff2: Vec<Complex<f64>> = vec![];
    let tau: i32 = 2;
    for i in 0..n {
        buff2.push(buff[((i - tau + n) % n) as usize]);
        println!(
            "{}, {}, {}",
            i,
            buff[i as usize].norm(),
            buff2[i as usize].norm()
        );
    }

    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(n as usize);

    fft.process(&mut buff);
    fft.process(&mut buff2);
    let mut R: Vec<Complex<f64>> = vec![];
    for i in 0..n as usize {
        R.push(buff[i] * buff2[i].conj() / (buff[i] * buff2[i].conj()).norm());
    }
    let inverse_fft = planner.plan_fft_inverse(n as usize);
    inverse_fft.process(&mut R);
    for i in 0..n as usize {
        println!("{}: {}", i, R[i].norm());
    }
    let (max_index, score) = R
        .iter()
        .map(|&c| c.norm())
        .enumerate()
        .max_by(|&(_, f1), &(_, f2)| f1.partial_cmp(&f2).unwrap())
        .unwrap();
    println!("score: {}, tau is {}", score, (n - max_index as i32) % n);
}
