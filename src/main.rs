mod parse_stars;
mod coordinates;

use nalgebra::{Vector3, Matrix3};
fn main() {
     // println!("Hello, world!");
     // if let Err(e) = coordinates::generate_database::generate_db(1.1, 1) {
     //     println!("{}", e);
     // }
    

    let M = Matrix3::from_columns(&[
        Vector3::new(1.0, 2.0, 3.0), Vector3::new(4.0, 5.0, 6.0), Vector3::new(7.0, 8.0, 9.0)
    ]);
    let v = Vector3::new(1.0, 2.0, 3.0);
    println!("{:?}", v.shape());

    println!("{:?}", M);
    println!("{:?}", v.transpose()*M);
}
