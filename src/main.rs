mod parse_stars;
mod coordinates;

fn main() {
    println!("Hello, world!");
    if let Err(e) = coordinates::generate_database::generate_db(1.1, 1) {
        println!("{}", e);
    }
}
