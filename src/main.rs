use std::{fs::File, path::Path};

use gmt_field_aberrations::Field;

use serde::Serialize;

#[derive(Debug, Serialize)]
pub struct Data {
    phase: Vec<f32>,
    zern: Vec<f64>,
    r: Vec<f64>,
    o: Vec<f64>,
}

fn main() -> color_eyre::Result<()> {
    env_logger::init();

    let file = File::open(Path::new(env!("CARGO_MANIFEST_DIR")).join("config.ron"))?;
    let field: Field = ron::de::from_reader(&file)?;
    println!("{:-<100}", "");
    println!("{field}");
    let zernp = field.zernike()?;
    println!("{zernp}");
    Ok(())
}
