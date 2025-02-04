use std::{fs::File, io::Write, path::Path};

use gmt_field_aberrations::Field;

fn main() -> color_eyre::Result<()> {
    let config = Path::new(env!("CARGO_MANIFEST_DIR")).join("default-config.rson");
    let mut file = File::create(config)?;
    let field = Field::default();
    let rson = rson_rs::ser::pretty::to_string(&field)?;
    println!("{rson}");
    write!(file, "{}", rson)?;
    Ok(())
}
