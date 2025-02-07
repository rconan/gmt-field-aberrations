use std::{fs::File, path::Path};

use gmt_field_aberrations::Field;

fn main() -> color_eyre::Result<()> {
    let config = Path::new(env!("CARGO_MANIFEST_DIR")).join("default-config.ron");
    let mut file = File::create(config)?;
    let field = Field::default();
    ron::ser::to_writer_pretty(&mut file, &field, Default::default())?;
    Ok(())
}
