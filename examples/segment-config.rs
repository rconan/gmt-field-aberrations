use std::{fs::File, path::Path};

use gmt_field_aberrations::Field;

fn main() -> color_eyre::Result<()> {
    let config = Path::new(env!("CARGO_MANIFEST_DIR")).join("segment-config.ron");
    let mut file = File::create(config)?;
    let field = Field::default().pupil_mode(gmt_field_aberrations::PupilMode::Segment {
        sid: 1,
        mirror: gmt_field_aberrations::Mirror::M1(Default::default()),
        zeroed: true,
    });

    ron::ser::to_writer_pretty(&mut file, &field, Default::default())?;
    Ok(())
}
