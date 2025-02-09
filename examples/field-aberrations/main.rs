use std::{fs::File, path::Path};

use gmt_field_aberrations::Field;

fn main() -> color_eyre::Result<()> {
    let azimuth = skyangle::SkyAngle::Degree(0.);
    let aberrations = (0..21)
        .map(|i| {
            let zenith = skyangle::SkyAngle::Arcminute(0.5 * i as f32);
            let field = Field::new(5).pointing(zenith, azimuth);
            field
                .zernike()
                .map(|f| (zenith.into_value(), f.coefficients().to_vec()))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("examples")
        .join("field-aberrations")
        .join("field-aberrations.pkl");
    serde_pickle::to_writer(&mut File::create(path)?, &aberrations, Default::default())?;
    Ok(())
}
