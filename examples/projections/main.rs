use std::fs::File;

use gmt_field_aberrations::{Field, PupilMode};

fn main() -> color_eyre::Result<()> {
    env_logger::init();

    // let zernp = Field::default().zernike()?;
    let zernp = Field::default().pupil_mode(PupilMode::m1(1)).zernike()?;
    serde_pickle::to_writer(
        &mut File::create("examples/projections/projection.pkl")?,
        &zernp,
        Default::default(),
    )?;
    // println!(
    //     "{:3.0} {:8.0?} {:4.0}/{:4.0}",
    //     z,
    //     zernp.coefficients(),
    //     zernp.opd_std(),
    //     zernp.std()
    // );
    Ok(())
}
