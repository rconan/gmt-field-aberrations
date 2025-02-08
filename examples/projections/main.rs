use std::fs::File;

use gmt_field_aberrations::{
    segment::{rbm::Rbm, units::Txyz},
    Mirror, PupilMode,
};
use serde::Serialize;
use skyangle::{Conversion, SkyAngle};

#[derive(Debug, Serialize)]
pub struct Data {
    phase: Vec<f32>,
    zern: Vec<f64>,
    r: Vec<f64>,
    o: Vec<f64>,
}

fn main() -> color_eyre::Result<()> {
    env_logger::init();

    let gmt_mode = PupilMode::Full;
    // let gmt_mode = PupilMode::Segment(1);
    let rbm = Rbm {
        t_xyz: [Txyz::Mu(0.), Txyz::Mu(0.), Txyz::Mu(0.)],
        r_xyz: [
            SkyAngle::Arcsecond(1f64 / 8.),
            SkyAngle::Arcsecond(0f64),
            SkyAngle::Arcsecond(0f64),
        ],
    };
    let mirror = Mirror::M2(rbm);

    let n_radial_order = 5;

    let a = 0f32;
    println!("{:?}", gmt_mode);
    println!("Azimuth: {}deg", a);
    println!(r#"z["]{:8.0?}"#, (1..=15).collect::<Vec<_>>());
    let z = 6;
    let zernp = match gmt_mode {
        PupilMode::Full => {
            mirror::field_zernike((z as f32).from_arcmin(), a.to_radians(), n_radial_order)?
        }
        PupilMode::Segment { sid, ref mirror } => segment::field_zernike(
            sid,
            (z as f32).from_arcmin(),
            a.to_radians(),
            n_radial_order,
            &mirror,
        )?,
    };
    serde_pickle::to_writer(
        &mut File::create("examples/projections/projection.pkl")?,
        &zernp,
        Default::default(),
    )?;
    println!(
        "{:3.0} {:8.0?} {:4.0}/{:4.0}",
        z,
        zernp.coefficients(),
        zernp.opd_std(),
        zernp.std()
    );
    Ok(())
}
