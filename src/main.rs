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

    let file = File::open(Path::new(env!("CARGO_MANIFEST_DIR")).join("config.rson"))?;
    let field: Field = rson_rs::de::from_reader(&file)?;
    println!("{:-<100}", "");
    println!("{field}");
    let zernp = field.zernike()?;
    println!("{zernp}");
    /* // let gmt_mode = PupilMode::Full;
    let rbm = Rbm {
        t_xyz: [Txyz::Mu(0.), Txyz::Mu(0.), Txyz::Mu(0.)],
        r_xyz: [
            SkyAngle::Arcsecond(1f64 / 8.),
            SkyAngle::Arcsecond(0f64),
            SkyAngle::Arcsecond(0f64),
        ],
    };
    let mirror = Mirror::M1(rbm);
    let gmt_mode = PupilMode::Segment { sid: 1, mirror };

    let n_radial_order = 5;

    // let a: Option<f32> = Some(0f32);
    println!("{:?}", gmt_mode);
    // println!("Azimuth: {}deg", a.unwrap_or_default());
    let z = 6.;
    println!(r#"Zenith: {}'""#, z);
    println!(r#"[o]z{:8.0?}"#, (1..=15).collect::<Vec<_>>());
    println!("{:-^155}", "");
    for a in [0f32, 120., 240.] {
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
        println!(
            "{:3.0} {:8.0?} {:4.0}/{:4.0}",
            a,
            zernp.coefficients(),
            zernp.opd_std(),
            zernp.std()
        );
    } */
    Ok(())
}
