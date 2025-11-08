use std::{fs::File, path::Path};

use gmt_field_aberrations::{
    segment::{rbm::Rbm, units::Txyz},
    zernike::OpdToZernike,
    Field, Mirror, PupilMode,
};

fn main() -> color_eyre::Result<()> {
    for i in 1..=7 {
        let sid = Some(i);
        // let rbm= Some(Rbm::t_x(Txyz::Mu(1f64)));
        let azimuth = skyangle::SkyAngle::Degree(90.);
        let aberrations = (-20..21)
            .map(|i| {
                let zenith = skyangle::SkyAngle::Arcminute(0.5 * i as f32);
                let field = Field::new(5)
                    .pointing((zenith, azimuth))
                    // .pupil_mode(PupilMode::Full);
                    .pupil_mode(
                        sid.map_or_else(|| PupilMode::Full, |sid| PupilMode::m2(sid).non_zeroed()),
                        // sid.map_or_else(
                        //     || PupilMode::Full,
                        //     |sid| {
                        //         PupilMode::segment(sid, Mirror::m2(Rbm::t_y(Txyz::Mu(1e3f64))))
                        //             .non_zeroed()
                        //     },
                        // ),
                    );
                field
                    .zernike(OpdToZernike::LeastSquareFit)
                    .map(|f| (zenith.into_value(), f.coefficients().to_vec()))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join("field-aberrations")
            .join(sid.map_or_else(
                || "field-aberrations.pkl".into(),
                |sid| format!("field-aberrations_s{sid}_90_lsq.pkl"),
            ));
        println!("saving field aberrations to {path:?}");
        serde_pickle::to_writer(&mut File::create(path)?, &aberrations, Default::default())?;
    }
    Ok(())
}
