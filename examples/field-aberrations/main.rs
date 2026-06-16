use std::{fs::File, path::Path};

use gmt_field_aberrations::{
    segment::{
        rbm::Rbm,
        units::{Rxyz, Txyz},
    },
    zernike::OpdToZernike,
    Field, Mirror, PupilMode,
};
use skyangle::Conversion;

const OPD2ZERN: OpdToZernike = OpdToZernike::Projection;

fn main() -> color_eyre::Result<()> {
    let segments_rbm =
        geotrans::Mirror::<geotrans::M2>::tiptilt_2_rigidbodymotions((1f64.from_arcsec(), 0f64));
    let mut rbms = segments_rbm.chunks(6);

    for sid in 1..=7 {
        // let sid = Some(i);
        // let rbm= Some(Rbm::t_x(Txyz::Mu(1f64)));
        let rbm = rbms.next().unwrap();
        println!(" {}:, {:?}", sid, rbm);
        let pm = PupilMode::segment(sid, Mirror::m2(rbm));
        // .non_zeroed();c
        let azimuth = skyangle::SkyAngle::Degree(0.);
        let aberrations = (-20..21)
            .map(|i| {
                let zenith = skyangle::SkyAngle::Arcminute(0.5 * i as f32);
                let field = Field::new(5)
                    .pointing((zenith, azimuth))
                    // .pupil_mode(PupilMode::Full);
                    .pupil_mode(
                        // sid.map_or_else(|| PupilMode::Full, |sid| PupilMode::m2(sid).non_zeroed()),
                        // sid.map_or_else(
                        //     || PupilMode::Full,
                        //     |sid| {
                        //         PupilMode::segment(sid, Mirror::m2(Rbm::r_y(Rxyz::Arcsecond(1f64))))
                        //             .non_zeroed()
                        //     },
                        // ),
                        pm.clone(),
                    );
                field
                    .zernike(OPD2ZERN)
                    .map(|f| (zenith.into_value(), f.coefficients().to_vec()))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join("field-aberrations")
            .join(
                /* sid.map_or_else(
                || "field-aberrations.pkl".into(),
                |sid| { */
                format!(
                    "field-aberrations_{pm}_{}.pkl",
                    match OPD2ZERN {
                        OpdToZernike::Projection => "prj",
                        OpdToZernike::LeastSquareFit => "lsq",
                    }
                ),
            );
        //     },
        // ));
        println!("saving field aberrations to {path:?}");
        serde_pickle::to_writer(&mut File::create(path)?, &aberrations, Default::default())?;
    }
    Ok(())
}
