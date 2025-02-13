use faer::{MatMut, MatRef};
use gmt_field_aberrations::{
    segment::{
        dof::Dof,
        rbm::Rbm,
        units::{Rxyz, Txyz},
    },
    Field, Mirror, PupilMode,
};
use skyangle::{Conversion, SkyAngle};

pub fn get_coefs(
    sid: i32,
    n_radial_order: u32,
    // rbm: impl Into<Rbm> + Clone,
    mirror: Mirror,
    azimuth: [i32; 3],
    skip: usize,
    take: usize,
    calib: &mut Vec<f64>,
    scale: f64,
) -> color_eyre::Result<()> {
    for az in azimuth {
        let pupil_mode = PupilMode::Segment {
            sid,
            mirror: mirror.clone(), // Mirror::M2(rbm.clone().into()),
            zeroed: true,
        };
        let gs = Field::new(n_radial_order)
            .pointing((SkyAngle::Arcminute(6.), SkyAngle::Degree(az as f32)))
            .pupil_mode(pupil_mode);
        // println!("{gs}");
        let zernp = gs.zernike()?;
        calib.extend(
            zernp
                .coefficients()
                .iter()
                .skip(skip)
                .take(take)
                .map(|x| scale * x),
        );
    }
    Ok(())
}

fn main() -> color_eyre::Result<()> {
    env_logger::init();

    let n_radial_order = 4;
    let sid = 7;
    let skip = 3;
    let take = 10 - skip;

    let mut calib = Vec::<f64>::new();
    let azimuth = [0, 120, 240];
    for dof in Dof::into_iter() {
        get_coefs(
            sid,
            n_radial_order,
            // dof,
            Mirror::M2(dof.into()),
            azimuth,
            skip,
            take,
            &mut calib,
            1e-9 / dof.as_f64(),
        )?;
    }
    let mat = MatRef::from_column_major_slice(&calib, take * azimuth.len(), Dof::len());
    println!("{mat:+8.6?}");
    let svd = mat.thin_svd().unwrap();
    let s: Vec<_> = svd.S().column_vector().iter().collect();
    dbg!(&s);
    // println!("{:.3?}", svd.U());
    let cond = s[0] / *s.last().unwrap();
    println!("Cond.: {cond}");

    // let rbm = Rbm::r_x(Rxyz::Arcsecond(1.));
    let mirror = Mirror::M1(Rbm::r_x(Rxyz::Arcsecond(1.))); // + Mirror::M2(Rbm::t_x(Txyz::Mu(-10.)));
    println!("{:.0?}", &mirror);
    let pupil_mode = PupilMode::Segment {
        sid,
        mirror: mirror.clone(),
        zeroed: true,
    };
    for az in azimuth {
        let zernp = Field::new(n_radial_order)
            .pointing((SkyAngle::Arcminute(6.), SkyAngle::Degree(az as f32)))
            .pupil_mode(pupil_mode.clone())
            .zernike()?;
        // let zernp = gs.zernike()?;
        println!("{zernp}");
    }

    let mut data = Vec::<f64>::new();
    get_coefs(
        sid,
        n_radial_order,
        mirror.clone(),
        azimuth,
        skip,
        take,
        &mut data,
        1e-9,
    )?;
    let n = data.len();
    let rhs = MatMut::from_column_major_slice_mut(&mut data, n, 1);
    let is: Vec<_> = s.iter().map(|x| x.recip()).collect();
    let imat = svd.V()
        * MatRef::from_column_major_slice(&is, is.len(), 1)
            .col(0)
            .as_diagonal()
        * svd.U().transpose();
    let x = imat * rhs;
    let x: Vec<_> = x.col(0).iter().collect();
    let rbm_e = Rbm::t_x(Txyz::to_mu(*x[0]))
        + Rbm::t_y(Txyz::to_mu(*x[1]))
        + Rbm::r_x(Rxyz::Arcsecond(x[2].to_arcsec()))
        + Rbm::r_y(Rxyz::Arcsecond(x[3].to_arcsec()));
    println!("{:.0?}", &rbm_e);

    let pupil_mode = PupilMode::Segment {
        sid,
        // mirror: Mirror::M2(rbm - rbm_e),
        mirror: mirror - Mirror::M2(rbm_e),
        zeroed: true,
    };
    for az in azimuth {
        let zernp = Field::new(n_radial_order)
            .pointing((SkyAngle::Arcminute(6.), SkyAngle::Degree(az as f32)))
            .pupil_mode(pupil_mode.clone())
            .zernike()?;
        // let zernp = gs.zernike()?;
        println!("{zernp}");
    }
    Ok(())
}
