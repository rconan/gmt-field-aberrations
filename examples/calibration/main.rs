use std::{fs::File, io::BufWriter, path::Path, vec};

use gmt_field_aberrations::{
    segment::{dof::Dof, rbm::Rbm},
    Mirror, Pointing, Probes, PupilMode,
};
use skyangle::SkyAngle;

fn calibrate(
    dof: vec::IntoIter<Dof>,
    pointing: Pointing,
    sid: i32,
    mirror: Mirror,
    n_radial_order: u32,
    skip_n_mode: usize,
    rxy_to_a23: Option<&[f64]>,
) -> Vec<f64> {
    dof.flat_map(|dof| {
        let p_plus = Probes::new(
            vec![pointing],
            n_radial_order,
            PupilMode::segment(
                sid,
                match mirror {
                    Mirror::M1(_) => Mirror::m1(dof),
                    Mirror::M2(_) => Mirror::m2(dof),
                    _ => unimplemented!(),
                },
            ),
        );
        let p_minus = Probes::new(
            vec![pointing],
            n_radial_order,
            PupilMode::segment(
                sid,
                match mirror {
                    Mirror::M1(_) => Mirror::m1(-dof),
                    Mirror::M2(_) => Mirror::m2(-dof),
                    _ => unimplemented!(),
                },
            ),
        );
        let (p_plus, p_minus) = if let Some(alphas) = rxy_to_a23 {
            let r_xy: Vec<_> = p_plus.projections[0].coefficients()[1..3]
                .iter()
                .zip(alphas)
                .map(|(a, c)| a / c)
                .collect();
            let r_xy = Rbm::r_x(SkyAngle::Radian(r_xy[0])) + Rbm::r_y(SkyAngle::Radian(r_xy[1]));
            let p_plus = Probes::new(
                vec![pointing],
                n_radial_order,
                PupilMode::segment(
                    sid,
                    match mirror {
                        Mirror::M1(_) => Mirror::both(dof, -r_xy.clone()),
                        Mirror::M2(_) => Mirror::m2(Rbm::from(dof) - r_xy.clone()),
                        _ => unimplemented!(),
                    },
                ),
            );
            let r_xy: Vec<_> = p_minus.projections[0].coefficients()[1..3]
                .iter()
                .zip(alphas)
                .map(|(a, c)| a / c)
                .collect();
            let r_xy = Rbm::r_x(SkyAngle::Radian(r_xy[0])) + Rbm::r_y(SkyAngle::Radian(r_xy[1]));
            let p_minus = Probes::new(
                vec![pointing],
                n_radial_order,
                PupilMode::segment(
                    sid,
                    match mirror {
                        Mirror::M1(_) => Mirror::both(-dof, -r_xy.clone()),
                        Mirror::M2(_) => Mirror::m2(Rbm::from(-dof) - r_xy.clone()),
                        _ => unimplemented!(),
                    },
                ),
            );
            (p_plus, p_minus)
        } else {
            (p_plus, p_minus)
        };

        p_plus
            .projections
            .iter()
            .zip(p_minus.projections.iter())
            .flat_map(|(p, m)| {
                p.coefficients()
                    .iter()
                    .zip(m.coefficients())
                    .skip(skip_n_mode)
                    .map(|(x, y)| 0.5 * (x - y) / dof.as_f64())
            })
            .collect::<Vec<_>>()
    })
    .collect()
}
fn main() -> color_eyre::Result<()> {
    env_logger::init();

    let sid = 1;
    let n_radial_order = 2;
    let pointing = Pointing::new(SkyAngle::Arcminute(0.), SkyAngle::Degree(0.));
    // let m1_data = calibrate(pointing, sid, Mirror::m1(Rbm::default()), n_radial_order);
    let m2_data = calibrate(
        Dof::r_xy(),
        pointing,
        sid,
        Mirror::m2(Rbm::default()),
        n_radial_order,
        1,
        None,
    );
    dbg!(&m2_data);
    let n_radial_order = 4;
    let pointing = Pointing::new(SkyAngle::Arcminute(6.), SkyAngle::Degree(0.));
    // let m1_data = calibrate(pointing, sid, Mirror::m1(Rbm::default()), n_radial_order);
    let m1_data = calibrate(
        Dof::into_iter(),
        pointing,
        sid,
        Mirror::m1(Rbm::default()),
        n_radial_order,
        1,
        Some(&m2_data[1..3]),
    );
    let m2_data = calibrate(
        Dof::into_iter(),
        pointing,
        sid,
        Mirror::m2(Rbm::default()),
        n_radial_order,
        1,
        Some(&m2_data[1..3]),
    );

    let file = File::create(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join("calibration")
            .join("m12_closed-loop_calibration.pkl"),
    )?;
    let mut buffer = BufWriter::new(file);
    serde_pickle::to_writer(&mut buffer, &(m1_data, m2_data), Default::default())?;

    Ok(())
}
