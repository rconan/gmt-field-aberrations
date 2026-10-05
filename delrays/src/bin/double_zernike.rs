use std::{
    io::{Write, stdout},
    thread,
    time::Instant,
};

use crseo::{CrseoError, FromBuilder, raytracing::Rays};
use delrays::{
    Gmt, Mesh, Set, Trace,
    inversion::{
        DoubleZernikes, Layout,
        lsq::{Coefficients, Design, terms},
    },
    zernikes::{
        AsZernikes, FieldZernike, SegmentsDoubleZernikes, Zernike,
        fmt::SegmentsDoubleZernikesFormat,
    },
};
use skyangle::Conversion;

/// Pupil Noll modes from `DZ_PUPIL_MODES`, e.g. "5-8" (default) or "1-10" or "4,5,6,7,8"
fn pupil_modes() -> Vec<usize> {
    let spec = std::env::var("DZ_PUPIL_MODES").unwrap_or_else(|_| "5-8".into());
    spec.split(',')
        .flat_map(|item| match item.split_once('-') {
            Some((a, b)) => (a.trim().parse::<usize>().expect("DZ_PUPIL_MODES")
                ..=b.trim().parse::<usize>().expect("DZ_PUPIL_MODES"))
                .collect::<Vec<_>>(),
            None => vec![item.trim().parse::<usize>().expect("DZ_PUPIL_MODES")],
        })
        .collect()
}
/// Radial order of the field fit from `DZ_FIELD_ORDER` (default: 3)
fn field_order() -> usize {
    std::env::var("DZ_FIELD_ORDER")
        .map(|s| s.parse().expect("DZ_FIELD_ORDER"))
        .unwrap_or(3)
}

fn main() -> anyhow::Result<()> {
    // field_mesh.plot();
    // println!("{field_mesh}");

    let pupil = pupil_modes();
    let field_order = field_order();
    let now = Instant::now();
    let mut results = vec![];
    for id in 1..=7 {
        let pupil = pupil.clone();
        results.push(thread::spawn(move || {
            let mut lock = stdout().lock();
            write!(lock, "{id}").unwrap();
            lock.flush().unwrap();

            let field_mesh = Mesh::disc(20f64, 2., None);
            let delaunay = Mesh::gmt_segment_with_id(id).unwrap();

            let mut gmt = Gmt::new().unwrap();

            let iter = field_mesh
                .vertex_iter()
                .map(|xy| (xy[0].hypot(xy[1]), xy[1].atan2(xy[0])));

            let field_zernikes = iter
                .map(|(zen, azi)| {
                    let rays_builder = Rays::builder().zenith(zen.from_arcmin()).azimuth(azi);
                    let mut rays: Rays = Rays::from_mesh(&delaunay, rays_builder)?;
                    let rays = rays.trace(&mut gmt);
                    // ENTRANCE=1 fits the segment OPD at the ray launch coordinates
                    // (entrance pupil) instead of the exit-pupil sphere intercepts
                    let opds = if std::env::var("ENTRANCE").is_ok() {
                        rays.opds_entrance(&delaunay)
                    } else {
                        rays.opds_centered(&delaunay).unwrap()
                    };
                    let modes = opds
                        .as_zernikes(4, &delaunay)
                        .reduce_into(pupil.iter().copied());
                    Ok(FieldZernike::new((zen, azi), modes))
                })
                .collect::<Result<Set<FieldZernike>, CrseoError>>()
                .unwrap();
            field_zernikes.as_zernikes(field_order + 1, &field_mesh)
        }));
    }

    let segments: SegmentsDoubleZernikes =
        results.into_iter().map(|res| res.join().unwrap()).collect();
    println!("\nElapsed time: {:.3?}", now.elapsed());

    let n_field = (field_order + 1) * (field_order + 2) / 2;
    println!(
        "{}",
        SegmentsDoubleZernikesFormat::from(&segments)
            .width(9)
            .precision(3)
            .field_modes(n_field)
    );

    // aberration coefficients Ω_klm (nm) from pupil astigmatism and coma
    let layout = if std::env::var("ENTRANCE").is_ok() {
        Layout::entrance()
    } else {
        Layout::exit()
    };
    if [5, 6, 7, 8].iter().all(|j| pupil.contains(j)) && field_order >= 3 {
        let omega = DoubleZernikes::from(&segments).scaled(1e9).invert(&layout);
        println!("Aberration coefficients, closed form (nm):\n{omega:.4}");
    }

    // least squares, with all the pupil and field modes, when the field fit reaches
    // radial order 5 (model k<=5, l<=8 if segment defocus or trefoil is available, else l<=6)
    if field_order >= 5 {
        let data = Coefficients::from(&segments).scaled(1e9);
        let lmax = if [4, 9, 10].iter().any(|j| pupil.contains(j)) {
            8
        } else {
            6
        };
        let design = Design::new(terms(5, lmax, 3), data.pupil(), data.field(), &layout);
        let solution = design.solve(&data)?;
        println!("Aberration coefficients, least squares (nm):\n{solution:.4}");
    }

    Ok(())
}
