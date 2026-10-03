use std::{
    io::{Write, stdout},
    thread,
    time::Instant,
};

use crseo::{CrseoError, FromBuilder, raytracing::Rays};
use delrays::{
    Gmt, Mesh, Set, Trace,
    inversion::{DoubleZernikes, Layout},
    zernikes::{
        AsZernikes, FieldZernike, SegmentsDoubleZernikes, Zernike,
        fmt::SegmentsDoubleZernikesFormat,
    },
};
use skyangle::Conversion;

const PUPIL_MODES: [usize; 4] = [5, 6, 7, 8];

fn main() -> anyhow::Result<()> {
    // field_mesh.plot();
    // println!("{field_mesh}");

    let now = Instant::now();
    let mut results = vec![];
    for id in 1..=7 {
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
                        .reduce_into(PUPIL_MODES);
                    Ok(FieldZernike::new((zen, azi), modes))
                })
                .collect::<Result<Set<FieldZernike>, CrseoError>>()
                .unwrap();
            field_zernikes.as_zernikes(4, &field_mesh)
        }));
    }

    let segments: SegmentsDoubleZernikes =
        results.into_iter().map(|res| res.join().unwrap()).collect();
    println!("\nElapsed time: {:.3?}", now.elapsed());

    println!("{}", SegmentsDoubleZernikesFormat::from(&segments).width(9).precision(3));

    // aberration coefficients Ω_klm (nm) from pupil astigmatism and coma
    let layout = if std::env::var("ENTRANCE").is_ok() {
        Layout::entrance()
    } else {
        Layout::exit()
    };
    let omega = DoubleZernikes::from(&segments).scaled(1e9).invert(&layout);
    println!("Aberration coefficients (nm):\n{omega:.4}");

    Ok(())
}
