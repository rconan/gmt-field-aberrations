use std::{fs::File, io::BufWriter};

use crseo::{gmt, utilities::MaskFilter, Builder, FromBuilder, Source};
use geotrans::{Segment, SegmentTrait, Transform, M1};
use gmt_field_aberrations::zernike::ZernikeBasis;
use skyangle::Conversion;

fn main() -> anyhow::Result<()> {
    env_logger::init();

    let sid = 7;
    let n_radial_order = 5;

    let za = 0f32;
    let aa = 0f32;

    let gramschmidt = false;

    println!("computing Zernike modes first {n_radial_order} radial order for segment #{sid} at {za}arcmin and {aa}deg of-axis");

    let segment = Segment::<M1>::new(sid)?;
    let [x0, y0, _] = Transform::to([0f64; 3], segment);
    log::info!("segment origin: {:.3?}", (x0, y0));

    let z = vec![za.from_arcmin()];
    let a = vec![aa.to_radians()];

    let mut gmt = gmt!();

    // segment chief rays (@ origin)
    let xyz_chief = {
        let mut src = Source::builder()
            .zenith_azimuth(z.clone(), a.clone())
            .rays_coordinates(vec![x0], vec![y0])
            .build()?;
        let rays = &mut src.as_raw_mut_ptr().rays;
        unsafe {
            gmt.m1.trace(rays);
            gmt.m2.trace(rays);
            rays.to_sphere1(-5.830, 2.197173);
        };

        src.rays().coordinates()
    };
    log::info!("segment chief ray at exit pupil: {:.3?}", &xyz_chief);

    // selecting the GMT segment
    gmt.keep(&[sid]);

    let mut src = Source::builder().zenith_azimuth(z, a).build()?;
    src.through(&mut gmt).xpupil();

    // building the zernike basis
    let xyz = src.rays().coordinates();
    let xy: Vec<_> = src
        .rays()
        .mask()
        .filter(xyz.chunks(3))
        .map(|c| [c[0] - xyz_chief[0], c[1] - xyz_chief[1]])
        .collect();
    log::info!("segment rays #: {}", xy.len());
    let zern = ZernikeBasis::builder(n_radial_order, &xy)
        .gramschmidt(gramschmidt)
        .build();

    let filename = format!(
        "zernike-modes_{n_radial_order}_gs{}_s{sid}_z{za:.1}a{aa:.1}.pkl",
        i32::from(gramschmidt)
    );
    let file = File::create(&filename)?;
    let mut buffer = BufWriter::new(file);
    serde_pickle::to_writer(&mut buffer, &zern, Default::default())?;
    println!("Zernike modes written to {filename}");
    Ok(())
}
