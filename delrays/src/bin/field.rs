use std::fs::File;

use crseo::{FromBuilder, raytracing::Rays};
use delrays::{FieldZernike, Gmt, Mesh, Trace, Zernike};
use skyangle::Conversion;
use triangle_rs::Delaunay;

fn main() -> anyhow::Result<()> {
    let field_mesh = Delaunay::disc(20f64, 2., None);
    field_mesh.plot();
    println!("{field_mesh}");

    let iter = field_mesh
        .vertex_iter()
        .map(|xy| (xy[0].hypot(xy[1]), xy[1].atan2(xy[0])));

    // segment ID set with environment variable SID (or set to 7 if not present)
    let delaunay = Delaunay::gmt_segment()?;
    println!("{}", delaunay);

    let mut gmt = Gmt::new()?;

    let mut field_zernikes = vec![];
    for (zen, azi) in iter {
        let rays_builder = Rays::builder().zenith(zen.from_arcmin()).azimuth(azi);
        let mut rays: Rays = Rays::from_mesh(&delaunay, rays_builder)?;
        let modes = rays
            .trace(&mut gmt)
            .opds()
            .as_zernikes(4)
            .reduce_into(&[5, 6]);
        field_zernikes.push(FieldZernike::new((zen, azi), modes));
    }
    serde_pickle::to_writer(
        &mut File::create("zerns.pkl")?,
        &field_zernikes,
        Default::default(),
    )?;

    Ok(())
}
