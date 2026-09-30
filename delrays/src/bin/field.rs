use std::fs::File;

use crseo::{CrseoError, FromBuilder, raytracing::Rays};
use delrays::{
    Gmt, Mesh, Set, Trace,
    zernikes::{AsZernikes, FieldZernike, Zernike},
};
use skyangle::Conversion;
use triangle_rs::Delaunay;

fn main() -> anyhow::Result<()> {
    let field_mesh = Delaunay::disc(20f64, 2., None);
    field_mesh.plot();
    // println!("{field_mesh}");

    let iter = field_mesh
        .vertex_iter()
        .map(|xy| (xy[0].hypot(xy[1]), xy[1].atan2(xy[0])));

    // segment ID set with environment variable SID (or set to 7 if not present)
    let delaunay = Delaunay::gmt_segment()?;
    // println!("{}", delaunay);

    let mut gmt = Gmt::new()?;

    let field_zernikes = iter
        .map(|(zen, azi)| {
            let rays_builder = Rays::builder().zenith(zen.from_arcmin()).azimuth(azi);
            let mut rays: Rays = Rays::from_mesh(&delaunay, rays_builder)?;
            let modes = rays
                .trace(&mut gmt)
                .opds()
                .as_zernikes(4)
                .reduce_into([5, 6]);
            Ok(FieldZernike::new((zen, azi), modes))
        })
        .collect::<Result<Set<FieldZernike>, CrseoError>>()?;

    serde_pickle::to_writer(
        &mut File::create("zerns.pkl")?,
        &field_zernikes,
        Default::default(),
    )?;

    let field_zern_coefs = field_zernikes.as_zernikes(4);
    for fzc in field_zern_coefs.into_iter() {
        println!("{:?}", fzc.jnm);
        let coef = fzc.coef.reduce_into(1..=8);
        coef.jnm()
            .iter()
            .zip(coef.coefficients().iter())
            .for_each(|(jnm, c)| println!(" {jnm:2?}: {:+6.0}", c * 1e9));
    }

    Ok(())
}
