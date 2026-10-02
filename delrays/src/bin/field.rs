//! # Computes Zernike coefficients accross a 20' field-of-view
//!
use std::{fs::File, time::Instant};

use crseo::{CrseoError, FromBuilder, raytracing::Rays};
use delrays::{
    Gmt, Mesh, Set, Trace,
    zernikes::{AsZernikes, FieldZernike, Zernike},
};
use skyangle::Conversion;

fn main() -> anyhow::Result<()> {
    let field_mesh = Mesh::disc(20f64, 2., None);
    field_mesh.plot();
    println!("{field_mesh}");
    let iter = field_mesh
        .vertex_iter()
        .map(|xy| (xy[0].hypot(xy[1]), xy[1].atan2(xy[0])));

    // segment ID set with environment variable SID (or set to 7 if not present)
    let delaunay = Mesh::gmt_segment()?;
    // println!("{}", delaunay);

    let mut gmt = Gmt::new()?;

    let now = Instant::now();
    let field_zernikes = iter
        .map(|(zen, azi)| {
            let rays_builder = Rays::builder().zenith(zen.from_arcmin()).azimuth(azi);
            let mut rays: Rays = Rays::from_mesh(&delaunay, rays_builder)?;
            let opds = rays.trace(&mut gmt).opds_centered(&delaunay).unwrap();
            let modes = opds.as_zernikes(4, &delaunay).reduce_into([5, 6]);
            Ok(FieldZernike::new((zen, azi), modes))
        })
        .collect::<Result<Set<FieldZernike>, CrseoError>>()?;
    println!("Elapsed time: {:.3?}", now.elapsed());

    serde_pickle::to_writer(
        &mut File::create("zerns.pkl")?,
        &field_zernikes,
        Default::default(),
    )?;

    let field_zern_coefs = field_zernikes.as_zernikes(4, &field_mesh);
    dbg!(field_zern_coefs.mean(&field_mesh));
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
