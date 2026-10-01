use std::fs::File;

use crseo::raytracing::Rays;
use delrays::{
    Gmt, Mesh, Trace,
    opd::Stats,
    zernikes::{AsZernikes, Zernike},
};
use skyangle::Conversion;

fn main() -> anyhow::Result<()> {
    // segment ID set with environment variable SID (or set to 7 if not present)
    let delaunay = Mesh::gmt_segment()?;
    delaunay.plot();
    println!("{}", delaunay);

    let mut rays: Rays = Rays::from_mesh(&delaunay, Default::default())?;
    dbg!(rays.chief_coordinates());

    let mut gmt = Gmt::new()?;
    rays.trace(&mut gmt);

    // let opds = rays.opds();
    let opds = rays.opds_centered(&delaunay).unwrap();
    println!("{}", Stats::from(&opds));
    serde_pickle::to_writer(&mut File::create("opds.pkl")?, &opds, Default::default())?;

    let zs = opds.as_zernikes(4);
    println!("{:+6.0?}", zs.speye());
    dbg!(zs.coefficients());

    Ok(())
}
