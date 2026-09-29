use std::fs::File;

use crseo::raytracing::Rays;
use delrays::{Gmt, Mesh, Stats, Trace};
use triangle_rs::Delaunay;

fn main() -> anyhow::Result<()> {
    // segment ID set with environment variable SID (or set to 7 if not present)
    let delaunay = Delaunay::gmt_segment()?;
    delaunay.plot();
    println!("{}", delaunay);

    let mut rays: Rays = Rays::from_mesh(&delaunay)?;
    dbg!(rays.chief_coordinates());

    let mut gmt = Gmt::new()?;
    rays.trace(&mut gmt);

    let opds = rays.opds();
    println!("{}", Stats::from(&opds));
    serde_pickle::to_writer(&mut File::create("opds.pkl")?, &opds, Default::default())?;

    Ok(())
}
