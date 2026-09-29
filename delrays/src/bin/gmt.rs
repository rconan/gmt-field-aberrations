use std::fs::File;

use crseo::raytracing::Rays;
use delrays::{Gmt, Mesh, Set, Stats, Trace};
use triangle_rs::Delaunay;

fn main() -> anyhow::Result<()> {
    let delaunay = Set::<Delaunay>::gmt();
    delaunay.plot();
    // println!("{}", delaunay);

    let mut rays = Set::<Rays>::from_mesh(&delaunay)?;

    let mut gmt = Gmt::new()?;
    rays.trace(&mut gmt);

    let opds = rays.opds();
    println!("{}", Stats::from(&opds));
    serde_pickle::to_writer(&mut File::create("opds.pkl")?, &opds, Default::default())?;

    Ok(())
}
