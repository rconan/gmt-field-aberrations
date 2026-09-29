use std::fs::File;

use crseo::{Builder, raytracing::Rays};
use delrays::{Gmt, Mesh, Set, Trace};
use triangle_rs::Delaunay;

fn main() -> anyhow::Result<()> {
    // segment ID set with environment variable SID (or set to 7 if not present)
    let delaunay = Set::<Delaunay>::gmt();
    delaunay.plot();
    // println!("{}", delaunay);

    // let mut rays: Rays = Rays::from_mesh(&delaunay).build()?;
    // dbg!(rays.chief_coordinates());

    // let mut gmt = Gmt::new()?;
    // rays.trace(&mut gmt);

    // let opds = rays.opds();
    // serde_pickle::to_writer(&mut File::create("opds.pkl")?, &opds, Default::default())?;

    Ok(())
}
