use std::fs::File;

use crseo::{Builder, raytracing::Rays};
use delrays::{Gmt, Mesh, Trace};
use plotters::prelude::*;
use triangle_rs::Delaunay;

fn main() -> anyhow::Result<()> {
    let rim_diameter = 8.365;
    let delta_rim = 1f64 / 4f64;
    let delaunay = Delaunay::disc(rim_diameter, delta_rim, Some([0., 8.71]));
    delaunay.plot();
    println!("{}", delaunay);

    let mut rays: Rays = Rays::from_mesh(&delaunay).build()?;
    dbg!(rays.chief_coordinates());

    let mut gmt = Gmt::new()?;
    rays.trace(&mut gmt);

    let opds = rays.opds();
    serde_pickle::to_writer(&mut File::create("opds.pkl")?, &opds, Default::default())?;

    let fig = SVGBackend::new("asm_facesheet.svg", (768, 768)).into_drawing_area();
    fig.fill(&WHITE).unwrap();


    Ok(())
}
