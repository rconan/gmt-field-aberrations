use std::{f64, fs::File};

use crseo::{FromBuilder, raytracing::Rays};
use delrays::{Gmt, Mesh, Stats, Trace};
use skyangle::Conversion;
use triangle_rs::Delaunay;

fn main() -> anyhow::Result<()> {
    // segment ID set with environment variable SID (or set to 7 if not present)
    let delaunay = Delaunay::gmt_segment()?;
    delaunay.plot();
    println!("{}", delaunay);

    let mut rays: Rays = Rays::from_mesh(
        &delaunay,
        Rays::builder()
            .zenith(3f64.from_arcmin())
            .azimuth(f64::consts::FRAC_PI_2),
    )?;
    dbg!(rays.chief_coordinates());

    let mut gmt = Gmt::new()?;
    rays.trace(&mut gmt);

    let opds = rays.opds();
    println!("{}", Stats::from(&opds));
    serde_pickle::to_writer(&mut File::create("opds.pkl")?, &opds, Default::default())?;

    Ok(())
}
