use std::{fs::File, io::BufWriter, path::Path, time::Instant};

use gmt_field_aberrations::{
    segment::{
        rbm::Rbm,
        units::{Rxyz, Txyz},
    },
    Mirror, Probes, PupilMode,
};
use skyangle::SkyAngle;
use triangle_rs::Builder;

fn main() -> color_eyre::Result<()> {
    env_logger::init();

    let radius = 7. * 0.5;
    let perimeter = 2. * std::f64::consts::PI * radius;
    let delta = 0.5f64;
    let n = (perimeter / delta).ceil() as usize;
    let nodes: Vec<_> = (0..n)
        .flat_map(|i| {
            let o = 2. * std::f64::consts::PI * i as f64 / n as f64;
            vec![radius * o.cos(), radius * o.sin()]
        })
        .collect();
    let tri = {
        let mut builder = Builder::new();
        builder.add_polygon(&nodes).set_switches("QDqa0.25");
        builder.build()
    };
    println!(
        "Surface area: {:.3}/{:.3}",
        std::f64::consts::PI * radius * radius,
        tri.area()
    );
    println!("{tri}");

    let field_angles: Vec<_> = tri
        .vertex_iter()
        .map(|xy| {
            let x = xy[0];
            let y = xy[1];
            let z = x.hypot(y);
            let a = y.atan2(x);
            (z, a)
        })
        .map(|(z, a)| (SkyAngle::Arcminute(z as f32), SkyAngle::Radian(a as f32)))
        .collect();
    let now = Instant::now();
    let probes = Probes::new(
        field_angles,
        4,
        PupilMode::segment(
            1,
            Mirror::M2(Rbm::t_y(Txyz::Mu(20.)) + Rbm::r_x(Rxyz::Arcsecond(1.))),
            // Mirror::Both(
            //     Rbm::r_x(Rxyz::Arcsecond(1.)),
            //     Rbm::r_x(Rxyz::Arcsecond(-8.)),
            // ),
        ),
    );
    println!("probed field in {:#?}", now.elapsed());

    let file = File::create(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join("segment-field-aberrations")
            .join("segment-field-aberrations.pkl"),
    )?;
    let mut buffer = BufWriter::new(file);
    serde_pickle::to_writer(&mut buffer, &probes, Default::default())?;
    Ok(())
}
