use std::{fs::File, io::BufWriter, path::Path, time::Instant};

use clap::{Parser, Subcommand};
use gmt_field_aberrations::{
    Mirror, Pointing, Probes, PupilMode,
    segment::{
        dof::Dof,
        rbm::Rbm,
        units::{Rxyz, Txyz},
    },
};
use skyangle::SkyAngle;
use triangle_rs::Builder;

#[derive(Parser)]
pub struct Cli {
    /// thread pool size
    #[arg(long, default_value_t = 10usize)]
    n_thread: usize,
    #[command(subcommand)]
    pupil: Pupil,
}
#[derive(Subcommand)]
pub enum Pupil {
    /// full pupil
    Full,
    /// segmented pupil
    Segments,
    /// single segment
    Segment {
        /// segment ID
        #[arg(long)]
        id: i32,
        /// translation along X [micron]
        #[arg(long)]
        tx: Option<f64>,
        /// translation along Y [micron]
        #[arg(long)]
        ty: Option<f64>,
        /// translation along Z [micron]
        #[arg(long)]
        tz: Option<f64>,
        /// rotation along X [arcsec]
        #[arg(long)]
        rx: Option<f64>,
        /// rotation along Y [arcsec]
        #[arg(long)]
        ry: Option<f64>,
        // /// rotation along Z [arcsec]
        // #[arg(long)]
        // rz: Option<f64>,
    },
}

fn main() -> color_eyre::Result<()> {
    env_logger::init();

    let cli = Cli::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.n_thread)
        .build_global()
        .unwrap();

    let radius = 20. * 0.5;
    let perimeter = 2. * std::f64::consts::PI * radius;
    let delta = 0.5f64 * 4.;
    let n = (perimeter / delta).ceil() as usize;
    let nodes: Vec<_> = (0..n)
        .flat_map(|i| {
            let o = 2. * std::f64::consts::PI * i as f64 / n as f64;
            vec![radius * o.cos(), radius * o.sin()]
        })
        .collect();
    let tri = {
        let mut builder = Builder::new();
        builder.add_polygon(&nodes).set_switches("QDqa1.0");
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
        .map(|(z, a)| Pointing::new(SkyAngle::Arcminute(z as f32), SkyAngle::Radian(a as f32)))
        .collect();

    let pupil_mode = match cli.pupil {
        Pupil::Full => todo!(),
        Pupil::Segments => todo!(),
        Pupil::Segment {
            id,
            tx,
            ty,
            tz,
            rx,
            ry,
            // rz,
        } => {
            let mut rbm: Rbm = Default::default();
            if let Some(v) = tx {
                rbm += Rbm::from(Dof::Tx(Txyz::Mu(v)));
            }
            if let Some(v) = ty {
                rbm += Rbm::from(Dof::Ty(Txyz::Mu(v)));
            }
            if let Some(v) = tz {
                rbm += Rbm::from(Dof::Tz(Txyz::Mu(v)));
            }
            if let Some(v) = rx {
                rbm += Rbm::from(Dof::Rx(Rxyz::Arcsecond(v)));
            }
            if let Some(v) = ry {
                rbm += Rbm::from(Dof::Ry(Rxyz::Arcsecond(v)));
            }
            // if let Some(v) = rz {
            //     rbm += Rbm::from(Dof::Rz(Rxyz::Arcsecond(v)));
            // }
            PupilMode::segment(id, Mirror::m2(dbg!(rbm))).non_zeroed()
        }
    };
    let now = Instant::now();
    let probes = Probes::new(field_angles, 4, pupil_mode);
    println!("probed field in {:#?}", now.elapsed());

    let file = File::create(Path::new(env!("CARGO_MANIFEST_DIR")).join("gmt-full-field.pkl"))?;
    let mut buffer = BufWriter::new(file);
    serde_pickle::to_writer(&mut buffer, &probes, Default::default())?;
    Ok(())
}
