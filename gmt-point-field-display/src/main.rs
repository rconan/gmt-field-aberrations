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

#[derive(Parser)]
pub struct Cli {
    /// zenith angle `[arcmin]`
    #[arg(short, long)]
    zenith: Vec<f64>,
    /// file name to save data to
    #[arg(short, long)]
    file: Option<String>,
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

    let mut zenith = cli.zenith.into_iter().map(|x| x as f32).cycle();
    let field_angles = vec![
        Pointing::new(
            SkyAngle::Arcminute(zenith.next().unwrap()),
            SkyAngle::Degree(0f32),
        ),
        Pointing::new(
            SkyAngle::Arcminute(zenith.next().unwrap()),
            SkyAngle::Degree(120f32),
        ),
        Pointing::new(
            SkyAngle::Arcminute(zenith.next().unwrap()),
            SkyAngle::Degree(240f32),
        ),
    ];

    const TLIM: f64 = 100f64;
    const RLIM: f64 = 10f64;
    let (id, mut rbm) = match cli.pupil {
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
            if let Some(_) = tx {
                rbm += Rbm::from(Dof::Tx(Txyz::Mu(-TLIM)));
            }
            if let Some(_) = ty {
                rbm += Rbm::from(Dof::Ty(Txyz::Mu(-TLIM)));
            }
            if let Some(_) = tz {
                rbm += Rbm::from(Dof::Tz(Txyz::Mu(-TLIM)));
            }
            if let Some(_) = rx {
                rbm += Rbm::from(Dof::Rx(Rxyz::Arcsecond(-RLIM)));
            }
            if let Some(_) = ry {
                rbm += Rbm::from(Dof::Ry(Rxyz::Arcsecond(-RLIM)));
            }
            // if let Some(v) = rz {
            //     rbm += Rbm::from(Dof::Rz(Rxyz::Arcsecond(v)));
            // }
            (id, rbm)
        }
    };
    let mut probes = vec![];
    let now = Instant::now();
    loop {
        let pupil_mode = PupilMode::segment(id, Mirror::m2(rbm.clone())).non_zeroed();
        probes.push(Probes::new(field_angles.clone(), 4, pupil_mode));
        match cli.pupil {
            Pupil::Full => todo!(),
            Pupil::Segments => todo!(),
            Pupil::Segment {
                tx,
                ty,
                tz,
                rx,
                ry,
                // rz,
                ..
            } => {
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
            }
        };
        if rbm.t_xyz.iter().find(|&&x| x.into_value() > TLIM).is_some()
            || rbm.r_xyz.iter().find(|&&x| x.into_value() > RLIM).is_some()
        {
            break;
        }
    }
    println!("probed field in {:#?}", now.elapsed());

    let file = File::create(
        Path::new(env!("CARGO_MANIFEST_DIR")).join(
            cli.file
                .unwrap_or("gmt-point_field-display.pkl".to_string()),
        ),
    )?;
    let mut buffer = BufWriter::new(file);
    serde_pickle::to_writer(&mut buffer, &probes, Default::default())?;
    Ok(())
}
