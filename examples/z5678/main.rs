use std::vec;

use faer::{MatMut, MatRef};
use gmt_field_aberrations::{Field, Mirror, PupilMode, Rbm};
use gmt_field_aberrations::{Rxyz, Txyz};
use skyangle::{Conversion, SkyAngle};

#[derive(Debug, Clone, Copy)]
pub enum Dof {
    Tx(Txyz),
    Ty(Txyz),
    Rx(Rxyz),
    Ry(Rxyz),
}
impl Dof {
    pub fn into_iter() -> vec::IntoIter<Dof> {
        vec![
            Dof::Tx(Txyz::Mu(1.)),
            Dof::Ty(Txyz::Mu(1.)),
            Dof::Rx(Rxyz::Arcsecond(1.)),
            Dof::Ry(Rxyz::Arcsecond(1.)),
        ]
        .into_iter()
    }
    pub fn as_f64(&self) -> f64 {
        match self {
            Dof::Tx(txyz) | Dof::Ty(txyz) => txyz.as_f64(),
            Dof::Rx(sky_angle) | Dof::Ry(sky_angle) => sky_angle.to_radians(),
        }
    }
    pub fn len() -> usize {
        4
    }
}
impl From<Dof> for Rbm {
    fn from(value: Dof) -> Self {
        match value {
            Dof::Tx(txyz) => Rbm::t_x(txyz),
            Dof::Ty(txyz) => Rbm::t_y(txyz),
            Dof::Rx(sky_angle) => Rbm::r_x(sky_angle),
            Dof::Ry(sky_angle) => Rbm::r_y(sky_angle),
        }
    }
}

pub fn get_coefs(
    sid: i32,
    n_radial_order: u32,
    rbm: impl Into<Rbm> + Clone,
    azimuth: [i32; 3],
    skip: usize,
    take: usize,
    calib: &mut Vec<f64>,
    scale: f64,
) -> color_eyre::Result<()> {
    for az in azimuth {
        let pupil_mode = PupilMode::Segment {
            sid,
            mirror: Mirror::M2(rbm.clone().into()),
            zeroed: true,
        };
        let gs = Field::new(n_radial_order)
            .pointing(SkyAngle::Arcminute(6.), SkyAngle::Degree(az as f32))
            .pupil_mode(pupil_mode);
        // println!("{gs}");
        let zernp = gs.zernike()?;
        calib.extend(
            zernp
                .coefficients()
                .iter()
                .skip(skip)
                .take(take)
                .map(|x| scale * x),
        );
    }
    Ok(())
}

fn main() -> color_eyre::Result<()> {
    env_logger::init();

    let n_radial_order = 4;
    let sid = 1;
    let skip = 4;
    let take = 4;

    let mut calib = Vec::<f64>::new();
    let azimuth = [0, 120, 240];
    for dof in Dof::into_iter() {
        get_coefs(
            sid,
            n_radial_order,
            dof,
            azimuth,
            skip,
            take,
            &mut calib,
            1e-9 / dof.as_f64(),
        )?;
    }
    let mat = MatRef::from_column_major_slice(&calib, take * azimuth.len(), Dof::len());
    println!("{mat:+8.6?}");
    let svd = mat.thin_svd().unwrap();
    let s: Vec<_> = svd.S().column_vector().iter().collect();
    dbg!(&s);
    let cond = s[0] / *s.last().unwrap();
    println!("Cond.: {cond}");

    let mut data = Vec::<f64>::new();
    get_coefs(
        sid,
        n_radial_order,
        Rbm::t_x(Txyz::Mu(10.)),
        azimuth,
        skip,
        take,
        &mut data,
        1e-9,
    )?;
    let n = data.len();
    let rhs = MatMut::from_column_major_slice_mut(&mut data, n, 1);
    let is: Vec<_> = s.iter().map(|x| x.recip()).collect();
    let imat = svd.V()
        * MatRef::from_column_major_slice(&is, is.len(), 1)
            .col(0)
            .as_diagonal()
        * svd.U().transpose();
    let x = imat * rhs;
    let x: Vec<_> = x.col(0).iter().collect();
    let rbm = Rbm::t_x(Txyz::to_mu(*x[0]))
        + Rbm::t_y(Txyz::to_mu(*x[1]))
        + Rbm::r_x(Rxyz::Arcsecond(x[2].to_arcsec()))
        + Rbm::r_y(Rxyz::Arcsecond(x[3].to_arcsec()));
    dbg!(&rbm);

    Ok(())
}
