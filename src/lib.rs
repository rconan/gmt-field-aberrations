mod field;
pub mod mirror;
pub mod segment;
mod zernike;
pub use field::Field;
use serde::{Deserialize, Serialize};
use skyangle::SkyAngle;
use zernike::{Projection, ZernikeBasis};

#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum Txyz {
    /// Microns
    Mu(f64),
    /// Nanometers
    Nm(f64),
}
impl Default for Txyz {
    fn default() -> Self {
        Self::Mu(0f64)
    }
}
impl Txyz {
    pub fn as_f64(&self) -> f64 {
        match self {
            Txyz::Mu(x) => *x * 1e-6,
            Txyz::Nm(x) => *x * 1e-9,
        }
    }
}
#[derive(Serialize, Deserialize, Debug, Default, Clone)]
pub struct Rbm {
    pub t_xyz: [Txyz; 3],
    pub r_xyz: [SkyAngle<f64>; 3],
}
#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum Mirror {
    M1(Rbm),
    M2(Rbm),
}
#[derive(Serialize, Deserialize, Default, Debug, Clone)]
pub enum PupilMode {
    #[default]
    Full,
    Segment {
        sid: i32,
        mirror: Mirror,
        zeroed: bool,
    },
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct CoefsFormat {
    pub width: usize,
    pub precision: usize,
}
impl Default for CoefsFormat {
    fn default() -> Self {
        Self {
            width: 8,
            precision: 0,
        }
    }
}
// pub fn zernike_basis(n_radial_order: u32, xy: &[[f64; 2]]) -> Vec<f64> {
//     let (mut r, o): (Vec<_>, Vec<_>) = xy
//         .iter()
//         .map(|x| (x[0].hypot(x[1]), x[1].atan2(x[0])))
//         .unzip();
//     let max_r = r
//         .iter()
//         .max_by(|x, y| x.partial_cmp(y).unwrap())
//         .copied()
//         .unwrap();
//     log::info!("Max. radius: {:.3}", max_r);
//     r.iter_mut().for_each(|r| *r /= max_r);

//     let (j, n, m) = jnm(n_radial_order);
//     let n_mode = j.len();
//     // println!("zernike mode #: {n_mode}");
//     let mut zern = vec![];
//     for (j, (n, m)) in j.into_iter().zip(n.into_iter().zip(m.into_iter())) {
//         zern.extend(
//             r.iter()
//                 .zip(o.iter())
//                 .map(|(r, o)| zernike(j, n, m, *r, *o)),
//         );
//     }
//     // dbg!((r.len(), zern.len(), r.len() * n_mode));
//     gram_schmidt(&zern, n_mode)
// }
