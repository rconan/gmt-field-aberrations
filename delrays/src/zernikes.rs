use faer::Mat;
use serde::{Deserialize, Serialize};

use crate::Set;

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct Mode {
    pub(crate) jnm: (usize, usize, usize),
    pub(crate) mode: Vec<f64>,
    pub(crate) coef: f64,
}

pub trait Zernike {
    fn new<'a>(xy: impl Iterator<Item = &'a [f64]>, n_radial_order: usize) -> Set<Mode> {
        let (mut r, o): (Vec<f64>, Vec<f64>) = xy
            .map(|xy| (xy[0].hypot(xy[1]), xy[1].atan2(xy[0])))
            .unzip();
        let max_r = r
            .iter()
            .cloned()
            .max_by(|x, y| x.partial_cmp(y).unwrap())
            .unwrap();
        r.iter_mut().for_each(|r| *r /= max_r);

        let (j, n, m) = zernike::jnm(n_radial_order as u32);
        j.into_iter()
            .zip(n.into_iter())
            .zip(m.into_iter())
            .map(|((j, n), m)| {
                let mode = r
                    .iter()
                    .zip(&o)
                    .map(|(&r, &o)| zernike::zernike(j, n, m, r, o))
                    .collect::<Vec<_>>();
                Mode {
                    jnm: (j as usize, n as usize, m as usize),
                    mode,
                    coef: 0f64,
                }
            })
            .collect()
    }
    fn pseudo_inverse(&self) -> Mat<f64>;
    fn coefficients(&self) -> Vec<f64>;
    fn reduce_into(self, j: &[usize]) -> Self;
}

impl Zernike for Set<Mode> {
    fn pseudo_inverse(&self) -> Mat<f64> {
        let ncols = self.len();
        let nrows = self[0].mode.len();
        let zerns: Vec<_> = self.iter().flat_map(|mode| mode.mode.to_vec()).collect();
        let mat = faer::MatRef::from_column_major_slice(&zerns, nrows, ncols);
        let svd = mat.svd().unwrap();
        svd.pseudoinverse()
    }

    fn coefficients(&self) -> Vec<f64> {
        self.iter().map(|mode| mode.coef).collect()
    }

    fn reduce_into(self, j: &[usize]) -> Self {
        j.iter()
            .flat_map(|j| self.iter().filter(|mode| mode.jnm.0 == *j))
            .cloned()
            .collect()
    }
}

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct FieldZernike {
    za: (f64, f64),
    modes: Set<Mode>,
}

impl FieldZernike {
    pub fn new(za: (f64, f64), modes: Set<Mode>) -> Self {
        Self { za, modes }
    }
}
