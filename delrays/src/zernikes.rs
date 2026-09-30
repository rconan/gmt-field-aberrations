use faer::{Mat, MatRef};
use serde::{Deserialize, Serialize};

use crate::Set;

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct Mode<C = f64> {
    pub jnm: (usize, usize, usize),
    pub(crate) mode: Vec<f64>,
    pub coef: C,
}

pub trait Zernike {
    type Coefs;
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
    fn jnm(&self) -> Vec<(usize, usize, usize)>;
    fn coefficients(&self) -> Self::Coefs;
    fn reduce_into(self, j: impl IntoIterator<Item = usize>) -> Self;
}

impl<C: Clone> Zernike for Set<Mode<C>> {
    type Coefs = Vec<C>;
    fn pseudo_inverse(&self) -> Mat<f64> {
        let ncols = self.len();
        let nrows = self[0].mode.len();
        let zerns: Vec<_> = self.iter().flat_map(|mode| mode.mode.to_vec()).collect();
        let mat = faer::MatRef::from_column_major_slice(&zerns, nrows, ncols);
        let svd = mat.svd().unwrap();
        svd.pseudoinverse()
    }

    fn jnm(&self) -> Vec<(usize, usize, usize)> {
        self.iter().map(|mode| mode.jnm.clone()).collect()
    }

    fn coefficients(&self) -> Self::Coefs {
        self.iter().map(|mode| mode.coef.clone()).collect()
    }

    fn reduce_into(self, j: impl IntoIterator<Item = usize>) -> Self {
        j.into_iter()
            .flat_map(|j| self.iter().filter(move |mode| mode.jnm.0 == j))
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

pub trait AsZernikes {
    type Into;
    fn as_zernikes(&self, n_radial_order: usize) -> Self::Into;
}

impl AsZernikes for Set<FieldZernike> {
    type Into = Set<Mode<Set<Mode>>>;
    fn as_zernikes(&self, n_radial_order: usize) -> Self::Into {
        let xy: Vec<_> = self
            .iter()
            .map(|field| {
                let (z, a) = field.za;
                let (s, c) = a.sin_cos();
                [z * c, z * s]
            })
            .collect();
        let mut zerns = Set::<Mode>::new(xy.iter().map(|xy| xy.as_slice()), n_radial_order);
        let n_mode = self.0[0].modes.len();
        let pinv = zerns.pseudo_inverse();
        let mut modes = vec![];
        for i in 0..n_mode {
            let c: Vec<_> = self.iter().map(|field| field.modes[i].coef).collect();
            let mat = MatRef::<f64>::from_column_major_slice(&c, c.len(), 1);
            let a = &pinv * mat;

            a.col(0).iter().zip(zerns.iter_mut()).for_each(|(c, z)| {
                z.coef = *c;
            });
            let Mode { jnm, mode, .. } = self[0].modes[i].clone();
            modes.push(Mode {
                jnm,
                mode,
                coef: zerns.clone(),
            });
        }
        modes.into_iter().collect()
    }
}
