use std::sync::Arc;

use serde::{Deserialize, Serialize};
use zernike::{gram_schmidt, jnm, zernike};

mod projection;
pub use projection::{OpdToZernike, Projection};

/// Zernike modal basis
#[derive(Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct ZernikeBasis {
    xy: Arc<[[f64; 2]]>,
    n_radial_order: u32,
    n_mode: usize,
    modes: Vec<f64>,
    sids: Option<Vec<i32>>,
}
impl ZernikeBasis {
    /// Creates a new Zernike basis instance with the number of radial orders
    /// and \[x,y\] coordinates
    pub fn new(n_radial_order: u32, xy: &[[f64; 2]], sids: Option<Vec<i32>>) -> Self {
        let (mut r, o): (Vec<_>, Vec<_>) = xy
            .iter()
            .map(|x| (x[0].hypot(x[1]), x[1].atan2(x[0])))
            .unzip();
        let max_r = r
            .iter()
            .max_by(|x, y| x.partial_cmp(y).unwrap())
            .copied()
            .unwrap();
        log::info!("Max. radius: {:.3}", max_r);
        r.iter_mut().for_each(|r| *r /= max_r);

        let (j, n, m) = jnm(n_radial_order);
        let n_mode = j.len();
        // println!("zernike mode #: {n_mode}");
        let mut zern = vec![];
        for (j, (n, m)) in j.into_iter().zip(n.into_iter().zip(m.into_iter())) {
            zern.extend(
                r.iter()
                    .zip(o.iter())
                    .map(|(r, o)| zernike(j, n, m, *r, *o)),
            );
        }
        // dbg!((r.len(), zern.len(), r.len() * n_mode));
        Self {
            xy: xy.into(),
            n_radial_order,
            modes: gram_schmidt(&zern, n_mode),
            n_mode,
            sids,
        }
    }
    pub fn builder(n_radial_order: u32, xy: &[[f64; 2]]) -> ZernikeBasisBuilder {
        ZernikeBasisBuilder {
            built: ZernikeBasis {
                xy: xy.into(),
                n_radial_order,
                ..Default::default()
            },
            ..Default::default()
        }
    }
    pub fn mat(&self) -> faer::MatRef<'_, f64> {
        faer::MatRef::from_column_major_slice(&self.modes, self.xy.len(), self.n_mode)
    }
}

#[derive(Debug, Default)]
pub struct ZernikeBasisBuilder {
    built: ZernikeBasis,
    gramschmidt: bool,
}
impl ZernikeBasisBuilder {
    pub fn gramschmidt(mut self, gramschmidt: bool) -> Self {
        self.gramschmidt = gramschmidt;
        self
    }
    pub fn build(self) -> ZernikeBasis {
        let zernb = self.built;
        let (mut r, o): (Vec<_>, Vec<_>) = zernb
            .xy
            .iter()
            .map(|x| (x[0].hypot(x[1]), x[1].atan2(x[0])))
            .unzip();
        let max_r = r
            .iter()
            .max_by(|x, y| x.partial_cmp(y).unwrap())
            .copied()
            .unwrap();
        log::info!("Max. radius: {:.3}", max_r);
        r.iter_mut().for_each(|r| *r /= max_r);

        let (j, n, m) = jnm(zernb.n_radial_order);
        let n_mode = j.len();
        // println!("zernike mode #: {n_mode}");
        let mut zern = vec![];
        for (j, (n, m)) in j.into_iter().zip(n.into_iter().zip(m.into_iter())) {
            zern.extend(
                r.iter()
                    .zip(o.iter())
                    .map(|(r, o)| zernike(j, n, m, *r, *o)),
            );
        }
        // dbg!((r.len(), zern.len(), r.len() * n_mode));
        ZernikeBasis {
            modes: if self.gramschmidt {
                gram_schmidt(&zern, n_mode)
            } else {
                zern
            },
            n_mode,
            sids: None,
            ..zernb
        }
    }
}

#[cfg(test)]
mod tests {
    use std::iter;

    use super::*;

    #[test]
    fn zernike() {
        let n = 64;
        let mut xy = vec![];
        for i in 0..n {
            let x = 2f64 * i as f64 / (n - 1) as f64 - 1f64;
            for j in 0..n {
                let y = 2f64 * j as f64 / (n - 1) as f64 - 1f64;
                let r = x.hypot(y);
                if r <= 1f64 {
                    xy.push([x, y]);
                }
            }
        }
        let zern = ZernikeBasis::builder(1, &xy).build();
        dbg!(&zern.modes[..5]);
    }
    #[test]
    fn project() {
        let n = 64;
        let mut xy = vec![];
        for i in 0..n {
            let x = 2f64 * i as f64 / (n - 1) as f64 - 1f64;
            for j in 0..n {
                let y = 2f64 * j as f64 / (n - 1) as f64 - 1f64;
                let r = x.hypot(y);
                if r <= 1f64 {
                    xy.push([x, y]);
                }
            }
        }
        let zern = ZernikeBasis::builder(1, &xy).build();
        dbg!(&zern.modes[..5]);

        let mut projection = Projection::new(zern);
        let opd: Vec<_> = iter::repeat(1e-9f32).take(xy.len()).collect::<Vec<_>>();
        projection.least_square_fit(opd).unwrap();
        dbg!(projection.coefficients());
    }
}
