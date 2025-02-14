use std::sync::Arc;

use serde::{Deserialize, Serialize};
use zernike::{gram_schmidt, jnm, zernike};

mod projection;
pub use projection::Projection;

/// Zernike modal basis
#[derive(Debug, Serialize, Deserialize, PartialEq)]
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
}
