use std::{fmt::Display, sync::Arc};

use serde::{Deserialize, Serialize};
use zernike::{gram_schmidt, jnm, zernike};

use crate::CoefsFormat;

/// Zernike modal basis
#[derive(Debug, Serialize, Deserialize)]
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
/// Projection onto a Zernike basis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Projection {
    basis: Arc<ZernikeBasis>,
    coefficients: Vec<f64>,
    opd: Arc<[f32]>,
    coefs_format: CoefsFormat,
}
impl Projection {
    /// Creates a new [Projection] instance
    pub fn new(basis: impl Into<Arc<ZernikeBasis>>) -> Self {
        Self {
            basis: basis.into(),
            coefficients: Vec::new(),
            opd: Default::default(),
            coefs_format: Default::default(),
        }
    }
    /// Sets the formatting parameters for the Zernike coefficients
    pub fn coefficients_formatting(mut self, coefs_format: CoefsFormat) -> Self {
        self.coefs_format = coefs_format;
        self
    }
    /// Projects the wavefront onto the Zernike modes
    pub fn project(&mut self, opd: impl Into<Arc<[f32]>>) -> &mut Self {
        self.opd = opd.into();
        let n = self.opd.len();
        assert_eq!(n, self.basis.modes.len() / self.basis.n_mode);
        self.coefficients = self
            .basis
            .modes
            .chunks(self.opd.len())
            .map(|z| {
                z.iter()
                    .zip(self.opd.iter())
                    .map(|(z, o)| z * (*o as f64) * 1e9)
                    .sum::<f64>()
            })
            .map(|x| x / (n as f64).sqrt())
            .collect();
        self
    }
    /// Returns the Zernike projections coefficients
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }
    /// Returns the wavefront standard deviation in nanometers
    pub fn opd_std(&self) -> f32 {
        if self.opd.is_empty() {
            return 0f32;
        }
        let (mut mean, mut mean_squared) =
            self.opd.iter().fold((0f32, 0f32), |(mut x, mut y), c| {
                x += c;
                y += c * c;
                (x, y)
            });
        let n = self.opd.len() as f32;
        mean /= n;
        mean_squared /= n;
        1e9 * (mean_squared - mean * mean).sqrt()
    }
    /// Returns the root sum squared of the coefficients
    pub fn rms(&self) -> f64 {
        self.coefficients.iter().map(|&x| x * x).sum::<f64>().sqrt()
    }
    /// Returns the root sum squared of the coefficients ommiting the piston mode
    pub fn std(&self) -> f64 {
        self.coefficients
            .iter()
            .skip(1)
            .map(|&x| x * x)
            .sum::<f64>()
            .sqrt()
    }
}
impl Display for Projection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let CoefsFormat {
            width: w,
            precision: p,
        } = self.coefs_format;
        writeln!(
            f,
            r#"Mode # {:w$.0?}"#,
            (1..=self.basis.n_mode).collect::<Vec<_>>(),
        )?;
        write!(
            f,
            "Coefs. {:w$.p$?} {:4.0}/{:4.0}",
            self.coefficients(),
            self.opd_std(),
            self.std()
        )
    }
}
