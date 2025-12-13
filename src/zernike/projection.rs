use std::{fmt::Display, sync::Arc};

use faer::linalg::svd::SvdError;
use serde::{Deserialize, Serialize};

use crate::CoefsFormat;

use super::ZernikeBasis;

#[derive(Debug, Clone)]
pub enum OpdToZernike {
    Projection,
    LeastSquareFit,
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
    pub fn least_square_fit(&mut self, opd: impl Into<Arc<[f32]>>) -> Result<&mut Self, SvdError> {
        self.opd = opd.into();
        let n = self.opd.len();
        assert_eq!(n, self.basis.modes.len() / self.basis.n_mode);
        let mut iter = self.opd.iter().map(|&x| 1e9 * x as f64);
        let opd = faer::Mat::from_fn(n, 1, |_, _| iter.next().unwrap());
        let svd = self.basis.mat().svd()?;
        let c = svd.pseudoinverse() * opd;
        self.coefficients = c.col_as_slice(0).to_vec();

        Ok(self)
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
