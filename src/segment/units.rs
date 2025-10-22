use std::ops::{AddAssign, Neg};

use serde::{Deserialize, Serialize};
use skyangle::SkyAngle;

/// Units for RBM translations
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
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
    /// Returns the translation value in meters
    pub fn as_f64(&self) -> f64 {
        match self {
            Txyz::Mu(x) => *x * 1e-6,
            Txyz::Nm(x) => *x * 1e-9,
        }
    }
    pub fn to_nm(v: f64) -> Txyz {
        Txyz::Nm(1e9 * v)
    }
    pub fn to_mu(v: f64) -> Txyz {
        Txyz::Mu(1e6 * v)
    }
    pub fn into_value(self) -> f64 {
        match self {
            Txyz::Mu(v) => v,
            Txyz::Nm(v) => v,
        }
    }
}
pub type Rxyz = SkyAngle<f64>;
impl Neg for Txyz {
    type Output = Self;

    fn neg(self) -> Self::Output {
        match self {
            Txyz::Mu(x) => Txyz::Mu(-x),
            Txyz::Nm(x) => Txyz::Nm(-x),
        }
    }
}
impl AddAssign for Txyz {
    fn add_assign(&mut self, rhs: Self) {
        match self {
            Txyz::Mu(v) => *v += rhs.into_value(),
            Txyz::Nm(v) => *v += rhs.into_value(),
        }
    }
}
