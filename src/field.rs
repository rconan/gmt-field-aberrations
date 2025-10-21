use std::fmt::Display;

use crseo::CrseoError;
use serde::{Deserialize, Serialize};
use skyangle::SkyAngle;

use crate::{zernike::Projection, CoefsFormat, PupilMode};

mod mirror;
mod segment;

#[derive(Debug, thiserror::Error)]
pub enum FieldError {
    #[error("failed to build CRSEO object")]
    Crseo(#[from] CrseoError),
    #[error("failed to create a segment")]
    Segment(#[from] geotrans::Error),
}
type Result<T> = std::result::Result<T, FieldError>;

/// Field position \[zenith,azimuth\]
#[derive(Serialize, Deserialize, Debug, Default, Clone, Copy)]
pub struct Pointing {
    /// zenith angle
    pub zenith: SkyAngle<f32>,
    /// azimuth angle
    pub azimuth: SkyAngle<f32>,
}
impl Pointing {
    /// Creates an on-axis [Pointing] instance
    pub fn on_axis() -> Self {
        Default::default()
    }
    /// Creates a [Pointing] instance
    pub fn new(zenith: SkyAngle<f32>, azimuth: SkyAngle<f32>) -> Self {
        Self { zenith, azimuth }
    }
    /// Return the pointing direction in cartesian coordinates
    pub fn cartesian(&self) -> (f64, f64) {
        let (s, c) = (self.azimuth.to_radians() as f64).sin_cos();
        let z = self.zenith.into_arcmin().into_value() as f64;
        (z * c, z * s)
    }
}
impl From<(SkyAngle<f32>, SkyAngle<f32>)> for Pointing {
    fn from((zenith, azimuth): (SkyAngle<f32>, SkyAngle<f32>)) -> Self {
        Self { zenith, azimuth }
    }
}

/// Field location
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Field {
    /// pointing angle
    pub pointing: Pointing,
    /// Zernike basis number of radial order
    pub n_radial_order: u32,
    /// exit pupil definition
    pub pupil_mode: PupilMode,
    /// Zernike coefficients formatting
    pub coefs_format: CoefsFormat,
}
impl Default for Field {
    fn default() -> Self {
        Field {
            pointing: Default::default(),
            n_radial_order: 5,
            pupil_mode: Default::default(),
            coefs_format: Default::default(),
        }
    }
}
impl Field {
    /// Creates a new field instance
    pub fn new(n_radial_order: u32) -> Self {
        Self {
            n_radial_order,
            ..Default::default()
        }
    }
    /// Sets the GMT pointing direction
    pub fn pointing(self, pointing: impl Into<Pointing>) -> Self {
        Self {
            pointing: pointing.into(),
            ..self
        }
    }
    /// Configures the pupil
    pub fn pupil_mode(self, pupil_mode: PupilMode) -> Self {
        Self { pupil_mode, ..self }
    }
    /// Returns the [Projections] of the exit pupil wavefront onto the Zernike basis
    pub fn zernike(&self) -> Result<Projection> {
        let z = self.pointing.zenith.to_radians();
        let a = self.pointing.azimuth.to_radians();
        Ok(match &self.pupil_mode {
            PupilMode::Full => self.mirror(z, a)?,
            PupilMode::Segment {
                sid,
                mirror,
                zeroed,
            } => self.segment(*sid, z, a, mirror, *zeroed)?,
        }
        .coefficients_formatting(self.coefs_format.clone()))
    }
    /// Field pointing direction
    pub fn at(&self) -> &Pointing {
        &self.pointing
    }
}
impl Display for Field {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            ron::ser::to_string_pretty(self, Default::default()).unwrap()
        )
    }
}
