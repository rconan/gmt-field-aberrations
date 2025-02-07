use std::fmt::Display;

use crseo::CrseoError;
use serde::{Deserialize, Serialize};
use skyangle::SkyAngle;

use crate::{mirror, segment, zernike::Projection, CoefsFormat, PupilMode};

#[derive(Debug, thiserror::Error)]
pub enum FieldError {
    #[error("failed to build CRSEO object")]
    Crseo(#[from] CrseoError),
    #[error("failed to create a segment")]
    Segment(#[from] geotrans::Error),
}
type Result<T> = std::result::Result<T, FieldError>;

/// Field location
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Field {
    /// zenith angle
    pub zenith: SkyAngle<f32>,
    /// azimuth angle
    pub azimuth: SkyAngle<f32>,
    /// Zernike basis radial order
    pub n_radial_order: u32,
    /// exit pupil definition
    pub pupil_mode: PupilMode,
    /// Zernike coefficients formatting
    pub coefs_format: CoefsFormat,
}
impl Default for Field {
    fn default() -> Self {
        Field {
            zenith: Default::default(),
            azimuth: Default::default(),
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
    pub fn pointing(self, zenith: SkyAngle<f32>, azimuth: SkyAngle<f32>) -> Self {
        Self {
            zenith,
            azimuth,
            ..self
        }
    }
    /// Configures the pupil
    pub fn pupil_mode(self, pupil_mode: PupilMode) -> Self {
        Self { pupil_mode, ..self }
    }
    /// Returns the [Projections] of the exit pupil wavefront onto the Zernike basis
    pub fn zernike(&self) -> Result<Projection> {
        let z = self.zenith.to_radians();
        let a = self.azimuth.to_radians();
        Ok(match &self.pupil_mode {
            PupilMode::Full => mirror::field_zernike(z, a, self.n_radial_order)?,
            PupilMode::Segment {
                sid,
                mirror,
                zeroed,
            } => segment::field_zernike(*sid, z, a, self.n_radial_order, mirror, *zeroed)?,
        }
        .coefficients_formatting(self.coefs_format.clone()))
    }
}
impl Display for Field {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", rson_rs::ser::pretty::to_string(self).unwrap())
    }
}
