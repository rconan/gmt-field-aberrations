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

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Field {
    zenith: SkyAngle<f32>,
    azimuth: SkyAngle<f32>,
    n_radial_order: u32,
    pupil_mode: PupilMode,
    coefs_format: CoefsFormat,
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
    pub fn new(n_radial_order: u32) -> Self {
        Self {
            n_radial_order,
            ..Default::default()
        }
    }
    pub fn pointing(self, zenith: SkyAngle<f32>, azimuth: SkyAngle<f32>) -> Self {
        Self {
            zenith,
            azimuth,
            ..self
        }
    }
    pub fn pupil_mode(self, pupil_mode: PupilMode) -> Self {
        Self { pupil_mode, ..self }
    }
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
