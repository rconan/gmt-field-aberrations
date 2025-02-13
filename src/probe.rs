use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use skyangle::SkyAngle;

use crate::{zernike::Projection, Field, PupilMode};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Probes {
    fields: Vec<Field>,
    projections: Vec<Projection>,
}
impl Probes {
    pub fn new(
        field_angles: Vec<(SkyAngle<f32>, SkyAngle<f32>)>,
        n_radial_order: u32,
        pupil_mode: PupilMode,
    ) -> Self {
        let (fields, projections): (Vec<Field>, Vec<Projection>) = field_angles
            .into_par_iter()
            .map(|(zenith, azimuth)| {
                let field = Field::new(n_radial_order)
                    .pointing(zenith, azimuth)
                    .pupil_mode(pupil_mode.clone());
                let zernp = field.zernike().unwrap();
                // let q = zernp.map(|zernp| (field, zernp));
                (field, zernp)
            })
            .unzip();
        Self {
            fields,
            projections,
        }
    }
}
