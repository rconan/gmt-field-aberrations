use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::{
    zernike::{OpdToZernike, Projection},
    Field, Pointing, PupilMode,
};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Probes {
    pub fields: Vec<Field>,
    pub projections: Vec<Projection>,
}
impl Probes {
    pub fn new(
        field_angles: Vec<Pointing>,
        n_radial_order: u32,
        pupil_mode: PupilMode,
        opd_to_zern: OpdToZernike,
    ) -> Self {
        let (fields, projections): (Vec<Field>, Vec<Projection>) = field_angles
            .into_par_iter()
            .map(|pointing| {
                let field = Field::new(n_radial_order)
                    .pointing(pointing)
                    .pupil_mode(pupil_mode.clone());
                let zernp = field.zernike(opd_to_zern.clone()).unwrap();

                // let q = zernp.map(|zernp| (field, zernp));
                (field, zernp)
            })
            .unzip();
        Self {
            fields,
            projections,
        }
    }
    pub fn get(&self, i: usize) -> Option<(&Field, &Projection)> {
        self.fields.get(i).zip(self.projections.get(i))
    }
    pub fn len(&self) -> usize {
        self.fields.len()
    }
}
