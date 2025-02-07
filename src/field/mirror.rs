use crseo::{gmt, utilities::MaskFilter, Builder, FromBuilder, Source};

use crate::{field::FieldError, Field, Projection, ZernikeBasis};

type Result<T> = std::result::Result<T, FieldError>;

impl Field {
    pub fn mirror(&self, z: f32, a: f32) -> Result<Projection> {
        let mut gmt = gmt!();
        // let mut src = source!();
        let mut src = Source::builder().zenith_azimuth(vec![z], vec![a]).build()?;
        src.through(&mut gmt).xpupil();
        let wfe_rms = src.wfe_rms_10e(-9);
        log::info!("WFE RMS: {:?}nm", wfe_rms);

        let xy: Vec<_> = src
            .rays()
            .mask()
            .filter(src.rays().coordinates().chunks(3))
            .map(|c| [c[0], c[1]])
            .collect();
        let zern = ZernikeBasis::new(self.n_radial_order, &xy);
        let opd: Vec<_> = src
            .rays()
            .mask()
            .filter(src.phase().iter())
            .cloned()
            .collect();
        let mut zernp = Projection::new(zern);
        zernp.project(&*opd);
        Ok(zernp)
    }
}
