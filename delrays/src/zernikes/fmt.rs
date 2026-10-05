use std::fmt::Display;

use super::{SegmentsDoubleZernikes, Zernike};

#[derive(Debug, Clone)]
pub struct SegmentsDoubleZernikesFormat<'a> {
    data: &'a SegmentsDoubleZernikes,
    scale: f64,
    width: usize,
    precision: usize,
    field_modes: usize,
}

impl<'a> SegmentsDoubleZernikesFormat<'a> {
    pub fn scale(mut self, scale: f64) -> Self {
        self.scale = scale;
        self
    }
    pub fn width(mut self, width: usize) -> Self {
        self.width = width;
        self
    }
    pub fn precision(mut self, precision: usize) -> Self {
        self.precision = precision;
        self
    }
    /// Number of field modes printed (default: 8)
    pub fn field_modes(mut self, field_modes: usize) -> Self {
        self.field_modes = field_modes;
        self
    }
}

impl<'a> From<&'a SegmentsDoubleZernikes> for SegmentsDoubleZernikesFormat<'a> {
    fn from(data: &'a SegmentsDoubleZernikes) -> Self {
        Self {
            data,
            scale: 1e9,
            width: 6,
            precision: 0,
            field_modes: 8,
        }
    }
}
impl<'a> Display for SegmentsDoubleZernikesFormat<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        'pupil: for j in 0.. {
            'field: for k in 0.. {
                for (i, field_zern_coefs) in self.data.iter().enumerate() {
                    let Some(fzc) = field_zern_coefs.get(j) else {
                        break 'pupil;
                    };
                    let coef = fzc.coef.clone().reduce_into(1..=self.field_modes);
                    let jnm = coef.jnm();
                    let coef = coef.coefficients();
                    let Some((jnm, c)) = jnm.get(k).zip(coef.get(k)) else {
                        break 'field;
                    };
                    if i == 0 {
                        if k == 0 {
                            writeln!(f, "{:?}", fzc.jnm)?;
                        }
                        write!(f, " {:2?}: ", jnm)?;
                    }
                    write!(f, "{:+1$.2$}", c * self.scale, self.width, self.precision)?;
                }
                writeln!(f)?;
            }
        }
        Ok(())
    }
}
