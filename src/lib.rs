use serde::{Deserialize, Serialize};

mod field;
mod mirror;
pub mod segment;
mod zernike;

pub use field::{Field, FieldError};
pub use mirror::Mirror;

/// GMT pupil either full or restrict to segment
#[derive(Serialize, Deserialize, Default, Debug, Clone)]
pub enum PupilMode {
    #[default]
    Full,
    Segment {
        /// segment id \[1,7\]
        sid: i32,
        /// segment [Rbm]
        mirror: Mirror,
        /// either remove the collimated (`true`) wavefront or not (`false`)
        zeroed: bool,
    },
}
impl PupilMode {
    pub fn m1(sid: i32) -> Self {
        Self::Segment {
            sid,
            mirror: Mirror::M1(Default::default()),
            zeroed: true,
        }
    }
    pub fn m2(sid: i32) -> Self {
        Self::Segment {
            sid,
            mirror: Mirror::M2(Default::default()),
            zeroed: true,
        }
    }
}
/// Zernike coefficients formatting
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct CoefsFormat {
    /// format width
    pub width: usize,
    /// format precision
    pub precision: usize,
}
impl Default for CoefsFormat {
    fn default() -> Self {
        Self {
            width: 8,
            precision: 0,
        }
    }
}
