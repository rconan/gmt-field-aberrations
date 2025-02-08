use serde::{Deserialize, Serialize};

mod field;
pub mod segment;
mod zernike;

pub use field::Field;
use segment::rbm::Rbm;

/// GMT mirror selection
#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum Mirror {
    M1(Rbm),
    M2(Rbm),
}
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
        /// either remove the collimated wavefront or not
        zeroed: bool,
    },
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
