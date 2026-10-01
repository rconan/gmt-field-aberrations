//! # Field aberrations for 2 mirror telescopes

pub mod delaunay;
pub mod gmt;
pub mod opd;
pub mod pickle;
mod set;
mod trace;
pub mod zernikes;

#[doc(inline)]
pub use delaunay::Mesh;
#[doc(inline)]
pub use gmt::Gmt;
pub use set::Set;
#[doc(inline)]
pub use trace::Trace;

#[derive(Debug, thiserror::Error)]
pub enum DelraysError {
    #[error("crseo failure")]
    Crseo(#[from] crseo::CrseoError),
}
