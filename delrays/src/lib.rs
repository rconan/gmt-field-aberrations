use crseo::{Builder, CrseoError, FromBuilder, raytracing::Rays};
use triangle_rs::Delaunay;

mod delaunay;
mod gmt;
mod opd;
pub use delaunay::Mesh;
pub use gmt::Gmt;
pub use opd::{Opd, Opds};

#[derive(Debug, thiserror::Error)]
pub enum DelraysError {
    #[error("crseo failure")]
    Crseo(#[from] crseo::CrseoError),
}

pub trait MaybeFrom<T> {
    type Error;
    fn maybe_from(value: T) -> Result<Self, Self::Error>
    where
        Self: Sized;
}
impl MaybeFrom<&Delaunay> for Rays {
    type Error = CrseoError;

    fn maybe_from(delaunay: &Delaunay) -> Result<Self, Self::Error> {
        Rays::builder()
            .xy(delaunay.vertex_iter().flatten().cloned().collect())
            .build()
    }
}

pub trait RayTracing {
    fn ray_tracing(self, rays: &mut Rays);
}
pub trait Trace {
    fn trace(&mut self, object: impl RayTracing);
    fn opds(&mut self) -> Opds;
}
impl Trace for Rays {
    fn trace(&mut self, object: impl RayTracing) {
        object.ray_tracing(self);
    }
    fn opds(&mut self) -> Opds {
        let opds = self.optical_path_difference();
        let xyz = self.coordinates();
        xyz.chunks(3).zip(opds.into_iter()).map(|(xyz, delta)| Opd {
            xyz: xyz.try_into().unwrap(),
            delta,
        }).collect()
    }
}
