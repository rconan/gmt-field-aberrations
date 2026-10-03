use crseo::{
    Builder, CrseoError,
    raytracing::{Rays, RaysBuilder},
};
use triangle_rs::Delaunay;

use crate::{
    Mesh, Set,
    opd::{Opd, Opds},
};

/// Ray tracing through an optical system
pub trait RayTracing {
    /// Draws rays through the system
    fn ray_tracing(&mut self, rays: &mut Rays);
}
/// [Rays] extensions
pub trait Trace {
    /// Type of the object used to initialize the [Rays] coordinates
    type From;
    /// Type of the object that contains the [Rays] optical path differences
    type OpdData;
    /// Creates a new object from a mesh and a [Rays] builder
    fn from_mesh(mesh: &Self::From, builder: RaysBuilder) -> Result<Self, CrseoError>
    where
        Self: Sized;
    /// Traces the [Rays]
    fn trace<T: RayTracing>(&mut self, object: &mut T) -> &mut Self;
    /// Retrieves the [Rays] OPD
    fn opds(&mut self) -> Self::OpdData;
    /// Retrieves the [Rays] OPD with cooordinates centered on mesh origin
    fn opds_centered(&mut self, mesh: &Mesh) -> Option<Self::OpdData>;
    /// Retrieves the [Rays] OPD at the ray launch (entrance-pupil) coordinates,
    /// centered on mesh origin
    ///
    /// Unlike [Trace::opds_centered], the coordinates do not depend on the field:
    /// the segment is the same disc for every field point, as assumed by the
    /// double-Zernike expansion of the wave aberration
    fn opds_entrance(&mut self, mesh: &Mesh) -> Self::OpdData;
}
impl Trace for Rays {
    type From = Delaunay;
    type OpdData = Opds;
    fn from_mesh(mesh: &Self::From, builder: RaysBuilder) -> Result<Self, CrseoError> {
        builder
            .xy(mesh.vertex_iter().flatten().cloned().collect())
            .build()
    }
    fn trace<T: RayTracing>(&mut self, object: &mut T) -> &mut Self {
        object.ray_tracing(self);
        self
    }
    fn opds(&mut self) -> Self::OpdData {
        let opds = self.optical_path_difference();
        let xyz = self.coordinates();
        xyz.chunks(3)
            .zip(opds.into_iter())
            .map(|(xyz, delta)| Opd {
                xyz: xyz.try_into().unwrap(),
                delta,
            })
            .collect()
    }

    fn opds_centered(&mut self, mesh: &Mesh) -> Option<Self::OpdData> {
        let mut opds = self.opds();
        let cidx = mesh.origin_vertex_position();
        let Some(Opd { xyz: origin, .. }) = opds.get(cidx).cloned() else {
            return None;
        };
        opds.iter_mut()
            .for_each(|Opd { xyz, .. }| xyz.iter_mut().zip(origin).for_each(|(x, o)| *x -= o));
        Some(opds)
    }

    fn opds_entrance(&mut self, mesh: &Mesh) -> Self::OpdData {
        let origin = mesh
            .vertex_iter()
            .nth(mesh.origin_vertex_position())
            .map(|xy| [xy[0], xy[1]])
            .unwrap();
        mesh.vertex_iter()
            .zip(self.optical_path_difference())
            .map(|(xy, delta)| Opd {
                xyz: [xy[0] - origin[0], xy[1] - origin[1], 0f64],
                delta,
            })
            .collect()
    }
}
impl Trace for Set<Rays> {
    type From = Set<Mesh>;
    type OpdData = Set<Opds>;
    fn from_mesh(mesh: &Self::From, builder: RaysBuilder) -> Result<Self, CrseoError>
    where
        Self: Sized,
    {
        mesh.iter()
            .map(|del| Rays::from_mesh(del, builder.clone()))
            .collect()
    }

    fn trace<T: RayTracing>(&mut self, object: &mut T) -> &mut Self {
        self.iter_mut().for_each(|rays| {
            rays.trace(object);
        });
        self
    }

    fn opds(&mut self) -> Self::OpdData {
        self.iter_mut().map(|rays| rays.opds()).collect()
    }

    fn opds_centered(&mut self, mesh: &Mesh) -> Option<Self::OpdData> {
        self.iter_mut()
            .map(|rays| rays.opds_centered(mesh))
            .collect()
    }

    fn opds_entrance(&mut self, mesh: &Mesh) -> Self::OpdData {
        self.iter_mut()
            .map(|rays| rays.opds_entrance(mesh))
            .collect()
    }
}
