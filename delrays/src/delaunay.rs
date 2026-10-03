//! # Delaunay triangulation meshes

use std::{
    env, f64,
    fmt::Display,
    num::ParseIntError,
    ops::{Deref, DerefMut},
};

use plotters::prelude::*;
use triangle_rs::{Builder, Delaunay};

use crate::Set;

const CLEAR_APERTURE_DIAMETER: f64 = 8.365;

#[derive(Debug, thiserror::Error)]
pub enum MeshError {
    #[error(
        "failed to parse segment ID from environment variable SID, it should be an integer in the range [1,7]"
    )]
    SidParsing(#[from] ParseIntError),
    #[error("expected segment id in the range [1,7], found {0}")]
    WrongSid(i32),
}

/// Specialized mesh builders
pub struct Mesh {
    del: Delaunay,
    origin: [f64; 2],
}
impl Deref for Mesh {
    type Target = Delaunay;

    fn deref(&self) -> &Self::Target {
        &self.del
    }
}
impl DerefMut for Mesh {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.del
    }
}
impl Display for Mesh {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.del.fmt(f)
    }
}
impl Mesh {
    /// Returns the weights of the mesh lump mass matrix
    ///
    /// The weights are the sum of 1/3 of the areas of the triangles
    /// that vertices belongs to
    pub fn lump_mass_matrix_weights(&self) -> Vec<f64> {
        let mut weights = vec![0f64; self.n_vertices()];
        let areas = self.triangle_areas();
        for (indices, area) in self.triangle_iter().zip(areas.into_iter()) {
            for &idx in indices {
                weights[idx] += area / 3f64;
            }
        }
        weights
    }
    /// Makes a mesh for a disc
    pub fn disc(diameter: f64, perimeter_pitch: f64, origin: Option<[f64; 2]>) -> Self {
        let [x0, y0] = origin.unwrap_or([0f64; 2]);
        let mut builder = Builder::new();
        let triangle_area = perimeter_pitch.powi(2) * 3f64.sqrt() / 4f64;
        let n_rim = (std::f64::consts::PI * diameter / perimeter_pitch).round();
        let outer_rim: Vec<_> = (0..n_rim as usize)
            .flat_map(|i| {
                let o = 2. * std::f64::consts::PI * i as f64 / n_rim;
                let (s, c) = o.sin_cos();
                let radius = 0.5 * diameter;
                vec![radius * c + x0, radius * s + y0]
            })
            .collect();
        builder.add_polygon(&outer_rim).add_nodes(&[x0, y0]);
        let del = builder
            .set_switches(&format!("QDpqa{}", triangle_area))
            .build();
        Self {
            del,
            origin: [x0, y0],
        }
    }
    /// Makes a mesh for GMT segment #`id`
    pub fn gmt_segment_with_id(id: i32) -> Result<Self, MeshError> {
        if !(id > 0 && id < 8) {
            return Err(MeshError::WrongSid(id));
        };
        let rim_diameter = CLEAR_APERTURE_DIAMETER;
        let delta_rim = 1f64 / 4f64;
        let origin = (id < 7).then_some({
            let o = (3 - 2 * (id - 1)) as f64 * f64::consts::FRAC_PI_6;
            let (s, c) = o.sin_cos();
            [8.71 * c, 8.71 * s]
        });
        Ok(Self::disc(rim_diameter, delta_rim, origin))
    }
    /// Makes a mesh for a GMT segment
    ///
    // segment ID is set with environment variable SID (or set to 7 if not present)
    pub fn gmt_segment() -> Result<Self, MeshError> {
        let id = if let Ok(sid) = env::var("SID") {
            sid.parse::<i32>()?
        } else {
            7
        };
        Self::gmt_segment_with_id(id)
    }
    /// Makes a set of meshes for the GMT segments
    pub fn gmt() -> Set<Self> {
        let rim_diameter = CLEAR_APERTURE_DIAMETER;
        let delta_rim = 1f64 / 4f64;
        let mut segment = vec![];
        for id in 1..=7 {
            let origin = (id < 7).then_some({
                let o = (3 - 2 * (id - 1)) as f64 * f64::consts::FRAC_PI_6;
                let (s, c) = o.sin_cos();
                [8.71 * c, 8.71 * s]
            });
            segment.push(Self::disc(rim_diameter, delta_rim, origin));
        }
        Set(segment)
    }
    /// Return the index of the vertex that is the origin of the mesh
    pub fn origin_vertex_position(&self) -> usize {
        let [x0, y0] = self.origin;
        self.vertex_iter()
            .position(|xy| xy[0] == x0 && xy[1] == y0)
            .unwrap()
    }
}

impl Mesh {
    // Plots the mesh
    pub fn plot(&self) {
        let fig = SVGBackend::new("mesh.svg", (768, 768)).into_drawing_area();
        fig.fill(&WHITE).unwrap();

        let xmin = self
            .x()
            .into_iter()
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        let xmax = self
            .x()
            .into_iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        let xrange = xmin.floor()..xmax.ceil();
        let ymin = self
            .y()
            .into_iter()
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        let ymax = self
            .y()
            .into_iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        let yrange = ymin.floor()..ymax.ceil();

        let mut chart = ChartBuilder::on(&fig)
            .set_label_area_size(LabelAreaPosition::Left, 40)
            .set_label_area_size(LabelAreaPosition::Bottom, 40)
            .margin(20)
            .build_cartesian_2d(xrange.clone(), yrange)
            .unwrap();
        let mut mesh = chart.configure_mesh();
        mesh.draw().unwrap();

        self.triangle_iter()
            .map(|t| {
                t.iter()
                    .map(|&i| (self.x()[i], self.y()[i]))
                    .collect::<Vec<(f64, f64)>>()
            })
            .into_iter()
            .for_each(|v| {
                chart
                    .draw_series(LineSeries::new(
                        v.iter().cycle().take(4).map(|(x, y)| (*x, *y)),
                        &BLACK,
                    ))
                    .unwrap();
            });
    }
}

impl Set<Mesh> {
    // Plots the meshes
    pub fn plot(&self) {
        let fig = SVGBackend::new("mesh.svg", (768, 768)).into_drawing_area();
        fig.fill(&WHITE).unwrap();

        let xmin = self
            .iter()
            .map(|del| {
                del.x()
                    .into_iter()
                    .min_by(|a, b| a.partial_cmp(b).unwrap())
                    .unwrap()
            })
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        let xmax = self
            .iter()
            .map(|del| {
                del.x()
                    .into_iter()
                    .max_by(|a, b| a.partial_cmp(b).unwrap())
                    .unwrap()
            })
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        let ymin = self
            .iter()
            .map(|del| {
                del.y()
                    .into_iter()
                    .min_by(|a, b| a.partial_cmp(b).unwrap())
                    .unwrap()
            })
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        let ymax = self
            .iter()
            .map(|del| {
                del.y()
                    .into_iter()
                    .max_by(|a, b| a.partial_cmp(b).unwrap())
                    .unwrap()
            })
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();

        let xrange = xmin.floor()..xmax.ceil();
        let yrange = ymin.floor()..ymax.ceil();

        let mut chart = ChartBuilder::on(&fig)
            .set_label_area_size(LabelAreaPosition::Left, 40)
            .set_label_area_size(LabelAreaPosition::Bottom, 40)
            .margin(20)
            .build_cartesian_2d(xrange.clone(), yrange)
            .unwrap();
        let mut mesh = chart.configure_mesh();
        mesh.draw().unwrap();

        for del in self.iter() {
            del.triangle_iter()
                .map(|t| {
                    t.iter()
                        .map(|&i| (del.x()[i], del.y()[i]))
                        .collect::<Vec<(f64, f64)>>()
                })
                .into_iter()
                .for_each(|v| {
                    chart
                        .draw_series(LineSeries::new(
                            v.iter().cycle().take(4).map(|(x, y)| (*x, *y)),
                            &BLACK,
                        ))
                        .unwrap();
                });
        }
    }
}
