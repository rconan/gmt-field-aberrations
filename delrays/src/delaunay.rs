use std::{env, f64, num::ParseIntError};

use plotters::prelude::*;
use triangle_rs::{Builder, Delaunay};

use crate::Set;

#[derive(Debug, thiserror::Error)]
pub enum MeshError {
    #[error(
        "failed to parse segment ID from environment variable SID, it should be an integer in the range [1,7]"
    )]
    SidParsing(#[from] ParseIntError),
    #[error("expected segment id in the range [1,7], found {0}")]
    WrongSid(i32),
}

pub trait Mesh {
    fn disc(diameter: f64, perimeter_pitch: f64, origin: Option<[f64; 2]>) -> Delaunay {
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
        builder
            .set_switches(&format!("Qpqa{}", triangle_area))
            .build()
    }
    fn gmt_segment() -> Result<Delaunay, MeshError> {
        let id = if let Ok(sid) = env::var("SID") {
            sid.parse::<i32>()?
        } else {
            7
        };
        if !(id > 0 && id < 8) {
            return Err(MeshError::WrongSid(id));
        };
        let rim_diameter = 8.365;
        let delta_rim = 1f64 / 4f64;
        let origin = (id < 7).then_some({
            let o = (3 - 2 * (id - 1)) as f64 * f64::consts::FRAC_PI_6;
            let (s, c) = o.sin_cos();
            [8.71 * c, 8.71 * s]
        });
        Ok(Self::disc(rim_diameter, delta_rim, origin))
    }
    fn gmt() -> Set<Delaunay> {
        let rim_diameter = 8.365;
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
    fn plot(&self);
}

impl Mesh for Delaunay {
    fn plot(&self) {
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

impl Mesh for Set<Delaunay> {
    fn plot(&self) {
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
