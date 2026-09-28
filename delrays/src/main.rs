use std::fs::File;

use crseo::{
    Builder, FromBuilder,
    raytracing::{Conic, Rays},
};
use plotters::prelude::*;
use serde::{Deserialize, Serialize};
use triangle_rs as mesh;

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct Opds {
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
}

fn main() -> anyhow::Result<()> {
    let mut builder = mesh::Builder::new();
    let rim_diameter = 8.365;
    let delta_rim = 1f64 / 4f64;
    let triangle_area = delta_rim.powi(2) * 3f64.sqrt() / 4f64;
    let n_rim = (std::f64::consts::PI * rim_diameter / delta_rim).round();
    let outer_rim: Vec<_> = (0..n_rim as usize)
        .flat_map(|i| {
            let o = 2. * std::f64::consts::PI * i as f64 / n_rim;
            let (s, c) = o.sin_cos();
            let radius = 0.5 * rim_diameter;
            vec![radius * c, radius * s]
        })
        .collect();
    builder.add_polygon(&outer_rim);
    // let switches = format!("pDqa{}", 0.075 );
    let delaunay = builder
        .set_switches(&format!("pqa{}", triangle_area))
        .build();
    println!("{}", delaunay);

    let mut rays = Rays::builder()
        .xy(delaunay.vertex_iter().flatten().cloned().collect())
        .build()?;
    let mut m1 = Conic::builder()
        .conic_cst(1.-0.9982857)
        .curvature_radius(36.)
        .build()?;
    let mut m2 = Conic::builder()
        .conic_cst(1.-0.71692784)
        .curvature_radius(-4.1639009)
        .origin([0f64, 0f64, 20.26247614])
        .build()?;

    m1.trace(&mut rays);
    m2.trace(&mut rays);
    rays.to_sphere(-5.83, 2.197173);
    let opds = rays.optical_path_difference();
    serde_pickle::to_writer(
        &mut File::create("opds.pkl")?,
        &Opds {
            x: delaunay.x(),
            y: delaunay.y(),
            z: opds,
        },
        Default::default(),
    )?;

    let fig = SVGBackend::new("asm_facesheet.svg", (768, 768)).into_drawing_area();
    fig.fill(&WHITE).unwrap();

    let xyrange = -6f64..6f64;

    let mut chart = ChartBuilder::on(&fig)
        .set_label_area_size(LabelAreaPosition::Left, 40)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .margin(20)
        .build_cartesian_2d(xyrange.clone(), xyrange)
        .unwrap();
    let mut mesh = chart.configure_mesh();
    mesh.draw().unwrap();

    delaunay
        .triangle_iter()
        .map(|t| {
            t.iter()
                .map(|&i| (delaunay.x()[i], delaunay.y()[i]))
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

    // let mut colors = colorous::TABLEAU10.iter().cycle();
    // let this_color = colors.next().unwrap().as_tuple();

    // chart
    //     .draw_series(nodes.chunks(2).map(|xy| (xy[0], xy[1])).map(|point| {
    //         Circle::new(
    //             point,
    //             2,
    //             RED.filled(), // RGBColor(this_color.0, this_color.1, this_color.2).filled(),
    //         )
    //     }))
    //     .unwrap();

    Ok(())
}
