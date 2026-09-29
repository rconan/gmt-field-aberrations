use std::fs::File;

use crseo::
    raytracing::Rays
;
use delrays::{Gmt, MaybeFrom, Mesh, Trace};
use plotters::prelude::*;
use triangle_rs::Delaunay;


fn main() -> anyhow::Result<()> {
    let rim_diameter = 8.365;
    let delta_rim = 1f64 / 4f64;
    let delaunay = Delaunay::disc(rim_diameter, delta_rim);
    println!("{}", delaunay);

    let mut rays: Rays = MaybeFrom::maybe_from(&delaunay)?;

    let mut gmt = Gmt::new()?;
    rays.trace(&mut gmt);

    let opds = rays.opds();
    serde_pickle::to_writer(&mut File::create("opds.pkl")?, &opds, Default::default())?;

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
