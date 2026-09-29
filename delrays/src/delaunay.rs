use plotters::prelude::*;
use triangle_rs::{Builder, Delaunay};

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
        builder.add_polygon(&outer_rim).add_nodes(&[x0,y0]);
        builder
            .set_switches(&format!("Qpqa{}", triangle_area))
            .build()
    }
    fn plot(&self);
}

impl Mesh for Delaunay {
    fn plot(&self) {
        let fig = SVGBackend::new("mesh.svg", (768, 768)).into_drawing_area();
        fig.fill(&WHITE).unwrap();

        let xmin = self.x().into_iter().min_by(|a,b| a.partial_cmp(b).unwrap()).unwrap();
        let xmax = self.x().into_iter().max_by(|a,b| a.partial_cmp(b).unwrap()).unwrap();
        let xrange = xmin.floor()..xmax.ceil();
        let ymin = self.y().into_iter().min_by(|a,b| a.partial_cmp(b).unwrap()).unwrap();
        let ymax = self.y().into_iter().max_by(|a,b| a.partial_cmp(b).unwrap()).unwrap();
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
