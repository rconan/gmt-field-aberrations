//! Aberration coefficients from a saved `double_zernike` table
//!
//! ```sh
//! ENTRANCE=1 cargo run --release --bin double_zernike > dz.txt
//! cargo run --release --bin dz_invert -- dz.txt            # closed form, entrance-pupil layout
//! cargo run --release --bin dz_invert -- dz.txt --exit     # exit-pupil layout
//! cargo run --release --bin dz_invert < dz.txt             # from stdin
//!
//! # least squares with all the modes of the table (e.g. DZ_PUPIL_MODES=1-10 DZ_FIELD_ORDER=5)
//! cargo run --release --bin dz_invert -- dz.txt --lsq               # model k<=5, l<=8, m<=3
//! cargo run --release --bin dz_invert -- dz.txt --lsq --model 5,6,3
//! ```

use std::io::Read;

use delrays::inversion::{
    DoubleZernikes, Layout,
    lsq::{Coefficients, Design, terms},
};

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let flag = |name: &str| args.iter().any(|a| a == name);
    let value = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1))
    };
    let model = match value("--model") {
        Some(m) => {
            let v = m
                .split(',')
                .map(|x| x.trim().parse::<u32>())
                .collect::<Result<Vec<_>, _>>()?;
            anyhow::ensure!(v.len() == 3, "--model expects k,l,m");
            (v[0], v[1], v[2])
        }
        None => (5, 8, 3),
    };
    let path = args
        .iter()
        .enumerate()
        .find(|(i, a)| !a.starts_with("--") && (*i == 0 || args[i - 1] != "--model"))
        .map(|(_, a)| a);
    let text = match path {
        Some(path) => std::fs::read_to_string(path)?,
        None => {
            let mut text = String::new();
            std::io::stdin().read_to_string(&mut text)?;
            text
        }
    };
    let layout = if flag("--exit") {
        Layout::exit()
    } else {
        Layout::entrance()
    };
    if flag("--lsq") {
        let data: Coefficients = text.parse()?;
        let design = Design::new(
            terms(model.0, model.1, model.2),
            data.pupil(),
            data.field(),
            &layout,
        );
        println!("{:.4}", design.solve(&data)?);
    } else {
        let dz: DoubleZernikes = text.parse()?;
        println!("{:.4}", dz.invert(&layout));
    }
    Ok(())
}
