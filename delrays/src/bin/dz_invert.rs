//! Aberration coefficients from a saved `double_zernike` table
//!
//! ```sh
//! ENTRANCE=1 cargo run --release --bin double_zernike > dz.txt
//! cargo run --release --bin dz_invert -- dz.txt            # entrance-pupil layout
//! cargo run --release --bin dz_invert -- dz.txt --exit     # exit-pupil layout
//! cargo run --release --bin dz_invert < dz.txt             # from stdin
//! ```

use std::io::Read;

use delrays::inversion::{DoubleZernikes, Layout};

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let exit = args.iter().any(|a| a == "--exit");
    let text = match args.iter().find(|a| !a.starts_with("--")) {
        Some(path) => std::fs::read_to_string(path)?,
        None => {
            let mut text = String::new();
            std::io::stdin().read_to_string(&mut text)?;
            text
        }
    };
    let dz: DoubleZernikes = text.parse()?;
    let layout = if exit {
        Layout::exit()
    } else {
        Layout::entrance()
    };
    println!("{:.4}", dz.invert(&layout));
    Ok(())
}
