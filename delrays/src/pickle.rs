use std::{fs::File, io, path::Path};

use serde::Serialize;

#[derive(Debug, thiserror::Error)]
pub enum PickleError {
    #[error(transparent)]
    IO(#[from] io::Error),
    #[error("failed to serialize object")]
    Serialize(#[from] serde_pickle::Error),
}
pub trait Pickle {
    fn pickle(&self, path: impl AsRef<Path> + std::fmt::Debug) -> Result<(), PickleError>
    where
        Self: Serialize + Sized,
    {
        let mut file = File::create(path.as_ref()).map_err(|e| {
            print!("failed to create file {:?}", path);
            e
        })?;
        serde_pickle::to_writer(&mut file, self, Default::default())?;
        Ok(())
    }
}
