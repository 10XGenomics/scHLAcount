use failure::{Error, ResultExt};
use serde::{de::DeserializeOwned, Serialize};
use std::{
    io::{BufRead, BufReader, BufWriter, Read, Seek, Write},
    path::Path,
};

/// Open a reader for a text or gzip file
pub fn open_file(p: impl AsRef<Path>) -> Result<Box<dyn BufRead + Send>, Error> {
    let p = p.as_ref();

    let mut file =
        std::fs::File::open(p).with_context(|_| format!("Error opening file: {:?}", p))?;

    let mut buf = vec![0u8; 4];
    file.read_exact(&mut buf[..])?;
    file.seek(std::io::SeekFrom::Start(0))?;

    if &buf[0..2] == &[0x1F, 0x8B] {
        let gz = flate2::read::MultiGzDecoder::new(file);
        let buf_reader = BufReader::with_capacity(1 << 17, gz);
        Ok(Box::new(buf_reader))
    } else {
        let buf_reader = BufReader::with_capacity(32 * 1024, file);
        Ok(Box::new(buf_reader))
    }
}