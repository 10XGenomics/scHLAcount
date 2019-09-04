use clap::{App, Arg};

use failure::{Error, format_err};
use simplelog::*;
use sprs::io::write_matrix_market;
use sprs::TriMat;
use terminal_size::{terminal_size, Width};

use std::collections::HashMap;
use std::fs::{create_dir_all, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::str::FromStr;
use std::string::String;

use log::{debug, info};
use human_panic::setup_panic;

mod mapping;
use mapping::{map_and_count_pseudo, map_and_count_sw, mapping_wrapper};

pub mod em;
use em::em_wrapper;

pub mod hla;
use hla::make_hla_index;

mod locus;
use locus::Locus;

pub mod config;
use config::Barcode;