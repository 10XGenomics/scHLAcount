extern crate bio;
extern crate clap;
extern crate debruijn;
extern crate debruijn_mapping;
#[macro_use]
extern crate failure;
#[macro_use]
extern crate human_panic;
extern crate itertools;
#[macro_use]
extern crate log;
extern crate regex;
extern crate rust_htslib;
#[macro_use]
extern crate serde;
extern crate simplelog;
extern crate smallvec;
extern crate sprs;
extern crate tempfile;
extern crate terminal_size;

use clap::{Arg, App};

use failure::Error;
use simplelog::*;
use sprs::io::write_matrix_market;
use sprs::TriMat;
use terminal_size::{Width, terminal_size};

use std::collections::HashMap;
use std::fs::{File,create_dir_all};
use std::process;
use std::path::{Path,PathBuf};
use std::str::FromStr;
use std::string::String;
use std::io::{BufRead, BufReader, BufWriter, Write};


mod mapping;
use mapping::{mapping_wrapper, map_and_count_pseudo, map_and_count_sw};

mod em;
use em::{em_wrapper};

mod hla;
use hla::make_hla_index;

mod locus;
use locus::Locus;

mod config;
use config::Barcode;

fn get_args() -> clap::App<'static, 'static> {
    App::new("scHLAcount")
    .set_term_width(if let Some((Width(w), _)) = terminal_size() { w as usize } else { 120 })
    .version("DEV")
    .author("Charlotte Darby <cdarby@jhu.edu> and Ian Fiddes <ian.fiddes@10xgenomics.com> and Patrick Marks <patrick@10xgenomics.com>")
    .about("HLA genotyping and allele-specific expression for single-cell RNA sequencing")
    // Required parameters
    .arg(Arg::with_name("bam")
         .short("b")
         .long("bam")
         .value_name("FILE")
         .help("Cellranger BAM file")
         .required(true))
    .arg(Arg::with_name("cell_barcodes")
         .short("c")
         .long("cell-barcodes")
         .value_name("FILE")
         .help("File with cell barcodes to be evaluated")
         .required(true))
    // Output parameters (optional)
    .arg(Arg::with_name("out_dir")
         .short("o")
         .long("out-dir")
         .value_name("OUTPUT_DIR")
         .default_value("hla-typer-results"))
    .arg(Arg::with_name("pseudoalignertmp")
         .long("pl-tmp")
         .value_name("PSEUDOALIGNER_TMP")
         .default_value(".")
         .help("Directory to write the pseudoaligner temporary files generated"))
    // Input parameters (optional)
    .arg(Arg::with_name("fastagenomic")
         .short("g")
         .long("fasta-genomic")
         .value_name("FILE")
         .help("Multi-FASTA file with genomic sequence of each allele")
         .default_value(""))
    .arg(Arg::with_name("fastacds")
         .short("f")
         .long("fasta-cds")
         .value_name("FILE")
         .help("Multi-FASTA file with CDS sequence of each allele")
         .default_value(""))
     .arg(Arg::with_name("hladbdir")
         .short("d")
         .long("hladb-dir")
         .value_name("PATH")
         .help("Directory of the IMGT-HLA database")
         .default_value(""))
     .arg(Arg::with_name("hlaindex")
         .short("i")
         .long("hla-index")
         .value_name("FILE")
         .help("debruijn_mapping pseudoalignment index file constructed from IMGT-HLA database")
         .default_value(""))
    // Configuration parameters (optional)
    .arg(Arg::with_name("region")
         .short("r")
         .long("region")
         .value_name("STRING")
         .help("Samtools-format region string of reads to use")
         .default_value("6:28510120-33480577"))
    .arg(Arg::with_name("log_level")
         .long("log-level")
         .possible_values(&["info", "debug", "error"])
         .default_value("error")
         .help("Logging level"))
    .arg(Arg::with_name("primary_alignments")
         .long("primary-alignments")
         .help("If specified, will use primary alignments only"))
    .arg(Arg::with_name("exact_count")
         .long("use-exact-count")
         .help("If specified, will use exact alignment to allele sequences to count moleucles (very slow!)"))
    .arg(Arg::with_name("unmapped")
         .long("unmapped")
         .help("If specified, will also use unmapped reads for genotyping (very slow!)"))
}


fn main() {
    setup_panic!();  // pretty panics for users
    let mut cli_args = Vec::new();
    for arg in std::env::args_os() {
        cli_args.push(arg.into_string().unwrap());
    }
    let res = _main(cli_args);
    if let Err(e) = res {
        println!("Failed with error: {}", e);
        process::exit(1);
    }
}

// constructing a _main allows for us to run regression tests way more easily
#[allow(clippy::cognitive_complexity)] 
fn _main(cli_args: Vec<String>) -> Result<(), Error> {
    let args = get_args().get_matches_from(cli_args);
    let fasta_file_gen = args.value_of("fastagenomic").unwrap_or_default();
    let fasta_file_cds = args.value_of("fastacds").unwrap_or_default();
    let hla_db_dir = args.value_of("hladbdir").unwrap_or_default();
    let pseudoaligner_tmpdir = args.value_of("pseudoalignertmp").unwrap_or_default();
    let hla_index = args.value_of("hlaindex").unwrap_or_default();
    let bam_file = args.value_of("bam").expect("You must provide a BAM file");
    let cell_barcodes = args.value_of("cell_barcodes").expect("You must provide a cell barcodes file");
    let region = args.value_of("region").unwrap_or_default();
    let out_dir = args.value_of("out_dir").unwrap_or_default();
    let primary_only = args.is_present("primary_alignments");
    let exact_count = args.is_present("exact_count");
    let use_unmapped = args.is_present("unmapped");
    let ll = args.value_of("log_level").unwrap();

    let ll = match ll {
        "info" => LevelFilter::Info,
        "debug" => LevelFilter::Debug,
        "error" => LevelFilter::Error,
        &_ => { return Err(format_err!("Log level must be 'info', 'debug', or 'error'")); }
    };
    let _ = SimpleLogger::init(ll, Config::default());
    
    check_inputs_exist(bam_file, cell_barcodes, out_dir)?;
    let region : Locus = Locus::from_str(region).expect("Failed to parse region string");   
    let (mut genomic, mut cds) : (PathBuf,PathBuf) = (PathBuf::from(fasta_file_gen), PathBuf::from(fasta_file_cds));
    let bam_file = PathBuf::from(bam_file);
    // If the CDS or genomic FASTA files were not provided, generate them
    if fasta_file_gen.is_empty() || fasta_file_cds.is_empty() {
        check_inputs_exist_pseudoaligner(pseudoaligner_tmpdir)?;
        let p_tmpdir : PathBuf = [pseudoaligner_tmpdir, "pseudoaligner"].iter().collect();
        if hla_db_dir.is_empty() {
            return Err(format_err!("Must provide either -d (database directory) or both -g and -f (FASTA sequences)."));
        }
        check_inputs_exist_hla_db(hla_db_dir)?;
        let (counts_file, hla_index1) : (PathBuf, PathBuf) =
        // If the index for pseudoalignment was not provided, generate it
        if hla_index.is_empty() {
            let db_fasta: PathBuf = [hla_db_dir, "hla_nuc.fasta"].iter().collect();
            let allele_status: PathBuf = [hla_db_dir, "Allele_status.txt"].iter().collect();
            let p_hla_index: PathBuf = [pseudoaligner_tmpdir,"hla_nuc.fasta.idx"].iter().collect();
            let hla_index_generated: PathBuf = make_hla_index(db_fasta,p_hla_index,allele_status).expect("Pseudoaligner index building failed");
            (mapping_wrapper(hla_index_generated.clone(), p_tmpdir, bam_file.clone(), Some(&region), use_unmapped).expect("Pseudoalignment mapping failed"), hla_index_generated)
        } else {
            check_inputs_exist_hla_idx(hla_index)?;
            (mapping_wrapper(PathBuf::from(hla_index), p_tmpdir, bam_file.clone(), Some(&region), use_unmapped).expect("Pseudoalignment mapping failed"), PathBuf::from(hla_index))
        };

        let cdsdb : PathBuf = [hla_db_dir, "hla_nuc.fasta"].iter().collect();
        let gendb : PathBuf = [hla_db_dir, "hla_gen.fasta"].iter().collect();
        if let Ok(i) = em_wrapper(hla_index1, counts_file, cdsdb, gendb, pseudoaligner_tmpdir) {
            genomic = i.0;
            cds = i.1;
        } else {
            return Err(format_err!("EM/pseudoalignment-consensus step failed"));
        }
    }
    check_inputs_exist_fasta(&cds, &genomic)?;
    
    let cell_barcodes = load_barcodes(&cell_barcodes).unwrap();
    let (entries, nrows, metrics, rownames) = if exact_count { 
            map_and_count_sw(bam_file, &out_dir, &cell_barcodes, &region, cds, genomic, primary_only)? 
        } else {
            map_and_count_pseudo(bam_file, &out_dir, &cell_barcodes, &region, cds, genomic, primary_only)?
        };
    
    
    info!("Initialized a {} features x {} cell barcodes matrix", nrows, cell_barcodes.len());
    
    let mut matrix : TriMat<usize> = TriMat::new((nrows, cell_barcodes.len()));
    
    info!("Number of alignments evaluated (with 'CB' and 'UB' tags): {}", metrics.num_reads);
    if primary_only {
        info!("Number of alignments skipped due to not being primary: {}", metrics.num_non_primary);
    } else {
        info!("Number of alignments that were not primary: {}", metrics.num_non_primary);
    }
    info!("Number of alignments skipped due to not being associated with a cell barcode in the list provided: {}", metrics.num_not_cell_bc);
    info!("Number of reads with no alignment score above threshold: {}", metrics.num_not_aligned);
    info!("Number of alignments to CDS sequence: {}", metrics.num_cds_align);
    info!("Number of alignments to genomic sequence: {}", metrics.num_gen_align);
        
    for e in entries {
        matrix.add_triplet(e.row as usize, e.column, e.value);
    }
    
    let d: PathBuf = [out_dir,"count_matrix.mtx"].iter().collect();
    write_matrix_market(d.to_str().unwrap(), &matrix)?;
    debug!("Wrote reference matrix file");
    
    let d: PathBuf = [out_dir,"summary.tsv"].iter().collect();
    let mut summary_file = BufWriter::new(File::create(d).unwrap());
    let m : sprs::CsMat<usize> = matrix.to_csr();
    for (row_ind, row_vec) in m.outer_iterator().enumerate() {
        let mut s = 0;
        for (_col_ind, &val) in row_vec.iter() {
            s += val;
        }
        debug!("{} - {} molecules", rownames[row_ind], s);
        writeln!(summary_file, "{}\t{}", rownames[row_ind], s)?;

    }
    Ok(())
}

/* Validate Input/Output Files/Paths */

pub fn check_inputs_exist(bam_file: &str, cell_barcodes: &str, out_dir: &str) 
                          -> Result<(), Error> {
    for path in [bam_file, cell_barcodes].iter() {
        if !Path::new(&path).exists() {
            return Err(format_err!("Input {:?} does not exist", path));
        }
    }
    // check for BAM/CRAM index
    let extension = Path::new(bam_file).extension().unwrap().to_str().unwrap();
    match extension {
        "bam" => {
            let bai = bam_file.to_owned() + ".bai";
            if !Path::new(&bai).exists() {
                return Err(format_err!("BAM index {} does not exist", bai));
            }
        }
        "cram" => {
            let crai = bam_file.to_owned() + ".crai";
            if !Path::new(&crai).exists() {
                return Err(format_err!("CRAM index {} does not exist", crai));
            }
        }
        &_ => {
            return Err(format_err!("BAM file did not end in .bam or .cram. Unable to validate"));
        }
    }
    if !Path::new(&out_dir).exists() {
        match create_dir_all(&out_dir) {
            Err(_e) => { return Err(format_err!("Couldn't create results directory at {}", out_dir)); }
            _ => { info!("Created output directory at {}", out_dir);}
        }
    } else { return Err(format_err!("Specified output directory {} already exists", out_dir)); }
    Ok(())
}

pub fn check_inputs_exist_pseudoaligner(path: &str) -> Result<(), Error> {
    if !Path::new(&path).exists() {
        match create_dir_all(&path) {
            Err(_e) => { return Err(format_err!("Couldn't create temp directory at {}", path)); }
            _ => { info!("Created temp directory at {}", path);}
        }
    } else { return Err(format_err!("Specified tempdir {} already exists", path)); }
    Ok(()) 
}

pub fn check_inputs_exist_hla_db(path: &str) -> Result<(), Error> {
    if !Path::new(&path).exists() {
        return Err(format_err!("IMGT-HLA database directory {} does not exist", path));
    }
    for file in [ "hla_gen.fasta", "hla_nuc.fasta", "Allele_status.txt" ].iter() {
        let p : PathBuf = [path,file].iter().collect();
        if !p.exists() {
            return Err(format_err!("IMGT-HLA database file {} does not exist", file));
        }
    }
    Ok(())
}
    
pub fn check_inputs_exist_fasta(fasta_cds: &PathBuf, fasta_gen: &PathBuf)
                                 -> Result<(), Error> {
    for path in [fasta_cds, fasta_gen].iter() {
        if !Path::new(path).exists() {
            return Err(format_err!("Input file {:?} does not exist", path));
        }
    }
    Ok(())
}

pub fn check_inputs_exist_hla_idx(path: &str) -> Result<(), Error> {
    if !Path::new(&path).exists() {
        return Err(format_err!("Pseudoalignment index {} does not exist. Omit parameter -i to generate automatically.", path));
    }
    Ok(())
}

/* Helper Functions */

pub fn load_barcodes(filename: impl AsRef<Path>) -> Result<HashMap<Barcode, u32>, Error> {
    let r = File::open(filename.as_ref())?;
    let reader = BufReader::with_capacity(32 * 1024, r);

    let mut bc_set = HashMap::new();

    for (i, l) in reader.lines().enumerate() {
        let cb = Barcode::from_slice(l?.as_bytes());
        bc_set.insert(cb, i as u32);
    }
    let num_bcs = bc_set.len();
    if num_bcs == 0 {
        return Err(format_err!("Loaded 0 barcodes. Is your barcode file gzipped or empty?"));
    }
    debug!("Loaded {} barcodes", num_bcs);
    Ok(bc_set)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;
    use sprs::io::read_matrix_market;

    #[test]
    fn test_allele_fasta() {
        //This test takes ~1min in compliation "debug" mode
        //21 alignments x 6 allele sequences
        let mut cmds = Vec::new();
        let tmp_dir = tempdir().unwrap();
        let pl_dir = tmp_dir.path().join("p");
        let pl_dir = pl_dir.to_str().unwrap();
        let result_dir = tmp_dir.path().join("r");
        let out_file = result_dir.join("count_matrix.mtx");
        let out_file = out_file.to_str().unwrap();
        let result_dir = result_dir.to_str().unwrap();
        for l in &["scHLAcount", 
                   "-b", "test/test.bam",
                   "-g", "test/genomic_ABC.fa", 
                   "-f", "test/cds_ABC.fa",
                   "-c", "test/barcodes1.tsv",
                   "-o", result_dir,
                   "--pl-tmp", pl_dir ] {

            cmds.push(l.to_string());
        }
        let res = _main(cmds);
        assert!(!res.is_err());
        let seen_mat: TriMat<usize> = read_matrix_market(out_file).unwrap();
        let expected_mat: TriMat<usize> = read_matrix_market("test/test_allele_fasta.mtx").unwrap();
        assert_eq!(seen_mat.to_csr(), expected_mat.to_csr());
    }
    
    #[test]
    fn test_call() {
        let mut cmds = Vec::new();
        let tmp_dir = tempdir().unwrap();
        let pl_dir = tmp_dir.path().join("p");
        let pl_dir = pl_dir.to_str().unwrap();
        let result_dir = tmp_dir.path().join("r");
        let out_file = result_dir.join("count_matrix.mtx");
        let out_file = out_file.to_str().unwrap();
        let result_dir = result_dir.to_str().unwrap();
        for l in &["scHLAcount", 
                   "-b", "test/test.bam",
                   "-d", "test/fake_db",
                   "-i", "test/fake_db/hla_nuc.fasta.idx", 
                   "-c", "test/barcodes0.tsv",
                   "-o", result_dir,
                   "--pl-tmp", pl_dir ] {
            cmds.push(l.to_string());
        }
        let res = _main(cmds);
        assert!(!res.is_err());
        let seen_mat: TriMat<usize> = read_matrix_market(out_file).unwrap();
        let expected_mat: TriMat<usize> = read_matrix_market("test/test_call.mtx").unwrap();
        assert_eq!(seen_mat.to_csr(), expected_mat.to_csr());
    }
}
