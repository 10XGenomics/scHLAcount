// Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

use bio::alignment::pairwise::banded::Aligner;
use bio::alignment::pairwise::Scoring;
use bio::alignment::sparse::*;
use bio::alphabets::dna;
use bio::io::fasta;
use debruijn::dna_string::DnaString;
use debruijn_mapping::config::{KmerType, READ_COVERAGE_THRESHOLD};
use debruijn_mapping::pseudoaligner::Pseudoaligner;
use debruijn_mapping::utils::{read_obj, write_obj};
use failure::Error;
use itertools::Itertools;
use rust_htslib::bam::{IndexedReader, Read, Record};

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use serde::{Serialize, Deserialize};
use log::{debug, info};

use crate::config::{
    Barcode, EqClass, Umi, ALLELE_CONSENSUS_THRESHOLD, GENE_CONSENSUS_THRESHOLD, MIN_SCORE_CALL,
    MIN_SCORE_COUNT_ALIGNMENT, MIN_SCORE_COUNT_PSEUDO, PROC_BC_SEQ_TAG, PROC_UMI_SEQ_TAG,
};
use crate::hla::{read_hla_cds, read_hla_cds_string, Allele, AlleleParser};
use crate::locus::Locus;

const K: usize = 6; // kmer match length
const W: usize = 20; // Window size for creating the band
const MATCH: i32 = 1; // Match score
const MISMATCH: i32 = -1; // Mismatch score
const GAP_OPEN: i32 = -3; // Gap open score
const GAP_EXTEND: i32 = -1; // Gap extend score

/* Structs */

pub struct BamSeqReader {
    reader: IndexedReader,
    tmp_record: Record,
}

impl BamSeqReader {
    pub fn new(reader: IndexedReader) -> BamSeqReader {
        BamSeqReader {
            reader,
            tmp_record: Record::new(),
        }
    }

    pub fn fetch(&mut self, locus: &Locus) {
        let tid = self.reader.header().tid(locus.chrom.as_bytes()).unwrap();
        self.reader.fetch(tid, locus.start, locus.end).unwrap();
    }

    pub fn fetch_str(&mut self, region: &[u8]) {
        self.reader.fetch_str(&region).unwrap();
    }
}

#[derive(Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct BamCrRead {
    pub key: BamCrReadKey,
    pub sequence: DnaString,
    pub primary: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct BamCrReadKey {
    pub barcode: Barcode,
    pub umi: Umi,
}

#[derive(Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct BusCount {
    barcode_id: u32,
    umi_id: u32,
    eq_class_id: u32,
}

#[derive(Serialize, Deserialize)]
pub struct EqClassDb {
    barcodes: HashMap<Barcode, u32>,
    umis: HashMap<Umi, u32>,
    pub eq_classes: HashMap<EqClass, u32>,
    counts: Vec<BusCount>,
}

#[derive(Debug)]
pub struct Scores {
    pub cell_index: u32,
    pub umi: Umi,
    pub max_score: i32,
    pub max_gene: Vec<u8>,
    pub max_allele: usize,
}

#[derive(Debug)]
pub struct MatrixEntry {
    pub row: usize,
    pub column: usize,
    pub value: usize,
}

#[derive(Clone, Copy)]
pub struct Metrics {
    pub num_reads: usize,
    pub num_non_primary: usize,
    pub num_not_cell_bc: usize,
    pub num_cds_align: usize,
    pub num_gen_align: usize,
    pub num_not_aligned: usize,
}

impl<'a> EqClassDb {
    pub fn new(_nitems: usize) -> EqClassDb {
        EqClassDb {
            barcodes: HashMap::new(),
            umis: HashMap::new(),
            eq_classes: HashMap::new(),
            counts: Vec::new(),
        }
    }

    pub fn count(&mut self, barcode: &Barcode, umi: &Umi, eq_class: &EqClass) {
        let barcode_id = match self.barcodes.get(barcode).cloned() {
            Some(id) => id,
            None => {
                let id = self.barcodes.len() as u32;
                self.barcodes.insert(barcode.clone(), id);
                id
            }
        };

        let umi_id = match self.umis.get(umi).cloned() {
            Some(id) => id,
            None => {
                let id = self.umis.len() as u32;
                self.umis.insert(umi.clone(), id);
                id
            }
        };

        let eq_class_id = match self.eq_classes.get(eq_class).cloned() {
            Some(id) => id,
            None => {
                let id = self.eq_classes.len() as u32;
                self.eq_classes.insert(eq_class.clone(), id);
                id
            }
        };

        let count = BusCount {
            barcode_id,
            umi_id,
            eq_class_id,
        };

        self.counts.push(count);
    }

    pub fn sort(&mut self) {
        self.counts.sort();
    }

    pub fn eq_class_counts(&mut self, nitems: usize) -> (crate::em::EqClassCounts, Vec<usize>) {
        let mut rev_map = HashMap::<u32, &EqClass>::new();
        let mut counts_reads: HashMap<EqClass, u32> = HashMap::new();
        let mut counts_umi: HashMap<EqClass, u32> = HashMap::new();
        let mut reads_explained = vec![0; nitems];
        self.sort();

        for (cls, id) in &self.eq_classes {
            rev_map.insert(*id, cls);
        }

        let mut uniq = 0;
        let mut total_reads = 0;
        let mut buf = Vec::new();

        for ((_bc, _umi), hits) in &self.counts.iter().group_by(|c| (c.barcode_id, c.umi_id)) {
            buf.clear();
            let mut flag: bool = false;
            for c in hits {
                if rev_map.contains_key(&c.eq_class_id) {
                    if flag {
                        intersect(&mut buf, rev_map[&c.eq_class_id].as_ref());
                    } else {
                        buf.extend(rev_map[&c.eq_class_id].as_ref());
                        flag = true;
                    }
                    total_reads += 1;
                    for t in rev_map[&c.eq_class_id].as_ref() {
                        reads_explained[*t as usize] += 1;
                    }
                    let eqclass = EqClass::from_slice(rev_map[&c.eq_class_id].as_ref());
                    let count = counts_reads.entry(eqclass).or_default();
                    *count += 1;
                }
            }
            if flag {
                buf.sort();
                buf.dedup();
                let eqclass = EqClass::from_slice(&buf);
                let count = counts_umi.entry(eqclass).or_default();
                *count += 1;
                uniq += 1;
            }
        }

        debug!(
            "mean reads per UMI: {}",
            f64::from(total_reads) / f64::from(uniq)
        );

        (
            crate::em::EqClassCounts {
                nitems,
                counts_reads,
                counts_umi,
                nreads: total_reads,
            },
            reads_explained,
        )
    }
}

impl Iterator for BamSeqReader {
    type Item = Result<BamCrRead, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let r = self.reader.read(&mut self.tmp_record);
            let mut p: bool = true;
            if let Err(e) = r {
                if e.is_eof() {
                    return None;
                } else {
                    return Some(Err(e.into()));
                }
            };
            if self.tmp_record.is_secondary() || self.tmp_record.is_supplementary() {
                p = false;
            }

            // Get original read sequence from record.
            let sequence = DnaString::from_acgt_bytes(&self.tmp_record.seq().as_bytes());

            let barcode = match self
                .tmp_record
                .aux(PROC_BC_SEQ_TAG)
                .map(|x| Barcode::from_slice(x.string()))
            {
                Some(bc) => bc,
                None => continue,
            };

            let umi = match self
                .tmp_record
                .aux(PROC_UMI_SEQ_TAG)
                .map(|x| Umi::from_slice(x.string()))
            {
                Some(bc) => bc,
                None => continue,
            };

            return Some(Ok(BamCrRead {
                sequence,
                key: BamCrReadKey { barcode, umi },
                primary: p,
            }));
        }
    }
}

/// Compute the intersection of v1 and v2 inplace on top of v1
/// v1 and v2 must be sorted
fn intersect<T: Eq + Ord>(v1: &mut Vec<T>, v2: &[T]) {
    if v1.is_empty() {
        return;
    }

    if v2.is_empty() {
        v1.clear();
    }

    let mut fill_idx1 = 0;
    let mut idx1 = 0;
    let mut idx2 = 0;

    while idx1 < v1.len() && idx2 < v2.len() {
        match v1[idx1].cmp(&v2[idx2]) {
            Ordering::Less => idx1 += 1,
            Ordering::Greater => idx2 += 1,
            Ordering::Equal => {
                v1.swap(fill_idx1, idx1);
                idx1 += 1;
                idx2 += 1;
                fill_idx1 += 1;
            }
        }
    }
}

pub fn mapping_wrapper(
    hla_index: PathBuf,
    outdir: PathBuf,
    bam: PathBuf,
    locus: Option<&Locus>,
    use_unmapped: bool,
) -> Result<PathBuf, Error> {
    info!("Reading database index from disk");
    let index: Pseudoaligner<KmerType> = read_obj(&hla_index)?;
    info!("Finished reading index!");
    info!("Mapping reads from BAM to database");
    let mut itr = BamSeqReader::new(IndexedReader::from_path(&bam)?);

    if let Some(l) = locus {
        debug!("Locus: {:?}", l);
        itr.fetch(l);
    }

    let mut hits_file = outdir.clone();
    hits_file.set_extension("counts.bin");

    let mut rid = 0;
    let mut some_aln = 0;
    let mut long_aln = 0;

    let mut eq_counts = EqClassDb::new(index.tx_names.len());

    for _rec in itr {
        rid += 1;
        let rec = _rec?;

        if let Some((eq_class, cov)) = index.map_read(&rec.sequence) {
            some_aln += 1;
            if cov > MIN_SCORE_CALL {
                long_aln += 1;
                eq_counts.count(
                    &rec.key.barcode,
                    &rec.key.umi,
                    &EqClass::from_slice(&eq_class),
                );
            }
        } else if let Some((eq_class, cov)) = index.map_read(&DnaString::from_acgt_bytes(
            dna::revcomp(&rec.sequence.to_ascii_vec()).as_slice(),
        )) {
            some_aln += 1;
            if cov > MIN_SCORE_CALL {
                long_aln += 1;
                eq_counts.count(
                    &rec.key.barcode,
                    &rec.key.umi,
                    &EqClass::from_slice(&eq_class),
                );
            }
        }
        if rid % 100_000 == 0 {
            info!("analyzed {} reads. Mapped {} with score at least {}, used {} with score at least {}", 
                rid, some_aln, READ_COVERAGE_THRESHOLD, long_aln, MIN_SCORE_CALL);
        }
    }

    info!(
        "analyzed {} reads. Mapped {} with score at least {}, used {} with score at least {}",
        rid, some_aln, READ_COVERAGE_THRESHOLD, long_aln, MIN_SCORE_CALL
    );
    if use_unmapped {
        let mut rid = 0;
        let mut some_aln = 0;
        let mut long_aln = 0;
        let mut itr = BamSeqReader::new(IndexedReader::from_path(&bam)?);
        itr.fetch_str(b"*");
        for _rec in itr {
            rid += 1;
            let rec = _rec?;
            if let Some((eq_class, cov)) = index.map_read(&rec.sequence) {
                some_aln += 1;
                if cov > MIN_SCORE_CALL {
                    long_aln += 1;
                    eq_counts.count(
                        &rec.key.barcode,
                        &rec.key.umi,
                        &EqClass::from_slice(&eq_class),
                    );
                }
            }

            if rid % 100_000 == 0 {
                info!("analyzed {} unaligned reads. Mapped {} with score at least {}, used {} with score at least {}", 
                    rid, some_aln, READ_COVERAGE_THRESHOLD, long_aln, MIN_SCORE_CALL);
            }
        }
        info!("analyzed {} unaligned reads. Mapped {} with score at least {}, used {} with score at least {}", 
                    rid, some_aln, READ_COVERAGE_THRESHOLD, long_aln, MIN_SCORE_CALL);
    }
    write_obj(&eq_counts, &hits_file)?;
    Ok(hits_file)
}
#[allow(clippy::cognitive_complexity)]
pub fn map_and_count_sw(
    bam: PathBuf,
    out_dir: &str,
    barcodes: &HashMap<Barcode, u32>,
    locus: &Locus,
    cds: PathBuf,
    genomic: PathBuf,
    primary_only: bool,
) -> Result<(Vec<MatrixEntry>, usize, Metrics, Vec<String>), Error> {
    let fa: fasta::Reader<File> = fasta::Reader::from_file(genomic)?;
    let (seqs_gen, tx_names_gen) = read_hla_cds_string(fa, HashSet::new(), false)?;
    let mut genomic_kmer_hashes: Vec<(&String, HashMapFx<&[u8], Vec<u32>>, &[u8])> = Vec::new();
    for (i, s) in seqs_gen.iter().enumerate() {
        genomic_kmer_hashes.push((&tx_names_gen[i], hash_kmers(s, K), s));
    }
    info!("Indexed genome sequences for SW alignment");

    let fa: fasta::Reader<File> = fasta::Reader::from_file(cds)?;
    let (seqs_cds, tx_names_cds) = read_hla_cds_string(fa, HashSet::new(), false)?;
    let mut cds_kmer_hashes: Vec<(&String, HashMapFx<&[u8], Vec<u32>>, &[u8])> = Vec::new();
    for (i, s) in seqs_cds.iter().enumerate() {
        cds_kmer_hashes.push((&tx_names_cds[i], hash_kmers(s, K), s));
    }
    info!("Indexed CDS sequences for SW alignment");

    // Sequences must be in the same order in CDS and GEN.
    // Have to assume they have the same set of names, though, for the next steps to work!
    genomic_kmer_hashes.sort_by_key(|x| x.0);
    cds_kmer_hashes.sort_by_key(|x| x.0);

    //Local alignment
    let scoring = Scoring {
        gap_open: GAP_OPEN,
        gap_extend: GAP_EXTEND,
        match_fn: |a: u8, b: u8| if a == b { MATCH } else { MISMATCH },
        match_scores: Some((MATCH, MISMATCH)),
        xclip_prefix: 0,
        xclip_suffix: 0,
        yclip_prefix: 0,
        yclip_suffix: 0,
    };
    let mut aligner = Aligner::with_scoring(scoring, K, W);

    let mut itr = BamSeqReader::new(IndexedReader::from_path(bam)?);
    itr.fetch(locus);

    // what equivalence class corresponds to a gene or an allele?
    let allele_parser = AlleleParser::new();
    let mut alleles_cds: Vec<(Allele, u32)> = Vec::new();
    for (i, t) in tx_names_cds.iter().enumerate() {
        if let Ok(a) = allele_parser.parse(t) {
            alleles_cds.push((a, i as u32));
        }
    }

    // Mapping of genes/alleles to row numbers in the matrix
    let d: PathBuf = [out_dir, "labels.tsv"].iter().collect();
    let mut labels_file = BufWriter::new(File::create(d).unwrap());
    let mut eq_class_to_gene: HashMap<Vec<u32>, (Vec<u8>, usize)> = HashMap::new();
    let mut genes_to_rownums: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut nrows: usize = 0;
    let mut rownames: Vec<String> = Vec::new();
    for (g, mut a_group) in &alleles_cds
        .iter()
        .sorted_by_key(|a| &a.0.gene)
        .group_by(|a| &a.0.gene)
    {
        genes_to_rownums.insert(g.clone(), nrows);
        let a1 = a_group.next().unwrap();
        eq_class_to_gene.insert(vec![a1.1], (g.clone(), 0));
        if let Some(a2) = a_group.next() {
            //two alleles : three rows
            eq_class_to_gene.insert(vec![a2.1], (g.clone(), 1));
            let mut x = vec![a1.1, a2.1];
            x.sort();
            eq_class_to_gene.insert(x.to_owned(), (g.clone(), 2));
            nrows += 3;

            writeln!(
                labels_file,
                "{}",
                String::from_utf8(a1.0.name.clone()).unwrap()
            )?;
            writeln!(
                labels_file,
                "{}",
                String::from_utf8(a2.0.name.clone()).unwrap()
            )?;
            writeln!(
                labels_file,
                "{}",
                String::from_utf8(g.clone().to_vec()).unwrap()
            )?;
            rownames.push(String::from_utf8(a1.0.name.clone()).unwrap());
            rownames.push(String::from_utf8(a2.0.name.clone()).unwrap());
            rownames.push(String::from_utf8(g.clone().to_vec()).unwrap());
        } else {
            // only one allele : one row
            nrows += 1;
            writeln!(
                labels_file,
                "{}",
                String::from_utf8(a1.0.name.clone()).unwrap()
            )?;
            rownames.push(String::from_utf8(a1.0.name.clone()).unwrap());
        }
    }
    labels_file.flush()?;

    let mut metrics: Metrics = Metrics {
        num_reads: 0,
        num_non_primary: 0,
        num_not_cell_bc: 0,
        num_not_aligned: 0,
        num_cds_align: 0,
        num_gen_align: 0,
    };

    let mut scores: Vec<Scores> = Vec::new();
    for _rec in itr {
        let rec = _rec?;
        metrics.num_reads += 1;
        if metrics.num_reads % 100_000 == 0 {
            info!(
                "analyzed {} reads. Mapped {} with score at least {}",
                metrics.num_reads,
                metrics.num_gen_align + metrics.num_cds_align,
                MIN_SCORE_COUNT_ALIGNMENT
            );
        }

        if !rec.primary {
            metrics.num_non_primary += 1;
            if primary_only {
                continue;
            }
        }

        // iterator only returns reads with CB and UB but it might not be in the list
        if !barcodes.contains_key(&rec.key.barcode) {
            metrics.num_not_cell_bc += 1;
            continue;
        }

        let mut max_score: i32 = 0;
        let mut cls: Vec<u32> = Vec::new();

        let seq = rec.sequence;
        let read = seq.to_ascii_vec();
        let read = read.as_slice();
        let read_rev = dna::revcomp(seq.to_ascii_vec());
        let read_rev = read_rev.as_slice();

        //align to cds sequences
        for (i, h) in cds_kmer_hashes.iter().enumerate() {
            //align fwd read
            let aln = aligner.custom_with_prehash(read, h.2, &h.1);
            if aln.score > max_score {
                cls = vec![i as u32];
                max_score = aln.score;
            } else if aln.score == max_score {
                cls.push(i as u32);
            }
            //align reverse read
            let aln = aligner.custom_with_prehash(read_rev, h.2, &h.1);
            if aln.score > max_score {
                cls = vec![i as u32];
                max_score = aln.score;
            } else if aln.score == max_score {
                cls.push(i as u32);
            }
        }

        //if necessary, align to genomic sequences
        if max_score < MIN_SCORE_COUNT_ALIGNMENT as i32 {
            cls = Vec::new();
            for (i, h) in genomic_kmer_hashes.iter().enumerate() {
                //align fwd read
                let aln = aligner.custom_with_prehash(read, h.2, &h.1);
                if aln.score > max_score {
                    cls = vec![i as u32];
                    max_score = aln.score;
                } else if aln.score == max_score {
                    cls.push(i as u32);
                }
                //align reverse read
                let aln = aligner.custom_with_prehash(read_rev, h.2, &h.1);
                if aln.score > max_score {
                    cls = vec![i as u32];
                    max_score = aln.score;
                } else if aln.score == max_score {
                    cls.push(i as u32);
                }
            }
            if max_score >= MIN_SCORE_COUNT_ALIGNMENT as i32 {
                metrics.num_gen_align += 1;
            }
        } else {
            metrics.num_cds_align += 1;
        }
        //println!("{}",max_score);
        if max_score >= MIN_SCORE_COUNT_ALIGNMENT as i32 {
            if let Some((max_gene, max_allele)) = eq_class_to_gene.get(&cls) {
                let s = Scores {
                    cell_index: *barcodes.get(&rec.key.barcode).unwrap(),
                    umi: rec.key.umi,
                    max_score,
                    max_gene: max_gene.clone(),
                    max_allele: *max_allele,
                };
                //debug!("{:?}",&s);
                scores.push(s);
            }
        } else {
            metrics.num_not_aligned += 1;
        }
    }
    info!(
        "analyzed {} reads. Mapped {} with score at least {}",
        metrics.num_reads,
        metrics.num_gen_align + metrics.num_cds_align,
        MIN_SCORE_COUNT_ALIGNMENT
    );
    info!(
        "{} reads aligned to a single-gene equivalence class",
        scores.len()
    );

    let entries = count(&scores, genes_to_rownums, nrows);
    Ok((entries, nrows, metrics, rownames))
}

#[allow(clippy::cognitive_complexity)]
pub fn map_and_count_pseudo(
    bam: PathBuf,
    out_dir: &str,
    barcodes: &HashMap<Barcode, u32>,
    locus: &Locus,
    cds: PathBuf,
    genomic: PathBuf,
    primary_only: bool,
) -> Result<(Vec<MatrixEntry>, usize, Metrics, Vec<String>), Error> {
    let fa: fasta::Reader<File> = fasta::Reader::from_file(genomic)?;
    let (seqs, tx_names) = read_hla_cds(fa, HashSet::new(), false)?;
    let tx_gene_map = HashMap::new();
    let genomic_index = debruijn_mapping::build_index::build_index::<
        debruijn_mapping::config::KmerType,
    >(&seqs, &tx_names, &tx_gene_map)?;
    info!("Built genome index for pseudoalignment");
    let fa: fasta::Reader<File> = fasta::Reader::from_file(cds)?;
    let (seqs, tx_names) = read_hla_cds(fa, HashSet::new(), false)?;
    let cds_index = debruijn_mapping::build_index::build_index::<debruijn_mapping::config::KmerType>(
        &seqs,
        &tx_names,
        &tx_gene_map,
    )?;
    info!("Built CDS index for pseudoalignment");

    let mut itr = BamSeqReader::new(IndexedReader::from_path(bam)?);
    itr.fetch(locus);

    // what equivalence class corresponds to a gene or an allele?
    let allele_parser = AlleleParser::new();
    let mut alleles: Vec<(Allele, u32)> = Vec::new();
    for (i, t) in tx_names.iter().enumerate() {
        if let Ok(a) = allele_parser.parse(t) {
            alleles.push((a, i as u32));
        }
    }
    // Mapping of genes/alleles to row numbers in the matrix
    let d: PathBuf = [out_dir, "labels.tsv"].iter().collect();
    let mut labels_file = BufWriter::new(File::create(d).unwrap());
    let mut eq_class_to_gene: HashMap<Vec<u32>, (Vec<u8>, usize)> = HashMap::new();
    let mut genes_to_rownums: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut nrows: usize = 0;
    let mut rownames: Vec<String> = Vec::new();
    for (g, mut a_group) in &alleles
        .iter()
        .sorted_by_key(|a| &a.0.gene)
        .group_by(|a| &a.0.gene)
    {
        genes_to_rownums.insert(g.clone(), nrows);
        let a1 = a_group.next().unwrap();
        eq_class_to_gene.insert(vec![a1.1], (g.clone(), 0));
        if let Some(a2) = a_group.next() {
            //two alleles : three rows
            eq_class_to_gene.insert(vec![a2.1], (g.clone(), 1));
            let mut x = vec![a1.1, a2.1];
            x.sort();
            eq_class_to_gene.insert(x.to_owned(), (g.clone(), 2));
            nrows += 3;

            writeln!(
                labels_file,
                "{}",
                String::from_utf8(a1.0.name.clone()).unwrap()
            )?;
            writeln!(
                labels_file,
                "{}",
                String::from_utf8(a2.0.name.clone()).unwrap()
            )?;
            writeln!(
                labels_file,
                "{}",
                String::from_utf8(g.clone().to_vec()).unwrap()
            )?;
            rownames.push(String::from_utf8(a1.0.name.clone()).unwrap());
            rownames.push(String::from_utf8(a2.0.name.clone()).unwrap());
            rownames.push(String::from_utf8(g.clone().to_vec()).unwrap());
        } else {
            // only one allele : one row
            nrows += 1;
            writeln!(
                labels_file,
                "{}",
                String::from_utf8(a1.0.name.clone()).unwrap()
            )?;
            rownames.push(String::from_utf8(a1.0.name.clone()).unwrap());
        }
    }
    labels_file.flush()?;

    let mut metrics: Metrics = Metrics {
        num_reads: 0,
        num_non_primary: 0,
        num_not_cell_bc: 0,
        num_not_aligned: 0,
        num_cds_align: 0,
        num_gen_align: 0,
    };

    let mut scores: Vec<Scores> = Vec::new();
    for _rec in itr {
        let rec = _rec?;
        metrics.num_reads += 1;
        if metrics.num_reads % 100_000 == 0 {
            info!(
                "analyzed {} reads. Mapped {} with score at least {}",
                metrics.num_reads,
                metrics.num_gen_align + metrics.num_cds_align,
                MIN_SCORE_COUNT_PSEUDO
            );
        }

        if !rec.primary {
            metrics.num_non_primary += 1;
            if primary_only {
                continue;
            }
        }

        // iterator only returns reads with CB and UB but it might not be in the list
        if !barcodes.contains_key(&rec.key.barcode) {
            metrics.num_not_cell_bc += 1;
            continue;
        }

        let mut aln = cds_index.map_read(&rec.sequence);
        let gen_aln = genomic_index.map_read(&rec.sequence);
        let seq_rev: DnaString =
            DnaString::from_acgt_bytes(dna::revcomp(&rec.sequence.to_ascii_vec()).as_slice());
        let cds_rev_aln = cds_index.map_read(&seq_rev);
        let gen_rev_aln = genomic_index.map_read(&seq_rev);

        if let Some((_, cov)) = aln {
            if cov > MIN_SCORE_COUNT_PSEUDO {
                metrics.num_cds_align += 1;
            }
        } else if let Some((_, cov)) = cds_rev_aln {
            if cov > MIN_SCORE_COUNT_PSEUDO {
                metrics.num_cds_align += 1;
                aln = cds_rev_aln;
            }
        } else if let Some((_, cov)) = gen_aln {
            if cov > MIN_SCORE_COUNT_PSEUDO {
                metrics.num_gen_align += 1;
                aln = gen_aln;
            }
        } else if let Some((_, cov)) = gen_rev_aln {
            if cov > MIN_SCORE_COUNT_PSEUDO {
                metrics.num_gen_align += 1;
                aln = gen_rev_aln;
            }
        } else {
            metrics.num_not_aligned += 1;
        }

        if let Some((cls, cov)) = aln {
            if let Some((max_gene, max_allele)) = eq_class_to_gene.get(&cls) {
                let s = Scores {
                    cell_index: *barcodes.get(&rec.key.barcode).unwrap(),
                    umi: rec.key.umi,
                    max_score: cov as i32,
                    max_gene: max_gene.clone(),
                    max_allele: *max_allele,
                };
                scores.push(s);
            }
        }
    }
    info!(
        "analyzed {} reads. Mapped {} with score at least {}",
        metrics.num_reads,
        metrics.num_gen_align + metrics.num_cds_align,
        MIN_SCORE_COUNT_PSEUDO
    );
    info!(
        "{} reads aligned to a single-gene equivalence class",
        scores.len()
    );

    let entries = count(&scores, genes_to_rownums, nrows);
    Ok((entries, nrows, metrics, rownames))
}

pub fn count(
    scores: &[Scores],
    genes_to_rownums: HashMap<Vec<u8>, usize>,
    nrows: usize,
) -> Vec<MatrixEntry> {
    let mut entries: Vec<MatrixEntry> = Vec::new();
    for (cell_index, cell_scores) in &scores
        .iter()
        .sorted_by_key(|s| &s.cell_index)
        .group_by(|s| &s.cell_index)
    {
        let mut matrix_row: Vec<usize> = vec![0; nrows];
        let mut parsed_scores: HashMap<&Umi, Vec<&Scores>> = HashMap::new();
        for score in cell_scores.into_iter() {
            parsed_scores
                .entry(&score.umi)
                .or_insert_with(Vec::new)
                .push(score);
        }
        for (_umi, all_scores) in parsed_scores.into_iter() {
            let mut nreads: f64 = 0.0;
            let mut gene_frq: HashMap<Vec<u8>, (usize, usize, usize)> = HashMap::new();
            for score in all_scores.into_iter() {
                match score.max_allele {
                    0 => {
                        gene_frq
                            .entry(score.max_gene.clone())
                            .or_insert((0, 0, 0))
                            .0 += 1
                    } //allele 1
                    1 => {
                        gene_frq
                            .entry(score.max_gene.clone())
                            .or_insert((0, 0, 0))
                            .1 += 1
                    } //allele 2
                    2 => {
                        gene_frq
                            .entry(score.max_gene.clone())
                            .or_insert((0, 0, 0))
                            .2 += 1
                    } //gene
                    _ => (),
                };
                nreads += 1.0;
            }
            let mut max_gene_this_umi: Vec<u8> = Vec::new();
            let mut max_allele_this_umi: usize = 0;
            let mut max_reads_this_umi: usize = 0;
            for (g, frq) in gene_frq.into_iter() {
                let nreads_g = frq.0 + frq.1 + frq.2;
                if nreads_g > max_reads_this_umi
                    && nreads_g as f64 > GENE_CONSENSUS_THRESHOLD * nreads
                {
                    max_gene_this_umi = g;
                    max_reads_this_umi = nreads_g;
                    if (frq.0 as f64 > ALLELE_CONSENSUS_THRESHOLD * nreads)
                        && (ALLELE_CONSENSUS_THRESHOLD * nreads > frq.1 as f64)
                    {
                        max_allele_this_umi = 0; //allele 1
                    } else if (frq.1 as f64 > ALLELE_CONSENSUS_THRESHOLD * nreads)
                        && (ALLELE_CONSENSUS_THRESHOLD * nreads > frq.0 as f64)
                    {
                        max_allele_this_umi = 1; //allele 2
                    } else {
                        max_allele_this_umi = 2; //gene
                    }
                }
            }
            if !max_gene_this_umi.is_empty() {
                matrix_row
                    [max_allele_this_umi + genes_to_rownums.get(&max_gene_this_umi).unwrap()] += 1;
            }
        }
        for (c, val) in matrix_row.iter().enumerate() {
            if *val > 0 {
                let e = MatrixEntry {
                    row: c,
                    column: *cell_index as usize,
                    value: *val,
                };
                entries.push(e);
            }
        }
    }
    entries
}
