// Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

/* Constants */
pub const MIN_SCORE_CALL: usize = 40;
pub const MIN_SCORE_COUNT_PSEUDO: usize = 60;
pub const MIN_SCORE_COUNT_ALIGNMENT: usize = 60;
pub const GENE_CONSENSUS_THRESHOLD: f64 = 0.5;
pub const ALLELE_CONSENSUS_THRESHOLD: f64 = 0.1;

pub const EM_ITERS: usize = 2000;
pub const EM_REL_TH: f64 = 5e-4;
pub const EM_ABS_TH: f64 = 5e-3;
pub const EM_CARE_TH: f64 = 1e-5;
pub const MIN_READS_CALL: usize = 100;
pub const HOMOZYGOUS_TH: f64 = 0.15;

pub const PROC_BC_SEQ_TAG: &[u8] = b"CB";
pub const PROC_UMI_SEQ_TAG: &[u8] = b"UB";

pub const PAIRS_TO_OUTPUT: usize = 5;
pub const WEIGHTS_TO_OUTPUT: usize = 10;

/* Types */
use smallvec::SmallVec;
pub type Barcode = SmallVec<[u8; 24]>;
pub type Umi = SmallVec<[u8; 16]>;
pub type EqClass = SmallVec<[u32; 4]>;
