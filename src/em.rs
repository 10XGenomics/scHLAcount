use debruijn_mapping::{utils, config, pseudoaligner};
use failure::Error;
use itertools::Itertools;

use std::collections::{HashMap,HashSet};
use std::fs::File;
use std::io::{Write,BufWriter};
use std::path::{PathBuf};

use crate::mapping::EqClassDb;
use crate::config::{EM_ITERS,EM_REL_TH,EM_ABS_TH,EM_CARE_TH,MIN_READS_CALL,EqClass, PAIRS_TO_OUTPUT, WEIGHTS_TO_OUTPUT};

#[derive(Default)]
pub struct EqClassCounts {
    pub nitems: usize,
    pub counts_reads: HashMap<EqClass, u32>,
    pub counts_umi: HashMap<EqClass, u32>,
    pub nreads: u32,
}

impl EqClassCounts {
/*    pub fn new() -> EqClassCounts {
        EqClassCounts {
            nitems: 0,
            counts_reads: HashMap::new(),
            counts_umi: HashMap::new(),
            nreads: 0,
        }
    }*/
    pub fn pair_reads_explained(&self, a: u32, b: u32) -> u32 {
        let mut exp = 0;
        for (cls, count) in self.counts_reads.iter() {
            if cls.contains(&a) || cls.contains(&b) {
                exp += *count;
            }
        }
        exp
    }
}

#[allow(clippy::needless_range_loop)]
impl EmProblem for EqClassCounts {
    fn init(&self) -> Vec<f64> {
        vec![1.0/(self.nitems as f64); self.nitems]
    }

    fn reg(&self, theta: &mut [f64]) {
        let mut norm = 0.0;

        for i in 0 .. theta.len() {

            // Clamp weight
            let mut v = theta[i];
            if v > 1.0 { v = 1.0 };
            if v < 0.0 { v = 1e-15 };
            theta[i] = v;
            norm += v;
        }

        let inv_norm = 1.0 / norm;
        for i in 0 .. theta.len() {
            theta[i] *= inv_norm;
        }
    }

    fn f(&self, theta1: &[f64], theta2: &mut [f64]) {
        let nitems = self.nitems;

        let mut total_counts = 0.0;

        for i in 0..theta2.len() {
            theta2[i] = 0.0;
        }

        for (class, count) in &self.counts_umi {

            let mut norm = 0.0;
            for tx in class {
                norm += theta1[*tx as usize];
            }

            for tx in class {
                let tx_count = theta1[*tx as usize] / norm * f64::from(*count);
                theta2[*tx as usize] += tx_count;
                total_counts += tx_count;
            }
        }

        let mut max_abs_diff = 0.0;
        let mut max_rel_diff = 0.0;

        for i in 0 .. nitems {
            let old_weights = theta1[i];
            let new_weights = theta2[i] / total_counts;

            let abs_diff = (old_weights - new_weights).abs();
            let rel_diff = abs_diff / old_weights;

            if abs_diff > max_abs_diff {
                max_abs_diff = abs_diff;
            }

            if new_weights > 1e-2 && rel_diff > max_rel_diff {
                max_rel_diff = rel_diff
            }

            theta2[i] = new_weights;
        }
    }

    fn likelihood(&self, theta: &[f64]) -> f64 {
        let mut ll = 0.0;

        for (class, count) in &self.counts_umi {

            // Exclude empty equivalence classes
            if class.is_empty() {
                continue;
            }

            let mut theta_sum = 0.0;
            for tx in class {
                theta_sum += theta[*tx as usize];
            }

            ll += f64::from(*count) * theta_sum.ln();
        }

        ll
    }
}

pub fn em(eq_classes: &EqClassCounts) -> Vec<f64> {
    let nitems = eq_classes.nitems;

    // initialize weights
    let mut weights = vec![1.0/(nitems as f64); nitems];
    
    let mut iters = 0;

    // Abundance required to 'care' about a relative change
    //let rel_care_thresh = 1e-3 / (nitems as f64);

    loop {

        let mut pseudocounts = vec![0.0; nitems];
        let mut total_counts = 0.0;

        for (class, count) in &eq_classes.counts_umi {

            let mut norm = 0.0;
            for tx in class {
                norm += weights[*tx as usize];
            }

            for tx in class {
                let tx_count = weights[*tx as usize] / norm * f64::from(*count);
                pseudocounts[*tx as usize] += tx_count;
                total_counts += tx_count;
            }
        }

        let mut max_abs_diff = 0.0;
        let mut max_rel_diff = 0.0;
        let mut simpsons = 0.0;

        for i in 0 .. nitems {
            let old_weights = weights[i];
            let new_weights = pseudocounts[i] / total_counts;

            let abs_diff = (old_weights - new_weights).abs();
            let rel_diff = abs_diff / old_weights;

            if abs_diff > max_abs_diff {
                max_abs_diff = abs_diff;
            }

            if new_weights > 1e-2 && rel_diff > max_rel_diff {
                max_rel_diff = rel_diff
            }

            weights[i] = new_weights;
            simpsons += new_weights*new_weights;
        }

        let ll = eq_classes.likelihood(&weights);
        debug!("iter: {}, ll: {}, div: {}, rel_diff: {}, abs_diff: {}", iters, ll, 1.0/simpsons, max_rel_diff, max_abs_diff);
        iters += 1;
        if (max_abs_diff < 0.00005 && max_rel_diff < 0.0001) || iters > 5000 {
            break;
        }

        
    }

    weights
}

/// Encapsulate an EM optimization problem so that it can run through an accelerated EM loop (SquareM).
pub trait EmProblem {

    // Make an initial estimate of the parameter vector. May be a naive estimate.
    fn init(&self) -> Vec<f64>;
    
    // Regularize a parameter vector -- fixup an inconsistencies in the parameters.
    // E.g. legalize values or enforce global constrains.
    fn reg(&self, theta: &mut[f64]);

    // Update the parameters -- one EM step
    fn f(&self, theta1: &[f64], theta_2: &mut [f64]);

    // Compute the likelihood of a parameter set
    fn likelihood(&self, theta: &[f64]) -> f64;
}



/// SquareM EM acceleration method. 
/// As described in:
/// Varadhan, Ravi, and Christophe Roland. 
/// "Simple and globally convergent methods for accelerating the convergence of any EM algorithm." 
/// Scandinavian Journal of Statistics 35.2 (2008): 335-353.
/// Takes an implementation of `EmProblem` and applies the accelerated EM algorithm.
pub fn squarem<T: EmProblem>(p: &T) -> Vec<f64> {


    // Array for holding theta
    let mut theta0 = p.init();
    let mut theta1 = theta0.clone();
    let mut theta2 = theta0.clone();
    let mut theta_sq = theta0.clone();

    let mut r = theta0.clone();
    let mut v = theta0.clone();
    
    let n = theta0.len();
    let mut iters = 0; 

    let rel_care_thresh = Some(EM_CARE_TH);

    loop {

        // Get theta1
        p.f(&theta0, &mut theta1);
        p.f(&theta1, &mut theta2);

        let mut rsq: f64 = 0.0;
        let mut vsq: f64 = 0.0;

        for i in 0..n {
            r[i] = theta1[i] - theta0[i];
            rsq += r[i].powi(2);

            v[i] = theta2[i] - theta1[i] - r[i];
            vsq += v[i].powi(2);
        }

        let mut alpha = -rsq.sqrt() / vsq.sqrt();
        let mut alpha_tries = 1;
        let mut lsq : f64;
        let mut l2 : f64;

        loop {
            
            let alpha_sq = alpha.powi(2);

            for i in 0..n {
                theta_sq[i] = theta0[i] - 2.0 * alpha * r[i] + alpha_sq * v[i]
            }

            p.reg(&mut theta_sq);

            lsq = p.likelihood(&theta_sq);
            l2 = p.likelihood(&theta2);

            if lsq > l2 || alpha_tries > 5 { 
                break;
            } else {
                alpha = (alpha + -1.0) / 2.0;
            }

            alpha_tries += 1;
        }


        let (max_rel_diff, max_abs_diff) = 
            if lsq > l2 {
                let diff = diffs(&theta0, &theta_sq, rel_care_thresh);
                std::mem::swap(&mut theta0, &mut theta_sq);
                diff
            } else {
                let diff = diffs(&theta0, &theta2, rel_care_thresh);
                std::mem::swap(&mut theta0, &mut theta2);
                diff
            };


        debug!("iter: {}, ll2: {}, llsq: {}, alpha_tries: {}, rel_diff: {}, abs_diff: {}", iters, l2, lsq, alpha_tries, max_rel_diff, max_abs_diff);
        iters += 1;

        if (max_abs_diff <EM_ABS_TH && max_rel_diff < EM_REL_TH) || iters > EM_ITERS {
            break;
        }
    }

    theta0
}

/// Compute the change in the parameter vectors, returning the largest relative and absolute change, respectively.
/// Only parameters with a value greater than rel_thresh (if set), are counted in the relative change check.
fn diffs(t1: &[f64], t2: &[f64], rel_thresh: Option<f64>) -> (f64, f64) {

    let mut max_abs_diff = 0.0;
    let mut max_rel_diff = 0.0;

    for i in 0 .. t1.len() {
        let old_weights = t1[i];
        let new_weights = t2[i];

        let abs_diff = (old_weights - new_weights).abs();
        let rel_diff = abs_diff / old_weights;

        if abs_diff > max_abs_diff {
            max_abs_diff = abs_diff;
        }

        if rel_thresh.map_or(true, |thresh| new_weights > thresh) && rel_diff > max_rel_diff {
            max_rel_diff = rel_diff
        }
    }

    (max_rel_diff, max_abs_diff)
}

#[allow(clippy::cognitive_complexity)]
pub fn em_wrapper(hla_index: PathBuf, hla_counts: PathBuf, cds_db: PathBuf, gen_db: PathBuf, outdir: &str) -> Result<(PathBuf,PathBuf), Error> {
    let index: pseudoaligner::Pseudoaligner<config::KmerType> = utils::read_obj(&hla_index)?;
    let mut eq_counts: EqClassDb = utils::read_obj(&hla_counts)?;
    let weights_name: PathBuf = [outdir, "weights.tsv"].iter().collect();
    let mut weights_file = BufWriter::new(File::create(weights_name)?);
    let pairs_name: PathBuf = [outdir, "pairs.tsv"].iter().collect();
    let mut pairs_file = BufWriter::new(File::create(pairs_name)?);
    let mut alleles_called: HashSet<String> = HashSet::new();
    
    let (eq_class_counts, reads_explained) = eq_counts.eq_class_counts(index.tx_names.len());
    let weights = squarem(&eq_class_counts);
    let weight_names : Vec<(usize, f64, &String, usize)> = 
            weights.
            into_iter().
            enumerate().
            map(|(i,w)| (i, w, &index.tx_names[i], reads_explained[i])).
            collect();
    let mut weight_names : Vec<(u32, f64, &String, &str, usize)> = 
            weight_names.
            into_iter().
            map(|(i,w,n,e)| (i as u32, w, n, n.split('*').next().unwrap(), e)).
            collect();
    weight_names.sort_by(|(_, wa,_, _, _), (_, wb, _, _, _)| (-wa).partial_cmp(&-wb).unwrap()); //sort by descending weight
    weight_names.sort_by_key(|x| x.3); //sort by gene name, stable wrt weight
    
    for (gene, weights) in &weight_names.iter().group_by(|x| (x.3)) {
        info!("Evaluating weights for gene {}",gene);
        let mut weights_written = 0;
        let weights : Vec<_> = weights.collect();
        for (_,w,n,_,exp) in &weights {
            if *exp > MIN_READS_CALL {
                writeln!(weights_file, "{}\t{}\t{}", n, w, f64::from(*exp as u32) / f64::from(eq_class_counts.nreads))?;
                weights_written += 1;
                if weights_written == WEIGHTS_TO_OUTPUT {break;}
            } else {break;}
        }
        if weights_written == 0 {
            continue;
        }
        let mut pairs : Vec<(&String,&String,f64)> = Vec::new();
        'outer: for i in 0..weights_written-1 {
            for j in i+1..weights_written {
                let exp = eq_class_counts.pair_reads_explained(weights[i].0, weights[j].0);
                if exp > MIN_READS_CALL as u32 {
                    pairs.push((weights[i].2, weights[j].2, f64::from(exp) / f64::from(eq_class_counts.nreads)));
                } else {break 'outer;}
            }
        }
        if !pairs.is_empty() {
            pairs.sort_by(|(_,_,wa), (_,_,wb)| (-wa).partial_cmp(&-wb).unwrap());
            for p in pairs.iter().take(std::cmp::min(PAIRS_TO_OUTPUT,pairs.len())) {
                writeln!(pairs_file, "{}\t{}\t{}",p.0, p.1, p.2 )?;
            }
            alleles_called.insert(pairs[0].0.clone());
            alleles_called.insert(pairs[0].1.clone());
        } else if !weight_names.is_empty() {
            alleles_called.insert(weights[0].2.clone());
        }
    }
    
    info!("Writing CDS of {} alleles to file", alleles_called.len());
    use bio::io::fasta;
    let gen_name : PathBuf = [outdir, "gen_pseudoaln.fasta"].iter().collect();
    let mut gen_file = fasta::Writer::to_file(gen_name.clone())?;
    let cds_name : PathBuf = [outdir, "cds_pseudoaln.fasta"].iter().collect();
    let mut cds_file = fasta::Writer::to_file(cds_name.clone())?;
    let cds_db_file = fasta::Reader::from_file(&cds_db)?;
    let gen_db_file = fasta::Reader::from_file(&gen_db)?;
    
    for result in cds_db_file.records() {
        let record = result?;
        let allele_str = record.desc().ok_or_else(|| format_err!("no HLA allele"))?;
        let allele_str = allele_str.split(' ').next().ok_or_else(||format_err!("no HLA allele"))?;
        if alleles_called.contains(allele_str) {
            info!("Writing CDS of {}",allele_str);
            cds_file.write_record(&record)?;
        }
    }
    cds_file.flush()?;
    info!("Writing genomic sequence of {} alleles to file", alleles_called.len());
    for result in gen_db_file.records() {
        let record = result?;
        let allele_str = record.desc().ok_or_else(|| format_err!("no HLA allele"))?;
        let allele_str = allele_str.split(' ').next().ok_or_else(||format_err!("no HLA allele"))?;
        if alleles_called.contains(allele_str) {
            info!("Writing genomic sequence of {}",allele_str);
            gen_file.write_record(&record)?;
        }
    }
    gen_file.flush()?;
    Ok((gen_name, cds_name))
}

#[cfg(test)]
mod test_em {
    use super::*;
    use config::EqClass;

    fn test_ds() -> EqClassCounts {
        let mut counts = HashMap::new();

        let eq_a = EqClass::from(vec![0]);
        let eq_ab = EqClass::from(vec![0,1]);

        let eq_c = EqClass::from(vec![2]);
        let eq_d = EqClass::from(vec![3]);

        counts.insert(eq_a, 1);
        counts.insert(eq_ab, 19);
        counts.insert(eq_c, 10);
        counts.insert(eq_d, 10);
        
        EqClassCounts { counts_umi : counts, counts_reads : HashMap::new(), nitems: 4, nreads: 0 }
    }


    fn test2_ds() -> EqClassCounts {
        let mut counts = HashMap::new();

        let eq_a = EqClass::from(vec![0]);
        let eq_ab = EqClass::from(vec![0,1]);

        let eq_c = EqClass::from(vec![2]);
        let eq_d = EqClass::from(vec![3]);

        let eq_e = EqClass::from(vec![4,5]);

        counts.insert(eq_a, 1);
        counts.insert(eq_ab, 19);
        counts.insert(eq_c, 10);
        counts.insert(eq_d, 10);
        counts.insert(eq_e, 20);
        
        EqClassCounts { counts_umi : counts, counts_reads : HashMap::new(), nitems: 6, nreads: 0 }
    }


    #[test]
    fn simple_inf() {

        let eq_c = test_ds();
        let res = em(&eq_c);

        println!("{:?}", res);
    }

    #[test]
    fn med_inf() {

        let eq_c = test2_ds();
        let res = em(&eq_c);

        println!("{:?}", res);
    }


    #[test]
    fn accel_inf() {

        let eq_c = test_ds();
        let res = squarem(&eq_c);

        println!("{:?}", res);
    }


}
