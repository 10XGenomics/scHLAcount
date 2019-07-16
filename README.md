# scHLAcount
Count HLA alleles in single-cell RNA-seq data  

```
scHLAcount DEV
Charlotte Darby <cdarby@jhu.edu> and Ian Fiddes <ian.fiddes@10xgenomics.com> and Patrick Marks <patrick@10xgenomics.com>
HLA genotyping and allele-specific expression for single-cell RNA sequencing

USAGE:
    sc_hla_count [FLAGS] [OPTIONS] --bam <FILE> --cell-barcodes <FILE>

FLAGS:
    -h, --help                  Prints help information
        --primary-alignments    If specified, will use primary alignments only
        --unmapped              If specified, will also use unmapped reads for genotyping (very slow!)
    -V, --version               Prints version information

OPTIONS:
    -b, --bam <FILE>                    Cellranger BAM file
    -c, --cell-barcodes <FILE>          File with cell barcodes to be evaluated
    -f, --fasta-cds <FILE>              Multi-FASTA file with CDS sequence of each allele [default: ]
    -g, --fasta-genomic <FILE>          Multi-FASTA file with genomic sequence of each allele [default: ]
    -d, --hladb-dir <PATH>              Directory of the IMGT-HLA database [default: ]
    -i, --hla-index <FILE>              debruijn_mapping pseudoalignment index file constructed from IMGT-HLA database [default: ]
        --log-level <log_level>         Logging level [default: error]  [possible values: info, debug, error]
    -o, --out-dir <OUTPUT_DIR>           [default: hla-typer-results]
        --pl-tmp <PSEUDOALIGNER_TMP>    Directory to write the pseudoaligner temporary files generated [default: .]
    -r, --region <STRING>               Samtools-format region string of reads to use [default: 6:28510120-33480577]
```  

