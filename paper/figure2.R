
library(Matrix)
library(Seurat)
library(ggplot2)
library(gridExtra)
options(bitmapType='cairo')

# Code to generate Figure 2 (Petti AML Subject 809653 DRB1 case study)

# Get the cell type labels
# R data object from https://zenodo.org/record/3066262#.XUnM5ZNKijh
AML3 <- readRDS("/mnt/park1/compbio/HLA/aml_seurat/809653.seurat.rds")
celltypes.paper <- AML3@meta.data$CellType
names(celltypes.paper) <- dimnames(AML3@data)[[2]]
rm(AML3)

# Data processed with Cell Ranger 2.1.1
CELLRANGERDIR <- "/mnt/analysis/marsoc/pipestances/H3TNHDSXX/SC_RNA_COUNTER_PD/64474/HEAD/outs/filtered_gene_bc_matrices_mex/GRCh38/"
# Output of scHLAcount using true genotype FASTA files
ANALYSISDIR <- "/mnt/park1/compbio/HLA/only-count-true-ref/AML3"

# Read Cell Ranger data
cells.data <- Read10X(data.dir = CELLRANGERDIR)
cells <- CreateSeuratObject(counts = cells.data)

# Add cell type labels as metadata column 
cells <- AddMetaData(object=cells,col.name="celltype",metadata=celltypes.paper)

# Perform pre-processing with Seurat
cells <- NormalizeData(cells)
cells <- FindVariableFeatures(cells, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(cells)
cells <- ScaleData(cells, features = all.genes)

cells <- RunPCA(cells, features = VariableFeatures(object = cells))
cells <- FindNeighbors(cells, dims = 1:10)
cells <- FindClusters(cells, resolution = 0.5)
cells <- RunTSNE(cells,check_duplicates = F)

# Read scHLAcount data
mat <- readMM(file = paste0(ANALYSISDIR, "/count_matrix.mtx"))
feature.names  <- read.delim(paste0(ANALYSISDIR,"/labels.tsv"), header = FALSE, stringsAsFactors = FALSE)
dimnames(mat)[[1]] <- feature.names$V1
dimnames(mat)[[2]] <- dimnames(cells)[[2]]

# Normalize counts and insert into matrix as metadata columns
c <- median(cells$nCount_RNA)
curr_gene <- ""
gene_names <- c()
has_two_alleles <- c()
n_seen = 1

for (i in seq(1,length(feature.names$V1))) {
  f <- feature.names$V1[i]
  s <- strsplit(f, "*", fixed=T)
  if (s[[1]][1] == curr_gene) {
    n_seen <- n_seen + 1
  } else {
    if (n_seen == 1 && curr_gene != "") {
      #one allele
      gene_names <- c(gene_names, curr_gene)
      cells <- AddMetaData(object=cells,col.name=paste0("gene",curr_gene,"sum"),metadata=log2(1+c*mat[i-1,]/cells$nCount_RNA))
    } else if (n_seen == 3) {
      #two alleles
      gene_names <- c(gene_names, curr_gene)
      has_two_alleles <- c(has_two_alleles, curr_gene)
      cells <- AddMetaData(object=cells,col.name=paste0("gene",curr_gene,"sum"),metadata=log2(1+c*(mat[i-1,]+mat[i-2,]+mat[i-3,])/cells$nCount_RNA))
      cells <- AddMetaData(object=cells,col.name=paste0("gene",curr_gene),metadata=log2(1+c*mat[i-1,]/cells$nCount_RNA))
      cells <- AddMetaData(object=cells,col.name=paste0("gene",curr_gene,"allele1"),metadata=log2(1+c*mat[i-3,]/cells$nCount_RNA))
      cells <- AddMetaData(object=cells,col.name=paste0("gene",curr_gene,"allele2"),metadata=log2(1+c*mat[i-2,]/cells$nCount_RNA))
      x <- paste0("gene",curr_gene,"allele1")
      y <- paste0("gene",curr_gene,"allele2")
      cells <- AddMetaData(object=cells,col.name=paste0("gene",curr_gene,"ratio"),metadata=cells@meta.data[x]/(cells@meta.data[x]+cells@meta.data[y]))
      n_seen = 1
    }
    curr_gene <- s[[1]][1]
  }
}
if (n_seen == 1 && curr_gene != "") {
  #one allele
  gene_names <- c(gene_names, curr_gene)
  cells <- AddMetaData(object=cells,col.name=paste0("gene",curr_gene,"sum"),metadata=log2(1+c*mat[i-1,]/cells$nCount_RNA))
} else if (n_seen == 3) {
  #two alleles
  gene_names <- c(gene_names, curr_gene)
  has_two_alleles <- c(has_two_alleles, curr_gene)
  cells <- AddMetaData(object=cells,col.name=paste0("gene",curr_gene,"sum"),metadata=log2(1+c*(mat[i-1,]+mat[i-2,]+mat[i-3,])/cells$nCount_RNA))
  cells <- AddMetaData(object=cells,col.name=paste0("gene",curr_gene),metadata=log2(1+c*mat[i-1,]/cells$nCount_RNA))
  cells <- AddMetaData(object=cells,col.name=paste0("gene",curr_gene,"allele1"),metadata=log2(1+c*mat[i-3,]/cells$nCount_RNA))
  cells <- AddMetaData(object=cells,col.name=paste0("gene",curr_gene,"allele2"),metadata=log2(1+c*mat[i-2,]/cells$nCount_RNA))
  x <- paste0("gene",curr_gene,"allele1")
  y <- paste0("gene",curr_gene,"allele2")
  cells <- AddMetaData(object=cells,col.name=paste0("gene",curr_gene,"ratio"),metadata=cells@meta.data[x]/(cells@meta.data[x]+cells@meta.data[y]))
}

CELLTYPE_PLOT <- DimPlot(object=cells,group.by='celltype', label=F, cells=rownames(cells@meta.data)[!is.na(cells$celltype)]) + ggtitle("(c) Cell types from Paulson et. al")

# Bug in Seurat related to plotting, fails when a NA is first: https://github.com/satijalab/seurat/issues/1853
# Temporary fix: reordering the rows so a NA doesn't come first
new.cell.order <- rownames(cells@meta.data)[order(cells@meta.data["geneDRB1ratio"])]
DRB1_ASE_PLOT <- FeaturePlot(cells, "geneDRB1ratio", reduction="tsne",cols=c('blue','red'), cells = new.cell.order) + ggtitle("(b) Fraction of molecules assigned to allele 01:03")

DRB1_EXPR_PLOT <- FeaturePlot(cells, "geneDRB1sum", reduction="tsne", cells = new.cell.order) + ggtitle("(a) HLA-DRB1 normalized expression")

COMBINED <- grid.arrange(DRB1_EXPR_PLOT,DRB1_ASE_PLOT,CELLTYPE_PLOT,
             layout_matrix=rbind(c(1,1,1,2,2,2),c(3,3,3,3,NA,NA)))

ggsave("figure2.pdf", COMBINED, width=10, height=8, units="in")

