
# See figure3.R for example processing workflow 

# PBMC (Normal cells)
load("/mnt/park1/compbio/HLA/paulson_robj/paulson_disc_pbmc.Rdata")

# Average normalized expression
mean(na.omit(PBMC$geneAsumnolog))
mean(na.omit(PBMC$geneBsumnolog))
mean(na.omit(PBMC$geneCsumnolog))

# Molecules assigned to allele 1
sum(na.omit(PBMC$geneAallele1))/(sum(na.omit(PBMC$geneAallele1)) + sum(na.omit(PBMC$geneAallele2)))
sum(na.omit(PBMC$geneBallele1))/(sum(na.omit(PBMC$geneBallele1)) + sum(na.omit(PBMC$geneBallele2)))
sum(na.omit(PBMC$geneCallele1))/(sum(na.omit(PBMC$geneCallele1)) + sum(na.omit(PBMC$geneCallele2)))


# Tumor (mixture of tumor and normal cells)
load("/mnt/park1/compbio/HLA/paulson_robj/paulson_disc_tumor.Rdata")

# Tumor cells

# Average normalized expression
mean(na.omit(tumor$geneAsumnolog[tumor$cellTypes1 == "Tumor"]))
mean(na.omit(tumor$geneBsumnolog[tumor$cellTypes1 == "Tumor"]))
mean(na.omit(tumor$geneCsumnolog[tumor$cellTypes1 == "Tumor"]))

# Molecules assigned to allele 1
sum(na.omit(tumor$geneAallele1[tumor$cellTypes1 == "Tumor"]))/(sum(na.omit(tumor$geneAallele1[tumor$cellTypes1 == "Tumor"])) + sum(na.omit(tumor$geneAallele2[tumor$cellTypes1 == "Tumor"])))
sum(na.omit(tumor$geneBallele1[tumor$cellTypes1 == "Tumor"]))/(sum(na.omit(tumor$geneBallele1[tumor$cellTypes1 == "Tumor"])) + sum(na.omit(tumor$geneBallele2[tumor$cellTypes1 == "Tumor"])))
sum(na.omit(tumor$geneCallele1[tumor$cellTypes1 == "Tumor"]))/(sum(na.omit(tumor$geneCallele1[tumor$cellTypes1 == "Tumor"])) + sum(na.omit(tumor$geneCallele2[tumor$cellTypes1 == "Tumor"])))

# Normal cells

# Average normalized expression
mean(na.omit(tumor$geneAsumnolog[tumor$cellTypes1 != "Tumor"]))
mean(na.omit(tumor$geneBsumnolog[tumor$cellTypes1 != "Tumor"]))
mean(na.omit(tumor$geneCsumnolog[tumor$cellTypes1 != "Tumor"]))

# Molecules assigned to allele 1
sum(na.omit(tumor$geneAallele1[tumor$cellTypes1 != "Tumor"]))/(sum(na.omit(tumor$geneAallele1[tumor$cellTypes1 != "Tumor"])) + sum(na.omit(tumor$geneAallele2[tumor$cellTypes1 != "Tumor"])))
sum(na.omit(tumor$geneBallele1[tumor$cellTypes1 != "Tumor"]))/(sum(na.omit(tumor$geneBallele1[tumor$cellTypes1 != "Tumor"])) + sum(na.omit(tumor$geneBallele2[tumor$cellTypes1 != "Tumor"])))
sum(na.omit(tumor$geneCallele1[tumor$cellTypes1 != "Tumor"]))/(sum(na.omit(tumor$geneCallele1[tumor$cellTypes1 != "Tumor"])) + sum(na.omit(tumor$geneCallele2[tumor$cellTypes1 != "Tumor"])))
