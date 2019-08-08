
# Figure S2
#==> hlaa.bed <==
#6	29942532	29945870

#==> hlab.bed <==
#6	31353875	31357179

#==> hlac.bed <==
#6	31268749	31272092

# samtools view -b /mnt/yard2/ian/sra-data/paulson/SRR7692286/outs/possorted_genome_bam.bam 6 | bedtools coverage -split -a hlaa.bed -b stdin -d | cut -f5 > acov5p.txt
# samtools view -b /mnt/yard2/ian/sra-data/paulson/SRR7692286/outs/possorted_genome_bam.bam 6:31353875-31357179 | bedtools coverage -split -a hlab.bed -b stdin -d | cut -f5 > bcov5p.txt
# samtools view -b /mnt/yard2/ian/sra-data/paulson/SRR7692286/outs/possorted_genome_bam.bam 6:31268749-31272092 | bedtools coverage -split -a hlac.bed -b stdin -d | cut -f5 > ccov5p.txt

# samtools view -b yard/paulson_merged.bam 6:29942532-29945870 | bedtools coverage -split -a hlaa.bed -b stdin -d | cut -f5 > acov3p.txt
# samtools view -b yard/paulson_merged.bam 6:31353875-31357179 | bedtools coverage -split -a hlab.bed -b stdin -d | cut -f5 > bcov3p.txt
# samtools view -b yard/paulson_merged.bam 6:31268749-31272092 | bedtools coverage -split -a hlac.bed -b stdin -d | cut -f5 > ccov3p.txt


acov3p <- read.delim("acov3p.txt",header=F)
bcov3p <- read.delim("bcov3p.txt",header=F)
ccov3p <- read.delim("ccov3p.txt",header=F)

acov5p <- read.delim("acov5p.txt",header=F)
bcov5p <- read.delim("bcov5p.txt",header=F)
ccov5p <- read.delim("ccov5p.txt",header=F)

par(mfrow=c(3,1),mar=c(2,2,2,1),cex.main=2,cex.axis=1.3)
plot(seq(1,3338),sapply(acov5p$V1, function(x) (x-min(acov5p))/(max(acov5p)-min(acov5p))), ylim=c(-0.1,1),type='l',lwd=2,main="HLA-A Read Coverage chr6:29942532-29945870 (+ strand)",ylab="")
lines(seq(1,3338),sapply(acov3p$V1, function(x) (x-min(acov3p))/(max(acov3p)-min(acov3p))), type='l',col='red',lwd=2)
legend("top",c("3' GEX", "5' GEX"),col=c("red","black"),lty=c(1,1),lwd=c(3,3),cex=2)
plot(seq(1,3304),sapply(bcov5p$V1, function(x) (x-min(bcov5p))/(max(bcov5p)-min(bcov5p))), ylim=c(-0.1,1),type='l',lwd=2,main="HLA-B Read Coverage chr6:31353875-31357179 (- strand)",ylab="",xlim=c(3304,0))
lines(seq(1,3304),sapply(bcov3p$V1, function(x) (x-min(bcov3p))/(max(bcov3p)-min(bcov3p))), type='l',col='red',lwd=2)
plot(seq(1,3343),sapply(ccov5p$V1, function(x) (x-min(ccov5p))/(max(ccov5p)-min(ccov5p))), ylim=c(-0.1,1),type='l',lwd=2,main="HLA-C Read Coverage chr6:31268749-31272092 (- strand)",ylab="",xlim=c(3343,0))
lines(seq(1,3343),sapply(ccov3p$V1, function(x) (x-min(ccov3p))/(max(ccov3p)-min(ccov3p))), type='l',col='red',lwd=2)
