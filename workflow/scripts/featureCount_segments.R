# This script calculates RPF expression on CDS
library(Rsubread)
library(dplyr)
library(mgsub)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Count RPFs (normalized in RPKM) on CDS for each gene, using `featureCounts`
## run all bams together
samples <- read.table(snakemake@input[["samples"]], header=T)
bamfiles <- paste0("./merged_bam/", as.vector(samples$sample),"_merged_dedup_sorted.bam")

## run one bam file
# bamfiles <- snakemake@input[["bamfile"]]
RPFcounts <- featureCounts(files=bamfiles, annot.ext=snakemake@input[['saf']],
    isGTFAnnotationFile=FALSE, fracOverlap=1,
    strandSpecific=snakemake@params[["strand"]], countMultiMappingReads=FALSE, juncCounts=TRUE, 
    nthreads=snakemake@threads[[1]])

write.table(RPFcounts$counts, file=snakemake@output[[1]], sep="\t", quote=F, row.names = TRUE, col.names = NA)


