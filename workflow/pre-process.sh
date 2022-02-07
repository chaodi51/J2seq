## de-multiplex to include UMI
#Use picard tools for demultiplexing as described here: 
#https://sfvideo.blob.core.windows.net/sitefinity/docs/default-source/user-guide-manual/demultiplexing-data-containing-unique-molecular-indexes-(umis).pdf?sfvrsn=99953207_26


unaligned_bam=/home/dic/weitzman_lab/projects/Price_J2_2020_12/data/raw/unaligned_bam
for i in `sed 1d sample_table.tsv | cut -f1`;do 
	echo $i;
	echo "
	#!/bin/bash

	## cat mapped fastqs
	#cat $unaligned_bam/L001/$i.fastq $unaligned_bam/L002/$i.fastq $unaligned_bam/L003/$i.fastq $unaligned_bam/L004/$i.fastq > raw_fq/$i.fastq
	#gzip raw_fq/$i.fastq

	## cat unmapped bams
	# merge the headers to maintain @RG header
	samtools view -H $unaligned_bam/L001/$i\_unmapped.bam > $i.foo1
	samtools view -H $unaligned_bam/L002/$i\_unmapped.bam | tail -1 > $i.foo2
	samtools view -H $unaligned_bam/L003/$i\_unmapped.bam | tail -1 > $i.foo3
	samtools view -H $unaligned_bam/L004/$i\_unmapped.bam | tail -1 > $i.foo4
	cat $i.foo1 $i.foo2 $i.foo3 $i.foo4 > unmapped_bam/$i.header
	rm -f $i.foo*

	samtools cat -h unmapped_bam/$i.header $unaligned_bam/L001/$i\_unmapped.bam $unaligned_bam/L002/$i\_unmapped.bam $unaligned_bam/L003/$i\_unmapped.bam $unaligned_bam/L004/$i\_unmapped.bam -o unmapped_bam/$i\_unaligned.bam
	" > qsub/$i.sh
	qsub -e qsub/$i.e -o qsub/$i.o -l h_vmem=10G -l mem_free=10G -pe smp 6 -V -cwd qsub/$i.sh
done


