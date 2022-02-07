# run STAR mapping with only uniquely mapped reads
rule star_map_multi:
    input:
        "raw_fq/{sample}.fastq.gz",
    output:
        "STAR_align_multi/{sample}.bam",
        "STAR_align_multi/{sample}.Log.final.out",
    params:
        genome_dir = STAR_index,
        annotation = gtf_file,
        mem = '40G',
        jobName = "star_map.{sample}" 
    threads: 4
    shell:
    # also output the unmapped reads, include more multi reads for TEtranscripts to handle multi-mappers
        "STAR --runThreadN {threads} --sjdbGTFfile {params.annotation} --genomeDir {params.genome_dir} " 
        "--outSAMmapqUnique 255 --outSJfilterReads Unique "
        "--outFileNamePrefix STAR_align_multi/{wildcards.sample}. --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes All "
        "--outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 "
        "--quantMode GeneCounts --twopassMode Basic "
        "--readFilesIn {input} --readFilesCommand gunzip -c && "
        "mv STAR_align_multi/{wildcards.sample}.Aligned.sortedByCoord.out.bam STAR_align_multi/{wildcards.sample}.bam"

rule samtools_index:
    input:
        "STAR_align_multi/{sample}.bam"
    output:
        "STAR_align_multi/{sample}.bam.bai"
    params: 
      mem = '10G',
      jobName = "samtools_index.{sample}"
    shell:
        "samtools index {input} {output}"
