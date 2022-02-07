# run STAR mapping with only uniquely mapped reads
rule star_map:
    input:
        "raw_fq/{sample}.fastq.gz",
    output:
        "STAR_align/{sample}.bam",
        "STAR_align/{sample}.Log.final.out",
    params:
        genome_dir = STAR_index,
        annotation = gtf_file,
        mem = '40G',
        jobName = "star_map.{sample}" 
    threads: 4
    shell:
    # also output the unmapped reads, use --outFilterMultimapNmax 10(default) next time
        "STAR --runThreadN {threads} --sjdbGTFfile {params.annotation} --genomeDir {params.genome_dir} " 
        "--outSAMmapqUnique 255 --outSJfilterReads Unique "
        "--outFileNamePrefix STAR_align/{wildcards.sample}. --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes All "
        "--quantMode GeneCounts --twopassMode Basic "
        "--readFilesIn {input} --readFilesCommand gunzip -c && "
        "mv STAR_align/{wildcards.sample}.Aligned.sortedByCoord.out.bam STAR_align/{wildcards.sample}.bam"

rule samtools_index:
    input:
        "STAR_align/{sample}.bam"
    output:
        "STAR_align/{sample}.bam.bai"
    params: 
      mem = '10G',
      jobName = "samtools_index.{sample}"
    shell:
        "samtools index {input} {output}"

rule multiqc_star:
    input: expand('STAR_align/{sample}.Log.final.out', sample=SAMPLES)
    output: '../results/multiqc/star/multiqc_report.html'
    params:
        mem = '6G',
        jobName = 'multiqc_fastqc',
        outdir = "../results/multiqc/star"    
    log: 'logs/multiqc/star/multiqc.log'
    shell:
        '''
        rm -r ../results/multiqc/star/
        multiqc STAR_align -o {params.outdir} &> {log}
        '''