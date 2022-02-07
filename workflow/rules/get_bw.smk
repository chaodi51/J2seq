## convert the bam file to bigwig format (signal normalized to RPM) for viewing the data on Genome Browse
rule get_bw:
    input: 
        "merged_bam/{sample}_merged_dedup_sorted.bam"
    output: 
        fwd = "bw_cpm/{sample}.fwd.bw",
        rev = "bw_cpm/{sample}.rev.bw"
    log: "logs/get_bw/{sample}.log"
    threads: 4
    params:
        mem = '5G',
        jobName = "get_bw.{sample}"
    shell:
        ## strand of the data is opposite to dUTP method, as the UMI adapter process flipped the strandedness 
        '''
        bamCoverage -b {input} -o {output.fwd} --exactScaling --normalizeUsing CPM -bs 1 \\
        --numberOfProcessors {threads} --filterRNAstrand reverse \\
        --blackListFileName ~/public/genomes/GENCODE/human/GRCh38_unified_blacklist.bed &> {log}
        bamCoverage -b {input} -o {output.rev} --exactScaling --normalizeUsing CPM -bs 1 \\
        --numberOfProcessors {threads} --filterRNAstrand forward \\
        --blackListFileName ~/public/genomes/GENCODE/human/GRCh38_unified_blacklist.bed &> {log}

        '''

rule rev_bw:
    input: 
        "bw_cpm/{sample}.rev.bw"
    output: 
        "bw_cpm/{sample}.rev.minus.bw"
    log: "logs/rev_bw/{sample}.log"
    threads: 1
    params:
        chrom_size="/home/dic/public/genomes/GENCODE/human_ad5/STAR_index_gencode_v35_ad5/chrNameLength.txt",
        mem = '5G',
        jobName = "rev_bw.{sample}"
    shell:
        ## change signals in rev bw to negative values to show on IGV 
        '''
        bigWigToBedGraph {input} {wildcards.sample}.bg
        bedSort {wildcards.sample}.bg {wildcards.sample}.sorted.bg
        awk 'BEGIN{{FS=OFS="\t"}} {{if($4>0) $4=-$4; print $0}}' {wildcards.sample}.sorted.bg > {wildcards.sample}.sorted.minus.bg
        bedGraphToBigWig {wildcards.sample}.sorted.minus.bg {params.chrom_size} {output} &> {log}
        rm -f {wildcards.sample}.bg {wildcards.sample}.sorted.bg {wildcards.sample}.sorted.minus.bg
        '''
