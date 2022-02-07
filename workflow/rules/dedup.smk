# deduplication by UMI
rule dedup:
    input:
        "merged_bam/{sample}_merged.bam"
    output:
        bam = "merged_bam/{sample}_merged_dedup.bam"
    params:
        mem = '50G', # Ad5input2 needs 100G
        jobName = "dedup_bam.{sample}"
    log: "logs/dedup_bam/{sample}.log"     
    threads: 4
    shell:
        '''
        umi_tools dedup -I {input} --extract-umi-method=tag --umi-tag=RX -S {output.bam} \\
        --output-stats=merged_bam/{wildcards.sample}_merged_dedup.stats --mapping-quality=255 --no-sort-output \\
        --log={log}
        '''

rule samtools_sort_index:
    input:
        "merged_bam/{sample}_merged_dedup.bam"
    output:
        sorted_bam = "merged_bam/{sample}_merged_dedup_sorted.bam",
        bai = "merged_bam/{sample}_merged_dedup_sorted.bam.bai"
    params: 
      mem = '10G',
      jobName = "samtools_sort_index.{sample}"
    log: "logs/dedusamtools_sort_indexp_bam/{sample}.log"       
    shell:
        '''
        samtools sort {input} -o {output.sorted_bam}
        samtools index {output.sorted_bam} {output.bai}
        '''
