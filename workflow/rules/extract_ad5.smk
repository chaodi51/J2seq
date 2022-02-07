# deduplication by UMI
rule extract_ad5:
    input:
        "merged_bam/{sample}_merged_dedup_sorted.bam"
    output:
        "merged_bam/ad5.{sample}_merged_dedup_sorted.bam"
    params:
        mem = '10G', 
        jobName = "extract_ad5.{sample}"
    log: "logs/extract_ad5/{sample}.log"     
    threads: 4
    shell:
        '''
            samtools view -h {input} Ad5 -b > {output} 2> {log}
        '''

rule samtools_index_ad5:
    input:
        "merged_bam/ad5.{sample}_merged_dedup_sorted.bam"
    output:
        bai = "merged_bam/ad5.{sample}_merged_dedup_sorted.bam.bai"
    params: 
      mem = '10G',
      jobName = "samtools_index_ad5.{sample}"
    log: "logs/samtools_index_ad5/{sample}.log"       
    shell:
        '''
        samtools index {input} {output.bai}
        '''
