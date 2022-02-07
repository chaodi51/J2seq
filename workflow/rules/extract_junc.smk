# deduplication by UMI
rule extract_junc:
    input:
        bam = "merged_bam/{sample}_merged_dedup_sorted.bam",
    output:
        all_junc = "merged_bam/{sample}_merged_dedup_sorted.junc.bed",
        ad5_junc = "merged_bam/ad5.{sample}_merged_dedup_sorted.junc.bed"
    params:
        mem = '10G', 
        jobName = "extract_junc.{sample}"
    log: "logs/extract_junc/{sample}.log"     
    threads: 1
    shell:
        '''
            regtools junctions extract -s 1 {input.bam} -o {output.all_junc} &> {log}
            # junc from Ad5
            cat {output.all_junc} |grep Ad5 > {output.ad5_junc} 2>> {log}
        '''

rule anno_junc:
    input:
        all_junc = "merged_bam/{sample}_merged_dedup_sorted.junc.bed",
        gtf = gtf_file,
        fa = genome_seq
    output:
        all_anno = "merged_bam/{sample}_merged_dedup_sorted.junc.anno",
        ad5_anno = "merged_bam/ad5.{sample}_merged_dedup_sorted.junc.anno"
    params:
        mem = '10G', 
        jobName = "anno_junc.{sample}"
    log: "logs/anno_junc/{sample}.log"     
    threads: 1
    shell:
        '''
            regtools junctions annotate -o {output.all_anno} {input.all_junc} {input.fa} {input.gtf} &> {log}
            # junc from Ad5
            cat {output.all_anno} |grep Ad5 > {output.ad5_anno} 2>> {log}
        '''

rule junc_stat:
    input:
        expand("merged_bam/ad5.{sample}_merged_dedup_sorted.junc.anno", sample=SAMPLES),
    output:
        "../results/Junctions/J2seq.ad5.junc.tsv"        
    log: "logs/junc_stat.log"
    params:
        mem = '5G',
        jobName = "junc_stat"
    shell:
        '''
        rm -f {output};
        
        echo "sample\ttotal_junc" > {output};
        for i in {input}; do
            name=`echo $i | sed 's/merged_bam\/ad5.//g;s/_merged_dedup_sorted.junc.anno//g'`;
            awk 'NR>1 && $5>1 {{sum+=$5}}; END {{print "'$name'\t"sum}}' $i > foo;

            cat {output} foo > foo1;
            mv foo1 {output};
        done
        rm -f foo foo1
        '''