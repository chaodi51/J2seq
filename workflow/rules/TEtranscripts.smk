# run STAR mapping with only uniquely mapped reads
# TEcount is better suited than TEtranscripts for usage in the cluster environment, 
# as each sample (e.g. replicates of an experiment) can be quantified on separate nodes. 
# strand should be opposite of illumina TruSeq dUTP-"first-strand cDNA library" ("reverse" for --stranded)
rule TEcount:
    input:
        "merged_bam/{sample}_merged_dedup_sorted.bam"
    output:
        "../results/TEcount/TEcount.{sample}.cntTable"
    log: "logs/TEcount/{sample}.log"
    params:
        gene = gtf_file,
        TE = TE_gtf,
        mem = '20G',
        jobName = "TEcount.{sample}" 
    threads: 8
    shell:
        '''
        TEcount --sortByPos --format BAM --mode multi -b {input} \\
        --GTF {params.gene} --TE {params.TE} \\
        --stranded forward --project "../results/TEcount/TEcount.{wildcards.sample}"
        '''
rule combine_tab:
    input:
        expand("../results/TEcount/TEcount.{sample}.cntTable", sample=SAMPLES)
    output:
        "../results/TEcount/J2seq.TEcount.tsv"        
    log: "logs/combie_tab.log"
    params:
        mem = '5G',
        jobName = "combie_tab"
    shell:
        '''
        rm -f {output};
        cut -f1  ../results/TEcount/TEcount.Ad5input1.cntTable | sed 's/"//g;s/gene\/TE/gene_TE/g' > {output};
        for i in {input}; do
            cut -f2 $i | sed 's/merged_bam\///g;s/_merged_dedup_sorted.bam//g' > foo;
            paste {output} foo > foo1;
            mv foo1 {output};
        done
        rm -f foo foo1
        '''