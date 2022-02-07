shell.prefix("source ~/.bash_profile; ")

import os

# snakemake -s rules/STAR_index.smk -c "qsub -l h_vmem={params.mem} -l mem_free={params.mem} -pe smp {threads} -V -cwd -e qsub/{params.jobName}.e -o qsub/{params.jobName}.o" -j -p

## STAR indexing the genome
rule star_index:
    input:
        fa = "/home/dic/public/genomes/GENCODE/human_ad5/gencode_v35_ad5.fa", # reference FASTA file
        gtf = "/home/dic/public/genomes/GENCODE/human_ad5/gencode_v35_ad5.gtf" # GTF file
    output:
        directory("/home/dic/public/genomes/GENCODE/human_ad5/STAR_index_gencode_v35_ad5") # the index folder
    threads: 6 # set the maximum number of available cores
    params: 
        mem = "10G",
        jobName = "star_index"
    shell:
        """
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output} \\
        --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --sjdbOverhang 100
        """

