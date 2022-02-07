rule fastq_screen:
    input: 
        get_fastq ## call the function in main Snakefile
    output:
        "../results/fastq_screen_output/{sample}.tagged_filter.fastq.gz",
    log: "logs/fastq_screen/{sample}.log"
    threads: 4
    params:
        fastq_screen_conf = fastq_screen_conf,
        outdir = fastq_screen_output,
        mem = '10G',
        jobName = "fastq_screen.{sample}" 
    shell:
        # keep mapped reads only from human and Ad5
        '''
            fastq_screen --aligner bowtie2  --subset 1000000 --threads {threads} --conf {params.fastq_screen_conf} \\
            --tag --filter 3300000000 --outdir {params.outdir} {input} &> {log}
        '''

rule multiqc_fastqscreen:
    input: expand('../results/fastq_screen_output/{sample}.tagged_filter.fastq.gz', sample=SAMPLES)
    output: '../results/multiqc/fastqscreen/multiqc_report.html'
    log: 'logs/multiqc/fastqscreen/multiqc.log'
    params:
        mem = '6G',
        jobName = "multiqc_fastqscreen", 
        outdir = "../results/multiqc/fastqscreen"    
    shell:
        '''
        multiqc ../results/fastq_screen_output -o {params.outdir} &> {log}
        '''
    # wrapper: '0.36.0/bio/multiqc'   