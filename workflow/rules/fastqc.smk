rule fastqc:
    input: 
        get_fastq ## call the function in main Snakefile
    output:
        qc = "../results/fastqc/{sample}_fastqc.html",
    log: "logs/fastqc/{sample}.log"
    threads: 6
    params:
        outdir = "../results/fastqc",
        mem = '6G',
        jobName = "fastqc.{sample}" 
    shell:
        "fastqc --threads {threads} -o {params.outdir} {input} &> {log} "

rule multiqc_fastqc:
    input: expand('../results/fastqc/{sample}_fastqc.html', sample=SAMPLES)
    output: '../results/multiqc/fastqc/multiqc_report.html'
    params:
        mem = '6G',
        jobName = "multiqc_fastqc",
        outdir = "../results/multiqc/fastqc"
    log: 'logs/multiqc/fastqc/multiqc.log'

    shell:
        '''
        rm -r {params.outdir}
        multiqc ../results/fastqc -o {params.outdir} &> {log}
        '''