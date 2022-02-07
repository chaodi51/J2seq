
rule featureCount:
    input:
        samples = "sample_contrast.tsv",
        gtf = gtf_file
    output:
        "../results/J2seq_featureCount.tsv",
    log:
        "logs/featureCount.log"   
    threads: 8     
    params: 
      mem = '10G',
      jobName = "featureCount",
      featureType = "gene",
      strand = 1
    script:
        "../scripts/featureCount.R"

## count reads on the sense strand of annotated mRNAs (using featureType=exon)
rule featureCount_ad5_fwd:
    input:
        samples = "sample_contrast.tsv",
        gtf = "/home/dic/public/genomes/GENCODE/Ad5/Ad5_v9.1.Simple.gtf"
    output:
        "../results/J2seq_featureCount_ad5_fwd.tsv",
    log:
        "logs/featureCount_ad5_fwd.log"   
    threads: 8     
    params: 
      mem = '10G',
      jobName = "featureCount_ad5_fwd",
      featureType = "exon",
      strand = 1
    script:
        "../scripts/featureCount.R"

## count reads on the antisense strand of annotated mRNAs
rule featureCount_ad5_rev:
    input:
        samples = "sample_contrast.tsv",
        gtf = "/home/dic/public/genomes/GENCODE/Ad5/Ad5_v9.1.Simple.gtf"
    output:
        "../results/J2seq_featureCount_ad5_rev.tsv",
    log:
        "logs/featureCount_ad5_rev.log"   
    threads: 8     
    params: 
      mem = '10G',
      jobName = "featureCount_ad5_rev",
      featureType = "exon",
      strand = 2
    script:
        "../scripts/featureCount.R"

rule featureCount_ad5_segments:
    input:
        samples = "sample_contrast.tsv",
        saf = "polyA_read_through_segments.saf"
    output:
        "../results/J2seq_featureCount_ad5_segments.tsv"
    log:
        "logs/featureCount_ad5_segments.log"   
    threads: 8     
    params: 
      mem = '10G',
      jobName = "featureCount_ad5_segments",
      strand = 1
    script:
        "../scripts/featureCount_segments.R"