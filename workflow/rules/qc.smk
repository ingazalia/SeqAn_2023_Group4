
rule fastqc_raw:
    input: lambda wildcards: samples.at[wildcards.sample, 'fq' + wildcards.number]
    output:
        html = "results/quality/raw/{sample}_{number}_fastqc.html",
        zip = "results/quality/raw/{sample}_{number}_fastqc.zip"
    log:
        "results/logs/fastqc/{sample}_{number}_raw.log"
    conda:
        "../envs/env.yaml"
    threads: 6
    wrapper:
        "v1.31.1/bio/fastqc"

rule multiqc_raw:
    input:
        expand("results/quality/raw/{sample}_{number}_fastqc.zip", sample=sample, number=['1', '2'])
    output:
        "results/quality/raw/multiqc_report.html"
    log: "results/logs/multi/raw.log"
    conda:
        "../envs/env.yaml"
    shell:
        "multiqc -o results/quality/raw results/quality/raw &> {log}"

rule trimgalore:
    input:
        fastq1 = lambda wildcards: samples.at[wildcards.sample, 'fq1'],
        fastq2 = lambda wildcards: samples.at[wildcards.sample, 'fq2']
    output:
        "results/trimmed/{sample}_R1_val_1.fq.gz",
        "results/trimmed/{sample}_R2_val_2.fq.gz",
        "results/trimmed/{sample}_R1.fastq.gz_trimming_report.txt",
        "results/trimmed/{sample}_R2.fastq.gz_trimming_report.txt"
    params:
        extra=config["trimgalore"]
    log:
        "results/logs/trim_galore/{sample}.log"
    threads: 4
    conda:
        "../envs/env.yaml"
    wrapper:
        "v1.31.1/bio/trim_galore/pe"

rule fastqc_trimming:
    input: "results/trimmed/{sample}_R{number}_val_{number}.fq.gz"
    output:
        html = "results/quality/trimming/{sample}_{number}_fastqc.html",
        zip = "results/quality/trimming/{sample}_{number}_fastqc.zip"
    log:
        "results/logs/fastqc/trimming/{sample}_{number}.log"
    conda:
        "../envs/env.yaml"
    threads: 6
    wrapper:
        "v1.31.1/bio/fastqc"

rule multiqc_trimming:
    input:
        expand("results/quality/trimming/{sample}_{number}_fastqc.zip", sample=sample, number=['1', '2'])
    output:
        "results/quality/trimming/multiqc_report.html"
    log: "results/logs/multi/trim.log"
    conda:
        "../envs/env.yaml"
    shell:
        "multiqc -o results/quality/trimming results/quality/trimming &> {log}"

rule sortmerna:
    input:
        trim1="results/trimmed/{sample}_R1_val_1.fq.gz",
        trim2="results/trimmed/{sample}_R2_val_2.fq.gz",
        ref_rRNA_1=config["SortMeRNA"]["rRNA_REF_1"],
        ref_rRNA_2=config["SortMeRNA"]["rRNA_REF_2"],
        ref_rRNA_3=config["SortMeRNA"]["rRNA_REF_3"]
    output:
        nonrRNA_1="results/sortmerna/{sample}_fwd.fq.gz",
        nonrRNA_2="results/sortmerna/{sample}_rev.fq.gz",
        rRNA_1="results/sortmerna/aligned/{sample}_fwd.fq.gz",
        rRNA_2="results/sortmerna/aligned/{sample}_rev.fq.gz"
    log:
        "results/logs/sortmerna/{sample}.log"
    params:
        numAlign = config["sortmeRNA"]["numAlign"],
        other = config["sortmeRNA"]["other"]
    threads: 4
    conda:
        "../envs/env.yaml"
    shell:
        """
            sortmerna --ref {input.ref_rRNA_1} --ref {input.ref_rRNA_2} --ref {input.ref_rRNA_3} --reads {input.trim1} \
            --reads {input.trim2} --other results/sortmerna/{wildcards.sample} --aligned results/sortmerna/aligned/{wildcards.sample} \
            --fastx --fast -a {threads} -num_alignments {params.numAlign} {params.other} --workdir results/sortmerna --out2 True \
            --kvdb results/sortmerna/{wildcards.sample} --idx-dir results/sortmerna/{wildcards.sample} \
            --readb results/sortmerna/{wildcards.sample} --paired_in &> {log}
        """

