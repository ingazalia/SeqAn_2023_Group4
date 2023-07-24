rule hisat_build:
    input:
        fasta = config["ref_FASTA"]
    output:
        genome = directory("results/genomes/")
    params:
        prefix = "results/genomes/refgenome.fasta"
    log:
        "results/logs/hisat_build/build.log"
    threads: 4
    conda:
        "../envs/env.yaml"
    wrapper:
        "v2.2.1/bio/hisat2/index"


rule HISAT:
    input:
        genome = "results/genomes/",
        s1 = "results/sortmerna/{sample}_fwd.fq.gz",
        s2 = "results/sortmerna/{sample}_rev.fq.gz"
    output:
        sam = temp("results/hisat/{sample}.sam")
    log:
        "results/logs/hisat/{sample}.log"
    threads: 4
    params: config["hisat2"]
    conda:
        "../envs/env.yaml"
    shell:
        """
        hisat2 {params} -x results/genomes/refgenome.fasta -p {threads} -1 {input.s1} -2 {input.s2} -S {output.sam} &> {log} 
        """

rule convertsamtobam:
    input: "results/hisat/{sample}.sam"
    output: "results/hisat/{sample}.bam"
    log:
        "results/logs/hisat/bam/{sample}.log"
    threads: 4
    conda:
        "../envs/env.yaml"
    shell:
        "samtools view -Sb {input} > {output} -@ {threads} 2> {log}"

rule bam_sorted:
    input: "results/hisat/{sample}.bam"
    output:
        bam_sorted = "results/bam_sorted/{sample}.bam",
        bai = "results/bam_sorted/{sample}.bam.bai"
    threads: 4
    log: "results/logs/bam_sorted/{sample}.log"
    conda:
        "../envs/env.yaml"
    shell:
        """
        samtools sort -@ {threads} -o {output.bam_sorted} {input} &> {log} && \
        samtools index {output.bam_sorted} {output.bai} 
        """