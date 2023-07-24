rule consensus:
    input:
        bam="results/bam_sorted/{sample}.bam" if config["ref_genome_present"] else "results/bam_sorted_denovo/{sample}.bam",
        bai="results/bam_sorted/{sample}.bam.bai" if config["ref_genome_present"] else "results/bam_sorted_denovo/{sample}.bam.bai"
    output:
        "results/consensus/{sample}.fasta"
    log:
        "results/logs/consensus/{sample}.log"
    conda:
        "../envs/env.yaml"
    params:
        depth = config["consensus"]["depth"],
        mode = config["consensus"]["mode"],# Select the consensus algorithm. "simple" frequency counting and the "bayesian" (Gap5) methods, with Bayesian being the default.
        min_MQ = config["consensus"]["min_MQ"], #Filters out reads with a mapping quality below INT. This defaults to zero.
        min_BQ = config["consensus"]["min_BQ"], #Filters out bases with a base quality below INT. This defaults to zero.
        other = config["consensus"]["other"]
    threads: 4
    shell:
        """
        samtools consensus -a -d {params.depth} -m {params.mode} --min-MQ {params.min_MQ} --min-BQ {params.min_BQ} \
        {params.other} -@ {threads} {input.bam} -o {output} &> {log}
         """

rule getfasta:
    input:
        genome = "results/consensus/{sample}.fasta",
        gff = config["ref_GFF"] if config["ref_genome_present"] else "results/ref_denovo/strain.gff"
    output:
        "results/cutgeneseqs/{sample}.fasta"
    log:
        "results/logs/cutgeneseqs/{sample}.log"
    params: config["getfasta"]
    conda:
        "../envs/env.yaml"
    threads:
        4
    shell:
        "bedtools getfasta -fi {input.genome} -bed {input.gff} -name {params} -fo {output} -s &> {log}"

checkpoint get_ref_genome: #checkpoint
    input:
        files=expand("results/cutgeneseqs/{sample}.fasta", sample=sample)
    params:
        threshold=config["threshold_N"]
    output:
        directory("results/genes")
    conda:
        "../envs/env.yaml"
    script:
        "../scripts/merge_genes.py"


def merged_all_input(wildcards):
    checkpoint_output = checkpoints.get_ref_genome.get(**wildcards).output[0]
    return expand("results/alignment/{gene}.fasta",
            gene=glob_wildcards(os.path.join(checkpoint_output,"{gene}.fasta")).gene)

rule alignment:
    input:
        "results/genes/{gene}.fasta"
    output:
        alignment="results/alignment/{gene}.fasta"
    log: "results/logs/alignment/{gene}.log"
    params:
        weight = config["mafft"]["weight"], # Weighting factor (valid of --globalpair, --localpair, --genafpair, --fastapair or --blastpair), default 2.7
        tree = config["mafft"]["tree"], # Guide tree, valid with 6mer distance, default 2
        maxiterate = config["mafft"]["maxiterate"], # number cycles, default 0
        other= config["mafft"]["other"]
    conda:
        "../envs/env.yaml"
    shell:
        """
            mafft --auto --weighti {params.weight} --retree {params.tree} --maxiterate {params.maxiterate} \
            {params.other} {input} > {output.alignment} 2> {log}
        """
