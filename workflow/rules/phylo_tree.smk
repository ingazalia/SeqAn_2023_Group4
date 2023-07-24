rule merged_all:
    input:
        merged_all_input,
        sample= config["samples"]
    output:
        "results/merged_alignment/alignment.fasta"
    conda:
        "../envs/env.yaml"
    script:
        "../scripts/merge_all.py"

rule iqtree:
    input:
        alignment = "results/merged_alignment/alignment.fasta"
    output:
        tree="results/merged_alignment/alignment.fasta.treefile",
        visualization="results/tree/tree.pdf"
    log:
        "results/logs/tree/iqtree.log",
        "results/logs/tree/vis.log"
    threads:
        4
    conda:
        "../envs/env.yaml"
    shell:
        "iqtree -s {input.alignment} -T {threads} &> {log[0]} && figtree -graphic PDF {output.tree} {output.visualization} &> {log[1]}"

