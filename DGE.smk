rule merge_counts:
    input:
        expand("results/featureCount/{sample}_counts.txt", sample=sample),
        sheet= config["samples_sheet"]
    output: "results/featureCount/merged.txt"
    script:
        "../scripts/featureCounts_merged.py"


rule diffexpression_R:
    input:
        "results/featureCount/merged.txt",
        config["samples_sheet"]
    output:
        "results/diffexpression/results_diffexpression.txt",
        "results/diffexpression/upregulatedgenes.txt",
        "results/diffexpression/upregulated_results.txt",
        "results/diffexpression/downregulationgenes.txt",
        "results/diffexpression/downregulated_results.txt",
        "results/diffexpression/MAplot.png",
        "results/diffexpression/pheatmap_down.png",
        "results/diffexpression/PCARNAseq.png"
    log: "results/logs/diffexpression/R.log"
    conda:
        "../envs/env.yaml"
    script:
        "../scripts/differentialexpression.R"
