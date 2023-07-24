# Sustainable Workflow for Shigella flexneri RNA analysis
Final project for the course Applied Sequence Analysis 2023
Inga Tummoszeit and Georg von Arnim
24/07/2023
------------
## Project outline
This work presents an automated workflow for analyzing bacterial RNA-Seq data of distinct organisms, focusing on differential gene expression and phylogenetic tree construction. The workflow incorporates paired-end sequencing data processing, quality control, and grouping into two categories. It offers both reference-based analysis using annotated genomes and reference-free analysis using de novo assembly from whole genome sequencing reads. The workflow is scalable, automated with Snakemake, and adheres to established standards of reproducibility, adaptability, and transparency.
## Dependencies and execution
To use the pipeline, conda and snakemake must be installed. Additional packages for the various bioinformatics tools are automatically downloaded into conda. 
The workflow can be executed in the SeqAn_2023_Group4 folder with: snakemake --use-conda -s workflow/Snakefile --core 6. The number of cores can be changed accordingly. Then all files in rule all in Workflow/Snakefile are created.
-----------------------------------

