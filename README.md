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
## Configuration file
The config file (config.yaml) can be found in the config folder, where every required input parameter and file should be specified. 
In addition, the user can set the mode of analysis (de-novo or reference-based), decide whether to do a metagenomic analysis using Kraken and set other settings like specific parameters for the different tools.
### Optional steps
Each of the optional steps can be set to True or False to indicate whether the particular step should be performed. In this workflow, the user can perform a reference-based mapping or a de novo guided reference mapping. When performing a de novo guided reference mapping, the user can also specify whether to perform decontamination and generate an kraken2 report based on which the genome will be downloaded for decontamination.
### Parameters
For most steps in the pipeline it is possible to adjust the parameters to specific needs. If this is possible for a tool, they can be added and changed in the config file under the tool's name. This includes software tools like Trim Galore where different adpater sequences for trimming can be specified, SortmeRNA, featureCount, samtools conensus and quast. Also for de novo-guided reference mapping different parameters can be changed.


