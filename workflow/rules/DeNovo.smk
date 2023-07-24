rule fastqc_raw_wgs:
    input:
        s1=lambda wildcards: config["wgsReads"][wildcards.wgs_sample]
    output:
        html = "results/quality/raw_wgs/{wgs_sample}_fastqc.html",
        zip = "results/quality/raw_wgs/{wgs_sample}_fastqc.zip"
    log:
        "results/logs/fastqc_wgs/{wgs_sample}.log"
    conda:
        "../envs/env.yaml"
    threads: 4
    wrapper:
        "v1.31.1/bio/fastqc"

rule multiqc_raw_wgs:
    input:
        expand("results/quality/raw_wgs/{wgs_sample}_fastqc.zip",wgs_sample=["wgs1","wgs2"])
    output:
        "results/quality/raw_wgs/multiqc_report.html"
    log: "results/logs/multi_wgs/raw.log"
    conda:
        "../envs/env.yaml"
    shell:
        "multiqc -o results/quality/raw_wgs results/quality/raw_wgs &> {log}"

#creation of reference genome
# rename files because of trim_galore -> not able to change outfiles
rule rename_files:
    input:
        config["wgsReads"]["wgs1"],
        config["wgsReads"]["wgs2"]
    output:
        temp("results/trimgalore_wgs/wgs1.fastq.gz"),
        temp("results/trimgalore_wgs/wgs2.fastq.gz")
    shell:
        """
        mkdir -p results/trimgalore_wgs  # Create the directory if it doesn't exist
        cp {input[0]} {output[0]}
        cp {input[1]} {output[1]}
        """

rule wgs_trimming:
    input:
        "results/trimgalore_wgs/wgs1.fastq.gz",
        "results/trimgalore_wgs/wgs2.fastq.gz"
    output:
        temp("results/trimgalore_wgs/wgs1_val_1.fq.gz"),
        temp("results/trimgalore_wgs/wgs2_val_2.fq.gz"),
        temp("results/trimgalore_wgs/wgs1.fastq.gz_trimming_report.txt"),
        temp("results/trimgalore_wgs/wgs2.fastq.gz_trimming_report.txt")
    log:
        "results/logs/WGS_deNovo/trimgalore.log"
    threads: 4
    params:
        extra = config["trimgalore_wgs"]
    conda:
        "../envs/env.yaml"
    wrapper:
        "v1.31.1/bio/trim_galore/pe"

rule fastqc_trimming_wgs:
    input: "results/trimgalore_wgs/wgs{index}_val_{index}.fq.gz"
    output:
        html = "results/quality/trim_wgs/wgs{index}_fastqc.html",
        zip = "results/quality/trim_wgs/wgs{index}_fastqc.zip"
    log:
        "results/logs/WGS_deNovo/wgs{index}_fastqc_trim.log"
    conda:
        "../envs/env.yaml"
    threads: 4
    wrapper:
        "v1.31.1/bio/fastqc"

rule multiqc_trimming_wgs:
    input:
        zip1 = expand("results/quality/trim_wgs/wgs{index}_fastqc.zip",index=["1","2"])
    output:
        "results/quality/trim_wgs/multiqc_report.html"
    log: "results/logs/WGS_deNovo/multiqc_trim.log"
    conda:
        "../envs/env.yaml"
    shell:
        "multiqc -o results/quality/trim_wgs results/quality/trim_wgs &> {log}"

rule assembly:
    input:
        fastq1 = "results/decontaminated/genomeseq_1_trim.fastq" if config["Decontamination"]["decontaminate"] else "results/trimgalore_wgs/wgs1_val_1.fq.gz",
        fastq2 = "results/decontaminated/genomeseq_2_trim.fastq" if config["Decontamination"]["decontaminate"] else "results/trimgalore_wgs/wgs2_val_2.fq.gz"
    output:
        mega = directory("results/assembly/genome"),
        meganame = "results/assembly/genome/final.contigs.fa"
    conda:
        "../envs/env.yaml"
    threads: 4 # for 4 my laptop is not able to execute
    params:
        min_count = config["megahit"]["min-count"], #mincount = minimum multiplicity for filtering (k_min+1)-mers, default 2
        kmin = config["megahit"]["kmin"], #kmin = minimum kmer size (<= 127), must be odd number, default 21
        kmax = config["megahit"]["kmax"], #kmax = maximum kmer size (<= 127), must be odd number, default 99
        kstep = config["megahit"]["kstep"], #kstep = increment of kmer size of each iteration (<= 28), must be even number, default 20
        mincontiglength = config["megahit"]["mincontiglength"], #minimum length of contigs to output
        other = config["megahit"]["other"]
    log:
        "results/logs/WGS_deNovo/assembly.log"
    shell:
        """
        megahit -1 {input.fastq1} -2 {input.fastq2} -o {output.mega} --min-count {params.min_count} \
        --k-min {params.kmin} --k-max {params.kmax} --k-step {params.kstep} --min-contig-len {params.mincontiglength} \
        -f {params.other} -t {threads} &> {log}
        """

## creation of annotation genome

rule prokka:
    input:
        genome = "results/assembly/genome/final.contigs.fa"
    output:
        annotation_file = "results/ref_denovo/strain.gff"
    params:
        mincontig = config["prokka"]["mincontig"], # Minimum contig size [NCBI needs 200] (default '1')
        evalue = config["prokka"]["evalue"], #Similarity e-value cut-off (default '1e-06')
        other = config["prokka"]["other"] #e.g. --fast, --norrna, --notrna, --rawproduct "https://github.com/tseemann/prokka"
    log:
        "results/logs/WGS_deNovo/prokka/.log"
    threads: 4
    conda:
        "../envs/env.yaml"
    shell:
        """
        prokka {input.genome} --force --outdir results/ref_denovo/ --prefix strain \
        --cpus {threads} --mincontiglen {params.mincontig} --evalue {params.evalue} {params.other} &> {log}
        """

rule quast:
    input:
        genome = "results/assembly/genome/final.contigs.fa",
        annotation_file = "results/ref_denovo/strain.gff"
    output:
        dir = directory("results/quality/quast"),
        file= "results/quality/quast/report.tsv"
    threads: 4
    conda:
        "../envs/env.yaml"
    params:
        m = config["quast"]["min_contig"], #Lower threshold for contig length
        l = config["quast"]["minalignmentlength"], #The minimum alignment length
        other= config["quast"]["other"] #additional parameters
    log:
        "results/logs/WGS_deNovo/quast.log"
    shell:
        """
        quast.py {input.genome} -o {output.dir} -g {input.annotation_file} -l {params.l} -m {params.m} \
        {params.other} -t {threads} &> {log}
        """

## mapping of reads to de novo genome
rule hisat_build_denovo:
    input:
        fasta = "results/assembly/genome/final.contigs.fa"
    output:
        genome = directory("results/genomes_denovo")
    params:
        prefix = "results/genomes_denovo/refgenome.fasta"
    log:
            "results/logs/WGS_deNovo/hisat_build_denovo.log"
    conda:
        "../envs/env.yaml"
    threads:
        4
    wrapper:
        "v2.2.1/bio/hisat2/index"

rule HISAT_denovo:
    input:
        genome = "results/genomes_denovo/",
        s1 = "results/sortmerna/{sample}_fwd.fq.gz",
        s2 = "results/sortmerna/{sample}_rev.fq.gz"
    output:
        sam=temp("results/hisat_denovo/{sample}.sam")
    log:
        "results/logs/hisat_denovo/{sample}.log"
    conda:
        "../envs/env.yaml"
    threads:
        4
    params: config["hisat2"]
    shell:
        "hisat2 {params} -p {threads} -x results/genomes_denovo/refgenome.fasta -1 {input.s1} -2 {input.s2} -S {output.sam} &> {log}"


rule convertsamtobam_denovo:
    input:
        "results/hisat_denovo/{sample}.sam"
    output:
        "results/hisat_denovo/{sample}.bam"
    log:
        "results/logs/convertsamtobam_denovo/{sample}.log"
    conda:
        "../envs/env.yaml"
    threads:
        4
    shell:
        "samtools view -Sb {input} --threads {threads} &> {output} 2> {log}"


rule sortbam_denovo:
    input:
        "results/hisat_denovo/{sample}.bam"
    output:
        bam=temp("results/bam_sorted_denovo/{sample}.bam"),
        bai=temp("results/bam_sorted_denovo/{sample}.bam.bai")
    log:
        sort="results/logs/sortbam_denovo/{sample}.log",
        index="results/logs/indexbam_denovo/{sample}.log"
    threads: 4
    conda:
        "../envs/env.yaml"
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input} &> {log.sort} && \
        samtools index {output.bam} {output.bai} &> {log.index}
        """


