
rule kraken2_report:
    input:
        krakendb = config["Decontamination"]["krakendatabase"],
        trim1 = "results/trimgalore_wgs/wgs1_val_1.fq.gz",
        trim2 = "results/trimgalore_wgs/wgs2_val_2.fq.gz",
    output:
        "results/kraken2/kraken2.report"
    log:
        "results/logs/kraken2/krakenreport.log"
    conda:
        "../envs/env.yaml"
    threads: 1
    shell:
        "kraken2 --db {input.krakendb} --paired {input.trim1} {input.trim2} --threads {threads} --report {output} > /dev/null 2> {log}"

rule multiqc_kraken:
    input:
        "results/kraken2/kraken2.report"
    output:
        "results/kraken2/multiqc_report.html"
    log:
        "results/logs/kraken2/multiqc.log"
    conda:
        "../envs/env.yaml"
    shell:
        "multiqc -o results/kraken2 results/kraken2 &> {log}"


rule trim_index:
    input:
        genome = config["Decontamination"]["contaminationgenome"],
    output:
        str(config["Decontamination"]["contaminationgenome"] + ".bwt")
    conda:
        "../envs/env.yaml"
    log:
        "results/logs/trim_index/index.log"
    shell:
        "bwa index {input.genome} &> {log}"


rule decontaminate_map:
    input:
        genome = config["Decontamination"]["contaminationgenome"],
        index = str(config["Decontamination"]["contaminationgenome"] + ".bwt"),
        trim1 = "results/trimgalore_wgs/wgs1_val_1.fq.gz",
        trim2 = "results/trimgalore_wgs/wgs2_val_2.fq.gz",
    output:
        temp("results/decon_mapped/genomeseq.sam")
    conda:
        "../envs/env.yaml"
    params:
        B = config["bwa-mem"]["B"],  #mismatch penalty
        k = config["bwa-mem"]["k"],  #minimum seed length
        w = config["bwa-mem"]["w"], #band width for banded alignment [100]
        d = config["bwa-mem"]["d"], #off-diagonal X-dropoff [100]
        r = config["bwa-mem"]["r"], #look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]
        y = config["bwa-mem"]["y"] ,#seed occurrence for the 3rd round seeding [20]
        c = config["bwa-mem"]["c"], #skip seeds with more than INT occurrences [500]
        D = config["bwa-mem"]["D"], #drop chains shorter than FLOAT fraction of the longest overlapping chain [0.50]
        W = config["bwa-mem"]["W"], #discard a chain if seeded bases shorter than INT [0]
        m = config["bwa-mem"]["m"], #perform at most INT rounds of mate rescues for each read [50]
        other = config["bwa-mem"]["other"]
    threads:
        4
    log:
        "results/logs/decon_map/genome.log"
    shell:
        """
        bwa mem -t {threads} -k {params.k} -B {params.B} -w {params.w} -d {params.d} -r {params.r} -y {params.y} \
        -c {params.c} -D {params.D} -W {params.W} -m {params.m} {params.other} {input.genome} {input.trim1} {input.trim2} > {output} 2> {log}
        """


rule contaminatedsamtobam:
    input:
        "results/decon_mapped/genomeseq.sam"
    output:
        "results/decon_mapped/genomeseq.bam"
    log:
        "results/logs/convertsamtobam/genomeseq.log"
    conda:
        "../envs/env.yaml"
    threads:
        4
    shell:
        "samtools view -Sb {input} --threads {threads} &> {output} 2> {log}"

rule sortcontamination:
    input:
        contamination = "results/decon_mapped/genomeseq.bam",
    output:
        cont_sorted = "results/decon_mapped_sorted/genomeseq.bam"
    conda:
        "../envs/env.yaml"
    log:
        "results/logs/sort_contamination/genomeseq.log",
        "results/logs/index_contamination/genomeseq.log"
    threads: 4
    shell:
        """
        samtools sort -@ {threads} -n {input.contamination} -o {output.cont_sorted} &> {log[0]}
        """

rule decontaminate_filter:
    input:
        contamination = "results/decon_mapped_sorted/genomeseq.bam",
    output:
        bam_unmapped = temp("results/decon_unmapped_bam/genomeseq.bam")
    conda:
        "../envs/env.yaml"
    log:
        "results/logs/decontaminate_filter/genomeseq.log"
    threads:
        4
    shell:
        "samtools view -b -f 12 -@ {threads} {input.contamination} > {output.bam_unmapped} 2> {log}"

rule decontaminate_create_fastq:
    input:
        contamination = "results/decon_unmapped_bam/genomeseq.bam"
    output:
        fasta1 = "results/decontaminated/genomeseq_1_trim.fastq",
        fasta2= "results/decontaminated/genomeseq_2_trim.fastq"
    conda:
        "../envs/env.yaml"
    log:
        "results/logs/decontaminate_create_fastq/genomeseq.log"
    threads: 4
    shell:
        "samtools fastq -@ {threads} -1 {output.fasta1} -2 {output.fasta2} {input.contamination} &> {log}"

rule fastqc_decon_wgs:
    input: "results/decontaminated/genomeseq_{number}_trim.fastq"
    output:
        html = "results/quality/decon_wgs/genomeseq_{number}_trim_fastqc.html",
        zip = "results/quality/decon_wgs/genomeseq_{number}_trim_fastqc.zip"
    log:
        "results/logs/fastqc_wgs/genomeseq_{number}.log"
    conda:
        "../envs/env.yaml"
    threads: 4
    wrapper:
        "v1.31.1/bio/fastqc"

rule multiqc_decon_wgs:
    input:
        zip = expand("results/quality/decon_wgs/genomeseq_{number}_trim_fastqc.zip", number=["1", "2"])
    output:
        "results/quality/decon_wgs/multiqc_report.html"
    log: "results/logs/multi_wgs/decon.log"
    conda:
        "../envs/env.yaml"
    shell:
        "multiqc -o results/quality/decon_wgs results/quality/decon_wgs &> {log}"
