
rule featureCount:
    input:
        bam = "results/bam_sorted/{sample}.bam" if config["ref_genome_present"] else "results/bam_sorted_denovo/{sample}.bam",
        gff = config["ref_GFF"] if config["ref_genome_present"] else "results/ref_denovo/strain.gff"
    output:
        counts = "results/featureCount/{sample}_counts.txt"
    threads: 4
    log:
        "results/logs/featureCount/{sample}.log"
    params:
        d = config["featureCount"]["Min_fragment"], # Minimum fragment/template length, 50 by default
        D = config["featureCount"]["Max_fragment"], # Maximum fragment/template length, 600 by default.
        Q = config["featureCount"]["Min_map_score"], # Minimum mapping quality score a read must satisfy, 0 by default
        fracOverlap = config["featureCount"]["FracOverlap"], # Minimum fraction overlapping bases in a reads that required for read assigment, float between 0 and 1, default 0
        fracFeature = config["featureCount"]["FracOverFeature"], # Minimum fraction of bases included in a feature that is required to overlap with a read. Value within range [0,1]. 0 by default.
        maxMOp = config["featureCount"]["MaxMOp"], # Specify the maximum number of ‘M’ operations (matches or mis-matches) allowed in a CIGAR string, default 10
        minOverlap = config["featureCount"]["MinOverlap"], # Minimum number of overlapping bases in a read that is required for read assignment. 1 by default.
        readEx3 = config["featureCount"]["ReadExtension3"], # Reads are extended downstream by < int > bases from their 3’ end. 0 by default.
        readEx5 = config["featureCount"]["ReadExtension5"], # Reads are extended upstream by < int > bases from their 5’ end. 0 by default
        other = config["featureCount"]["other"] # e.g. -B, -C, -f, -j etc
    conda:
        "../envs/env.yaml"
    shell:
        """
        featureCounts -a {input.gff} -F GFF -t CDS -g ID -p -T {threads} -d {params.d} -D {params.D} -Q {params.Q} \
         --fracOverlap {params.fracOverlap} --fracOverlapFeature {params.fracFeature} --maxMOp {params.maxMOp} --minOverlap {params.minOverlap} --readExtension3 {params.readEx3} --readExtension5 {params.readEx5} {params.other} -o {output.counts} {input.bam} &> {log}
        """

