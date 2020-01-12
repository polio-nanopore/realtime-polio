rule map_to_cns:
    input:
        basecalls=config["output_path"]+"/binned_{sample}/{analysis_stem}.fastq",
        consensus= config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus.fasta"
    output:
        config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus.bam"
    shell:
        "minimap2 -ax map-ont {input.ref} {input.reads} | samtools view -b - > {output}"

rule index_map_to_cns:
    input:
        config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus.bam"
    output:
        config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus.bam.bai"
    shell:
        "samtools index {input}"

rule call_cns_variants:
    input:
        consensus= config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus.fasta",
        map = config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus.bam",
        index = config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus.bam.bai"
    params:
        dir = config["output_path"] + "/binned_{sample}/medaka_variant/{analysis_stem}"
    output:
        config["output_path"] + "/binned_{sample}/medaka_variant/{analysis_stem}/round_1_phased.vcf"
    shell:
        "medaka_variant -f {input.consensus} -i {input.map} -o {params.dir} -m r941_trans || touch {output}"
