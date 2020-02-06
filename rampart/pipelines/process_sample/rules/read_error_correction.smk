
rule trim_to_mapped_region:
    input:
        reads = config["output_path"]+"/binned_{sample}/{analysis_stem}.fastq"
        mapping = config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/mapped.racon4.paf"
    output:
        config["output_path"]+"/binned_{sample}/{analysis_stem}.trimmed.fastq"
    run:
        coord_dict = {}
        with open(str(input.mapping),"r") as f:
            for l in f:
                l = l.rstrip('\n')
                tokens = l.split()
                start,end = tokens[2:4]
                coord_dict[tokens[0]]=(int(start),int(end))
                
        records = []
        for record in SeqIO.parse(str(input.reads), "fastq"):
            start,end = coord_dict[record.id]
            record = record[start:end]
            records.append(record)
        with open(output[0], "w") as outfile:
            SeqIO.write(records, outfile, "fastq")

rule minimap2_for_variant:
    input:
        reads=rules.trim_to_mapped_region.output,
        ref=config["output_path"] + "/binned_{sample}/{analysis_stem}.consensus.fasta"
    output:
        config["output_path"] + "/binned_{sample}/error_correction/{analysis_stem}/mapped.bam"
    shell:
        "minimap2 -ax map-ont {input.ref} {input.reads} | samtools sort -o {output}"

rule index_bam:
    input:
        reads = rules.trim_to_mapped_region.output,
        mapping_file = rules.minimap2_for_variant.output
    output:
        config["output_path"] + "/binned_{sample}/error_correction/{analysis_stem}/mapped.bam.bai"
    shell:
        """
        samtools index {input} -o {output} 
        """

rule index_cns:
    input:
        config["output_path"] + "/binned_{sample}/{analysis_stem}.consensus.fasta"
    output:
        config["output_path"] + "/binned_{sample}/{analysis_stem}.consensus.fasta.fai"
    shell:
        "samtools faidx {input}"

rule call_variants:
    input:
        consensus = config["output_path"] + "/binned_{sample}/{analysis_stem}.consensus.fasta",
        ref_index = config["output_path"] + "/binned_{sample}/{analysis_stem}.consensus.fasta.fai",
        mapping=config["output_path"] + "/binned_{sample}/error_correction/{analysis_stem}/mapped.bam",
        mapping_index = config["output_path"] + "/binned_{sample}/error_correction/{analysis_stem}/mapped.bam.bai"
    output:
        config["output_path"] + "/binned_{sample}/error_correction/{analysis_stem}/variants.vcf"
    threads:
        2
    shell:
        """
        bcftools mpileup -Ou --max-depth 10000 -f {input.consensus} {input.mapping} | \
        bcftools call -vmO v --ploidy 10 -o {output}
        """

rule call_variants:
    input:
        expand(config["output_path"] + "/binned_{{sample}}/error_correction/{analysis_stem}/variants.vcf", analysis_stem=config["analysis_stem"])
    output:
        config["output_path"] + "/binned_{sample}/variant_calls/{sample}.vcf"
    threads:
        2
    shell:
        "cat {input} > {output}"



