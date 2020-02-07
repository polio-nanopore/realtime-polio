rule trim_to_mapped_region:
    input:
        reads = config["output_path"]+"/binned_{sample}/{analysis_stem}.fastq",
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
            if record.id in coord_dict:
                start,end = coord_dict[record.id]
                record = record[start:end]
                records.append(record)
            else:
                print(f"{record.id} not in the paf file.")
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
        mapping_file = rules.minimap2_for_variant.output
    output:
        config["output_path"] + "/binned_{sample}/error_correction/{analysis_stem}/mapped.bam.bai"
    shell:
        "samtools index {input}"

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
        bcftools mpileup -Ou -a FORMAT/AD --max-depth 100000 -f {input.consensus} {input.mapping} | \
        bcftools call -vmO v --ploidy 1 -o {output}
        """

rule assess_variants:
    input:
        config["output_path"] + "/binned_{sample}/error_correction/{analysis_stem}/variants.vcf"
    output:
        config["output_path"] + "/binned_{sample}/error_correction/{analysis_stem}/variant_sequences.fasta"
    run:
        
    params:
        sample = "{sample}",
        output_path= config["output_path"],
        path = workflow.current_basedir
    output:
        cns = config["output_path"] + "/consensus_sequences/{sample}.fasta",
        reports = config["output_path"] + "/reports/{sample}.report.md",
        vars = config["output_path"] + "/binned_{sample}/variant_calls/{sample}.vcf"
    run:
        analysis_stem = stems.fetch(params.sample)
        print(params.sample, analysis_stem)

        if analysis_stem != "":
            print("Passing {} for {} into processing pipeline.".format(analysis_stem, params.sample))
            config["analysis_stem"]= analysis_stem
            shell("snakemake --nolock --snakefile {params.path}/../process_sample/Snakefile "
                        "--configfile {input.config} "
                        "--config "
                        "analysis_stem={config[analysis_stem]} "
                        "output_path={params.output_path} "
                        "sample={params.sample} "
                        "--rerun-incomplete")
        else:
            shell("touch {output.cns} && touch {output.reports}")

rule cat_variants:
    input:
        expand(config["output_path"] + "/binned_{{sample}}/error_correction/{analysis_stem}/variants.vcf", analysis_stem=config["analysis_stem"])
    output:
        config["output_path"] + "/binned_{sample}/variant_calls/{sample}.vcf"
    threads:
        2
    shell:
        "cat {input} > {output}"



