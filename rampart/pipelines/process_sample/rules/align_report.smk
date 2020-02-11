# rule align_cns_to_ref:
#     input:
#        fasta = config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus.fasta",
#        ref = config["output_path"] + "/binned_{sample}/{analysis_stem}.fasta"
#     params:
#         temp_file = config["output_path"] + "/binned_{sample}/temp.cns_ref_aln.fasta"
#     output:
#         config["output_path"] + "/binned_{sample}/{analysis_stem}.consensus_to_ref.aln.fasta"
#     shell:
#         "cat {input.ref} {input.fasta} > {params.temp_file} && "
#         "mafft {params.temp_file} > {output} && "
#         "rm {params.temp_file}"

rule generate_report:
    input:
        cns = config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus.fasta",
        ref = config["output_path"] + "/binned_{sample}/{analysis_stem}.fasta"
    params:
        path_to_script = workflow.current_basedir,
        sample = "{sample}_{analysis_stem}"
    output:
        config["output_path"] + "/binned_{sample}/report/{analysis_stem}.report.md"
    shell:
        """
        python {params.path_to_script}/make_report.py \
        -i {input.cns} \
        -r {input.ref} \
        -o {output} \
        --sample {params.sample}
        """

rule gather_reports:
    input:
        expand(config["output_path"] + "/binned_{{sample}}/report/{analysis_stem}.report.md", analysis_stem=config["analysis_stem"])
    output:
        config["output_path"] + "/reports/{sample}.report.md"
    shell:
        "cat {input} > {output}"
