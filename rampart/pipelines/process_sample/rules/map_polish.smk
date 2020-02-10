rule files:
    params:
        ref=config["output_path"] + "/binned_{sample}/{analysis_stem}.fasta",
        reads=config["output_path"]+"/binned_{sample}/{analysis_stem}.fastq"

rule minimap2_racon0:
    input:
        reads=rules.files.params.reads,
        ref=rules.files.params.ref
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/mapped.paf"
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule racon1:
    input:
        reads=rules.files.params.reads,
        fasta=rules.files.params.ref,
        paf= rules.minimap2_racon0.output
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon1.fasta"
    shell:
        "racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output}"

rule mafft1:
    input:
       fasta = rules.racon1.output,
       ref = rules.files.params.ref
    params:
        temp_file = config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/temp.racon1.fasta"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon1.aln.fasta"
    shell:
        "cat {input.ref} {input.fasta} > {params.temp_file} && "
        "mafft {params.temp_file} > {output} && "
        "rm {params.temp_file}"

rule clean1:
    input:
        aln = rules.mafft1.output,
        cns = rules.racon1.output
    params:
        path_to_script = workflow.current_basedir,
        seq_name = "{analysis_stem}"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon1.clean.fasta"
    shell:
        "python {params.path_to_script}/clean.py "
        "--alignment_with_ref {input.aln} "
        "--name {params.seq_name} "
        "--output_seq {output} "
        "--polish_round 1"

rule minimap2_racon1:
    input:
        reads=rules.files.params.reads,
        ref= rules.clean1.output
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/mapped.racon1.paf"
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule racon2:
    input:
        reads=rules.files.params.reads,
        fasta= rules.clean1.output,
        paf= rules.minimap2_racon1.output
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon2.fasta"
    shell:
        "racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output}"


rule mafft2:
    input:
       fasta = rules.racon2.output,
       ref = rules.files.params.ref
    params:
        temp_file = config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/temp.racon2.fasta"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon2.aln.fasta"
    shell:
        "cat {input.ref} {input.fasta} > {params.temp_file} && "
        "mafft {params.temp_file} > {output} && "
        "rm {params.temp_file}"

rule clean2:
    input:
        aln = rules.mafft2.output,
        cns = rules.racon2.output
    params:
        path_to_script = workflow.current_basedir,
        seq_name = "{analysis_stem}"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon2.clean.fasta"
    shell:
        "python {params.path_to_script}/clean.py "
        "--alignment_with_ref {input.aln} "
        "--name {params.seq_name} "
        "--output_seq {output} "
        "--polish_round 2"

rule minimap2_racon2:
    input:
        reads=rules.files.params.reads,
        ref= rules.clean2.output
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/mapped.racon2.paf"
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule racon3:
    input:
        reads=rules.files.params.reads,
        fasta= rules.clean2.output,
        paf= rules.minimap2_racon2.output
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon3.fasta"
    shell:
        "racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output}"

rule mafft3:
    input:
       fasta = rules.racon3.output,
       ref = rules.files.params.ref
    params:
        temp_file = config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/temp.racon3.fasta"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon3.aln.fasta"
    shell:
        "cat {input.ref} {input.fasta} > {params.temp_file} && "
        "mafft {params.temp_file} > {output} && "
        "rm {params.temp_file}"

rule clean3:
    input:
        aln = rules.mafft3.output,
        cns = rules.racon3.output
    params:
        path_to_script = workflow.current_basedir,
        seq_name = "{analysis_stem}"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon3.clean.fasta"
    shell:
        "python {params.path_to_script}/clean.py "
        "--alignment_with_ref {input.aln} "
        "--name {params.seq_name} "
        "--output_seq {output} "
        "--polish_round 3"


rule minimap2_racon3:
    input:
        reads=rules.files.params.reads,
        ref= rules.clean3.output
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/mapped.racon3.paf"
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule racon4:
    input:
        reads=rules.files.params.reads,
        fasta= rules.clean3.output,
        paf= rules.minimap2_racon3.output
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon4.fasta"
    shell:
        "racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output}"

rule mafft4:
    input:
       fasta = rules.racon4.output,
       ref = rules.files.params.ref
    params:
        temp_file = config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/temp.racon4.fasta"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon4.aln.fasta"
    shell:
        "cat {input.ref} {input.fasta} > {params.temp_file} && "
        "mafft {params.temp_file} > {output} && "
        "rm {params.temp_file}"

rule clean4:
    input:
        aln = rules.mafft4.output,
        cns = rules.racon4.output
    params:
        path_to_script = workflow.current_basedir,
        seq_name = "{analysis_stem}"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon4.clean.fasta"
    shell:
        "python {params.path_to_script}/clean.py "
        "--alignment_with_ref {input.aln} "
        "--name {params.seq_name} "
        "--output_seq {output} "
        "--polish_round 4"

rule minimap2_racon4:
    input:
        reads=rules.files.params.reads,
        ref= rules.clean4.output
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/mapped.racon4.paf"
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule medaka:
    input:
        basecalls=rules.files.params.reads,
        draft= rules.clean4.output
    params:
        outdir=config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}"
    output:
        config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus.fasta"
    threads:
        2
    shell:
        "medaka_consensus -i {input.basecalls} -d {input.draft} -o {params.outdir} -t 2 || touch {output}"

# rule mafft5:
#     input:
#        fasta = rules.medaka.output,
#        ref = rules.files.params.ref
#     params:
#         temp_file = config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/temp.medaka.fasta"
#     output:
#         config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/medaka.aln.fasta"
#     shell:
#         "cat {input.ref} {input.fasta} > {params.temp_file} && "
#         "mafft {params.temp_file} > {output} && "
#         "rm {params.temp_file}"

# rule clean5:
#     input:
#         aln = rules.mafft5.output,
#         cns = rules.medaka.output
#     params:
#         path_to_script = workflow.current_basedir,
#         seq_name = "{analysis_stem}"
#     output:
#         config["output_path"] + "/binned_{sample}/{analysis_stem}.consensus.fasta"
#     shell:
#         "python {params.path_to_script}/clean.py "
#         "--alignment_with_ref {input.aln} "
#         "--name {params.seq_name} "
#         "--output_seq {output} "
#         "--polish_round medaka"

rule gather_files:
    input:
        expand(config["output_path"] + "/binned_{{sample}}/medaka/{analysis_stem}/consensus.fasta", analysis_stem=config["analysis_stem"])
    params:
        sample = "{sample}"
    output:
        config["output_path"] + "/consensus_sequences/{sample}.fasta"
    run:
        with open(str(output[0]),"w") as fw:
            for fasta in input:
                for record in SeqIO.parse(str(fasta),"fasta"):
                    fw.write(f">{params.sample} {record.description}\n{record.seq}\n")

