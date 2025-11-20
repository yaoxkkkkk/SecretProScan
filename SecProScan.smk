import os
import subprocess

configfile: "PhyPredi_config.yaml"

# Function to check if a file is gzipped
def is_gzipped(file):
    return file.endswith(".gz")

# Function to get the base name without extensions
def get_basename(file):
    if is_gzipped(file):
        file = os.path.splitext(file)[0]  # Remove .gz
    return os.path.splitext(os.path.basename(file))[0]

# Check if protein_file is provided
if "protein_file" in config and config["protein_file"]:
    ref_basename = get_basename(config["protein_file"])
    protein_file = config["protein_file"]
elif "genome_file" in config and config["genome_file"] and "annotation_file" in config and config["annotation_file"]:
    ref_basename = get_basename(config["genome_file"])
    genome_file = config["genome_file"]
    annotation_file = config["annotation_file"]
else:
    raise ValueError("Either protein_file or both genome_file and annotation_file must be provided in the config file.")

rule all:
    input:
        f"results/{ref_basename}_phytocytokine.txt",
        f"results/{ref_basename}_phytocytokine.fasta"

# If protein_file is provided, use the following rules
if "protein_file" in config and config["protein_file"]:
    rule file_treatment:
        input:
            protein_file=config["protein_file"]
        output:
            temp(f"{ref_basename}_treated.fasta")
        log:
            "logs/file_treatment.log"
        shell:
            r"""
            seqkit replace \
            -p "\s.+" \
            {input.protein_file} \
            | sed 's/\*$//' \
            1> {output} \
            2> {log}
            """

    rule length_filter:
        input:
            protein_file=f"{ref_basename}_treated.fasta"
        output:
            f"{ref_basename}_shortpep.fasta"
        params:
            max_length=config["max_length"]
        log:
            "logs/length_filter.log"
        shell:
            """
            seqkit seq \
            -M {params.max_length} \
            {input.protein_file} \
            1> {output} \
            2> {log}
            """

# If genome_file and annotation_file are provided, use the following rules
elif "genome_file" in config and config["genome_file"] and "annotation_file" in config and config["annotation_file"]:
    rule decompress_files:
        input:
            genome_file=config["genome_file"],
            annotation_file=config["annotation_file"]
        output:
            genome_file=temp(f"{ref_basename}.genome"),
            annotation_file=temp(f"{ref_basename}.annotation")
        run:
            if is_gzipped(input.genome_file):
                shell("gunzip -c {input.genome_file} > {output.genome_file}")
            else:
                shell("cp {input.genome_file} {output.genome_file}")
            if is_gzipped(input.annotation_file):
                shell("gunzip -c {input.annotation_file} > {output.annotation_file}")
            else:
                shell("cp {input.annotation_file} {output.annotation_file}")

    rule primary_transcripts_extract:
        input:
            annotation_file=f"{ref_basename}.annotation"
        output:
            primary_only_file=temp(f"{ref_basename}.primaryTranscriptOnly.gff3")
        container:
            config["singularity"]["AGAT_image"]
        log:
            "logs/primary_transcripts_extract.log"
        shell:
            """
            agat_sp_keep_longest_isoform.pl \
            -gff {input.annotation_file} \
            -o {output} \
            2> {log}
            """

    rule primary_sequence_extract:
        input:
            genome_file=f"{ref_basename}.genome",
            annotation_file=f"{ref_basename}.primaryTranscriptOnly.gff3"
        output:
            protein_file=f"{ref_basename}.protein_primaryTranscriptOnly.fasta"
        container:
            config["singularity"]["AGAT_image"]
        log:
            "logs/primary_sequence_extract.log"
        shell:
            """
            agat_sp_extract_sequences.pl \
            -g {input.annotation_file} \
            -f {input.genome_file} \
            -t cds \
            -p \
            --cfs \
            --cis \
            -o {output} \
            2> {log}
            """

    rule length_filter:
        input:
            protein_file=f"{ref_basename}.protein_primaryTranscriptOnly.fasta"
        output:
            f"{ref_basename}_shortpep.fasta"
        params:
            max_length=config["max_length"]
        log:
            "logs/length_filter.log"
        shell:
            """
            seqkit seq \
            -M {params.max_length} \
            {input.protein_file} \
            1> {output} \
            2> {log}
            """

# Common rules for both cases
rule Prediction_Phobius:
    input:
        protein_file=f"{ref_basename}_shortpep.fasta"
    output:
        f"results/Phobius/{ref_basename}_Phobius.txt"
    log:
        "logs/Prediction_Phobius.log"
    shell:
        """
        phobius.pl \
        -short {input.protein_file} \
        | sed '1s/SEQENCE //' \
        | awk 'BEGIN {{OFS="\\t"}} {{print $1, $2, $3, $4}}' \
        1> {output} \
        2> {log}
        """

rule NoTM_shortpep_extraction:
    input:
        protein_file=f"{ref_basename}_shortpep.fasta",
        Phobius_file=f"results/Phobius/{ref_basename}_Phobius.txt"
    output:
        NoTM_shortpep_list=temp(f"{ref_basename}_shortpep_NoTM.list"),
        NoTM_shortprotein_file=f"{ref_basename}_shortpep_NoTM.fasta"
    shell:
        """
        awk -F "\\t" '$2 == "0" {{print $1}}' {input.Phobius_file} > {output.NoTM_shortpep_list}

        seqkit grep -f {output.NoTM_shortpep_list} {input.protein_file} 1> {output.NoTM_shortprotein_file}
        """

rule Prediction_signalP:
    input:
        protein_file=f"{ref_basename}_shortpep_NoTM.fasta"
    output:
        "results/signalP/output.gff3"
    params:
        "results/signalP/"
    conda:
        config["conda_env"]["signalP_env"]
    log:
        "logs/Prediction_signalP.log"
    shell:
        """
        signalp6 \
        --fastafile {input} \
        --output_dir {params} \
        --format txt \
        --organism eukarya \
        --mode slow-sequential \
        2> {log}
        """

rule Prediction_Predisi:
    input:
        protein_file=f"{ref_basename}_shortpep_NoTM.fasta"
    output:
        f"results/Predisi/{ref_basename}_Predisi.txt"
    params:
        predisi_folder=config["predisi_folder"]
    log:
        "logs/Prediction_Predisi.log"
    shell:
        """
        java -cp {params.predisi_folder} \
        JSPP {params.predisi_folder}/matrices/eukarya.smx {input} {output} \
        2> {log}
        """

rule Phytocytokine_candidate_extraction:
    input:
        Phobius_file=f"results/Phobius/{ref_basename}_Phobius.txt",
        signalP_file="results/signalP/output.gff3",
        Predisi_file=f"results/Predisi/{ref_basename}_Predisi.txt"
    output:
        Phobius_candidate=f"results/Phobius/{ref_basename}_Phobius_candidate.txt",
        signalP_candidate=f"results/signalP/{ref_basename}_signalP_candidate.txt",
        Predisi_candidate=f"results/Predisi/{ref_basename}_Predisi_candidate.txt"
    shell:
        """
        awk -F "\\t" '$3 == "Y" {{print $1}}' {input.Phobius_file} > {output.Phobius_candidate}
        awk -F "\\t" '$3 == "signal_peptide" {{print $1}}' {input.signalP_file} > {output.signalP_candidate}
        awk -F "\\t" '$3 == "Y" {{print $1}}' {input.Predisi_file} > {output.Predisi_candidate}
        """

rule Extract_common_candidates:
    input:
        Phobius_candidate=f"results/Phobius/{ref_basename}_Phobius_candidate.txt",
        signalP_candidate=f"results/signalP/{ref_basename}_signalP_candidate.txt",
        Predisi_candidate=f"results/Predisi/{ref_basename}_Predisi_candidate.txt"
    output:
        common_candidates=f"results/{ref_basename}_phytocytokine.txt"
    shell:
        """
        awk 'FNR==NR {{a[$1]++; next}} {{a[$1]++}} END {{for (i in a) if (a[i] >= 2) print i}}' {input.Phobius_candidate} {input.signalP_candidate} {input.Predisi_candidate} > {output.common_candidates}
        """

rule Phytocytokine_sequence_extraction:
    input:
        common_candidates=f"results/{ref_basename}_phytocytokine.txt",
        protein_file=f"{ref_basename}_shortpep_NoTM.fasta"
    output:
        f"results/{ref_basename}_phytocytokine.fasta"
    shell:
        """
        seqkit grep -f {input.common_candidates} {input.protein_file} 1> {output}
        """