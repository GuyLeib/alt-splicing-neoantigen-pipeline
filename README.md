# alt-splicing-neoantigen-pipeline

This repository contains a set of Python scripts designed to identify alternative splicing-derived neoantigens from RNA-seq data. The pipeline leverages SplAdder for splice graph construction and novel event identification, custom scripts for processing SplAdder output and translating peptides, and pVACbind for neoantigen prediction.

# Overview
The pipeline consists of three main stages, executed by individual Python scripts:

run_spladder.py: Builds splice graphs and identifies alternative splicing events using SplAdder.

create_splicing_antigens.py: Processes SplAdder's confirmed events, filters based on gene expression, reconstructs spliced protein sequences, and generates FASTA files of candidate neoantigen peptides.

run_pvacbind_all_events.py: Performs MHC binding and Immunogenicity predictions on the generated peptides using pVACbind (via Docker).

How to run - examples: 

python run_spladder.py \
    -b /path/to/your/bam_files_directory \
    -o /path/to/your/output_results_directory/spladder_output \
    -g /path/to/your/reference.gtf \
    --workers 10 \

python create_splicing_antigens.py \
    -d /path/to/your/output_results_directory/spladder_output \
    -o /path/to/your/output_results_directory/splicing_peptides_fasta \
    -e /path/to/your/gene_expression.tsv.gz \
    -g /path/to/your/reference.gtf \
    -r /path/to/your/reference.fasta \
    --db-path /path/to/your/annotation.db \
    --min-expression 1.0 \
    --psi-threshold 0.5 \
    --event-types exon_skip alt_3prime alt_5prime mutex_exons intron_retention \
    --min-peptide-length 8 \
    --workers 8

python run_pvacbind_all_events.py \
    -i /path/to/your/output_results_directory/splicing_peptides_fasta \
    -o /path/to/your/output_results_directory/pvacbind_results \
    -a /path/to/your/hla_alleles.json \
    --docker-image griffithlab/pvactools \
    --iedb-install-dir /opt/iedb \
    --algorithms NetMHCpan NetMHC DeepImmuno \
    --workers 4 \
    --threads-per-sample 12 \
