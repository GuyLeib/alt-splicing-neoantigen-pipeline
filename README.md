# alt-splicing-neoantigen-pipeline

This repository contains a set of Python scripts designed to identify alternative splicing-derived neoantigens from RNA-seq data. The pipeline leverages SplAdder for splice graph construction and novel event identification, custom scripts for processing SplAdder output and translating peptides, and pVACbind for neoantigen prediction.

# Overview
The pipeline consists of three main stages, executed by individual Python scripts:

run_spladder.py: Builds splice graphs and identifies alternative splicing events using SplAdder.

create_splicing_antigens.py: Processes SplAdder's confirmed events, filters based on gene expression, reconstructs spliced protein sequences, and generates FASTA files of candidate neoantigen peptides.

run_pvacbind_all_events.py: Performs MHC binding and Immunogenicity predictions on the generated peptides using pVACbind (via Docker).

