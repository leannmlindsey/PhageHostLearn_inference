#!/bin/bash
#
# PhageHostLearn inference runner
# Modify the variables below for your environment and data paths.
#

# ---- Input data paths ----
PHAGES_FASTA="/path/to/all_phages.fasta"        # Multi-FASTA with all phage genomes
BACTERIA_LIST="/path/to/bacteria_paths.txt"      # One host genome FASTA path per line
KAPTIVE_DB="/path/to/Klebsiella_k_locus_primary_reference.gbk"

# ---- Output path (intermediate files and results go here) ----
OUTPUT_PATH="/path/to/output"

# ---- Software paths ----
PHANOTATE_PATH="/path/to/phanotate.py"
HMMER_PATH="/path/to/hmmer"

# ---- Model files (defaults are in the code directory) ----
PFAM_PATH="RBPdetect_phageRBPs.hmm"
XGB_RBPDETECT_PATH="RBPdetect_xgb_hmm.json"
XGB_MODEL_PATH="phagehostlearn_esm2_xgb.json"

# ---- Output suffix ----
SUFFIX="inference"

# ---- Run ----
python phagehostlearn_inference.py \
    --phages_fasta "$PHAGES_FASTA" \
    --bacteria_list "$BACTERIA_LIST" \
    --output_path "$OUTPUT_PATH" \
    --kaptive_db "$KAPTIVE_DB" \
    --phanotate_path "$PHANOTATE_PATH" \
    --hmmer_path "$HMMER_PATH" \
    --pfam_path "$PFAM_PATH" \
    --xgb_rbpdetect_path "$XGB_RBPDETECT_PATH" \
    --xgb_model_path "$XGB_MODEL_PATH" \
    --suffix "$SUFFIX"
