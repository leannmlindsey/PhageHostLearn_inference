#!/bin/bash
#
# PhageHostLearn inference runner
# Modify the variables below for your environment and data paths.
#

# ---- Data paths ----
DATA_PATH="./data"
PHAGES_PATH="./data/phage_genomes"
BACTERIA_PATH="./data/bacteria_genomes"
KAPTIVE_DB="./data/Klebsiella_k_locus_primary_reference.gbk"

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
    --data_path "$DATA_PATH" \
    --phages_path "$PHAGES_PATH" \
    --bacteria_path "$BACTERIA_PATH" \
    --kaptive_db "$KAPTIVE_DB" \
    --phanotate_path "$PHANOTATE_PATH" \
    --hmmer_path "$HMMER_PATH" \
    --pfam_path "$PFAM_PATH" \
    --xgb_rbpdetect_path "$XGB_RBPDETECT_PATH" \
    --xgb_model_path "$XGB_MODEL_PATH" \
    --suffix "$SUFFIX"
