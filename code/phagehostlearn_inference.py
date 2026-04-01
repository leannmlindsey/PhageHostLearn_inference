"""
PhageHostLearn inference script
Converted from phagehostlearn_inference.ipynb for HPC (Biowulf) usage.

Usage:
    python phagehostlearn_inference.py \
        --phages_path /path/to/phage_genomes \
        --bacteria_path /path/to/bacteria_genomes \
        --output_path /path/to/output \
        --kaptive_db /path/to/Klebsiella_k_locus_primary_reference.gbk \
        --phanotate_path /path/to/phanotate.py \
        --hmmer_path /path/to/hmmer \
        --suffix inference
"""

import argparse
import pickle
import numpy as np
import pandas as pd
from xgboost import XGBClassifier

import phagehostlearn_processing as phlp
import phagehostlearn_features as phlf


def parse_args():
    parser = argparse.ArgumentParser(description='PhageHostLearn inference pipeline')
    parser.add_argument('--phages_path', type=str, required=True,
                        help='Path to folder containing phage genome FASTA files')
    parser.add_argument('--bacteria_path', type=str, required=True,
                        help='Path to folder containing bacterial genome FASTA files')
    parser.add_argument('--output_path', type=str, required=True,
                        help='Path to output directory (intermediate files and results written here)')
    parser.add_argument('--kaptive_db', type=str, required=True,
                        help='Path to Kaptive K-locus reference database (.gbk)')
    parser.add_argument('--phanotate_path', type=str, required=True,
                        help='Path to phanotate.py')
    parser.add_argument('--hmmer_path', type=str, default='/usr/local/bin',
                        help='Path to HMMER installation')
    parser.add_argument('--pfam_path', type=str, default='RBPdetect_phageRBPs.hmm',
                        help='Path to RBP HMM profile')
    parser.add_argument('--xgb_rbpdetect_path', type=str, default='RBPdetect_xgb_hmm.json',
                        help='Path to RBP detection XGBoost model')
    parser.add_argument('--xgb_model_path', type=str, default='phagehostlearn_esm2_xgb.json',
                        help='Path to trained PhageHostLearn XGBoost model')
    parser.add_argument('--suffix', type=str, default='inference',
                        help='Suffix for output files')
    return parser.parse_args()


def main():
    args = parse_args()

    path = args.output_path
    suffix = args.suffix

    # Create output directory if it doesn't exist
    import os
    os.makedirs(path, exist_ok=True)

    # ---- Step 1: Data processing ----
    print('Step 1/4: Running PHANOTATE for phage gene calling...')
    phlp.phanotate_processing(path, args.phages_path, args.phanotate_path, data_suffix=suffix)

    print('Step 2/4: Computing protein embeddings for RBP detection...')
    phlp.compute_protein_embeddings(path, data_suffix=suffix)

    print('Step 3/4: Detecting phage RBPs...')
    gene_embeddings_file = path + '/phage_protein_embeddings' + suffix + '.csv'
    phlp.phageRBPdetect(path, args.pfam_path, args.hmmer_path, args.xgb_rbpdetect_path,
                        gene_embeddings_file, data_suffix=suffix)

    print('Step 4/4: Processing bacterial genomes with Kaptive...')
    phlp.process_bacterial_genomes(path, args.bacteria_path, args.kaptive_db, data_suffix=suffix)

    # ---- Step 2: Feature construction ----
    print('Computing ESM-2 embeddings for RBPs...')
    phlf.compute_esm2_embeddings_rbp(path, data_suffix=suffix)

    print('Computing ESM-2 embeddings for K-loci...')
    phlf.compute_esm2_embeddings_loci(path, data_suffix=suffix)

    print('Constructing feature matrices...')
    rbp_embeddings_path = path + '/esm2_embeddings_rbp' + suffix + '.csv'
    loci_embeddings_path = path + '/esm2_embeddings_loci' + suffix + '.csv'
    features_esm2, groups_bact = phlf.construct_feature_matrices(
        path, suffix, loci_embeddings_path, rbp_embeddings_path, mode='test')

    # ---- Step 3: Predict and rank ----
    print('Loading model and making predictions...')
    xgb = XGBClassifier()
    xgb.load_model(args.xgb_model_path)
    scores_xgb = xgb.predict_proba(features_esm2)[:, 1]

    # Build interaction score matrix
    groups_bact = np.asarray(groups_bact)
    loci_embeddings = pd.read_csv(loci_embeddings_path)
    rbp_embeddings = pd.read_csv(rbp_embeddings_path)
    bacteria = list(loci_embeddings['accession'])
    phages = list(set(rbp_embeddings['phage_ID']))

    score_matrix = np.zeros((len(bacteria), len(phages)))
    for i, group in enumerate(list(set(groups_bact))):
        scores_this_group = scores_xgb[groups_bact == group]
        score_matrix[i, :] = scores_this_group
    results = pd.DataFrame(score_matrix, index=bacteria, columns=phages)
    results.to_csv(path + '/prediction_results' + suffix + '.csv')

    # Rank phages per bacterium
    ranked = {}
    for group in list(set(groups_bact)):
        scores_this_group = scores_xgb[groups_bact == group]
        ranked_phages = [(x, y) for y, x in sorted(zip(scores_this_group, phages), reverse=True)]
        ranked[bacteria[group]] = ranked_phages

    with open(path + '/ranked_results' + suffix + '.pickle', 'wb') as f:
        pickle.dump(ranked, f)

    # Print summary
    print('\n=== Results Summary ===')
    print(f'Bacteria: {len(bacteria)}, Phages: {len(phages)}')
    top = 5
    scores = np.zeros((len(ranked.keys()), top))
    for i, acc in enumerate(ranked.keys()):
        topscores = [round(y, 3) for (x, y) in ranked[acc]][:top]
        scores[i, :] = topscores
    summary = pd.DataFrame(scores, index=list(ranked.keys()),
                           columns=[f'top_{k+1}' for k in range(top)])
    print(summary)

    print(f'\nResults saved to:')
    print(f'  {path}/prediction_results{suffix}.csv')
    print(f'  {path}/ranked_results{suffix}.pickle')


if __name__ == '__main__':
    main()
