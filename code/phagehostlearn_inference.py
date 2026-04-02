"""
PhageHostLearn inference script
Converted from phagehostlearn_inference.ipynb for HPC (Biowulf) usage.

Usage:
    python phagehostlearn_inference.py \
        --phages_fasta /path/to/all_phages.fasta \
        --bacteria_list /path/to/bacteria_paths.txt \
        --output_path /path/to/output \
        --kaptive_db /path/to/Klebsiella_k_locus_primary_reference.gbk \
        --phanotate_path /path/to/phanotate.py \
        --hmmer_path /path/to/hmmer \
        --suffix inference

Inputs:
    --phages_fasta:  A single multi-FASTA file with all phage genomes (separated by >).
                     Each record's header becomes the phage ID.
    --bacteria_list: A text file with one path per line, each pointing to an individual
                     host genome FASTA file. The filename (minus .fasta) becomes the strain ID.
"""

import argparse
import os
import pickle
import numpy as np
import pandas as pd
from Bio import SeqIO
from xgboost import XGBClassifier

import phagehostlearn_processing as phlp
import phagehostlearn_features as phlf


def split_phage_fasta(phages_fasta, output_dir):
    """Split a multi-FASTA file into individual FASTA files, one per phage."""
    os.makedirs(output_dir, exist_ok=True)
    count = 0
    for record in SeqIO.parse(phages_fasta, 'fasta'):
        out_file = os.path.join(output_dir, record.id + '.fasta')
        SeqIO.write(record, out_file, 'fasta')
        count += 1
    print(f'  Split {count} phage genomes into {output_dir}')
    return output_dir


def link_bacteria_genomes(bacteria_list, output_dir):
    """Symlink individual bacterial genome FASTAs into a single directory."""
    os.makedirs(output_dir, exist_ok=True)
    count = 0
    with open(bacteria_list, 'r') as f:
        for line in f:
            fasta_path = line.strip()
            if not fasta_path:
                continue
            basename = os.path.basename(fasta_path)
            # Ensure .fasta extension for downstream compatibility
            name, ext = os.path.splitext(basename)
            if ext != '.fasta':
                basename = name + '.fasta'
            link_path = os.path.join(output_dir, basename)
            if not os.path.exists(link_path):
                os.symlink(os.path.abspath(fasta_path), link_path)
            count += 1
    print(f'  Linked {count} bacterial genomes into {output_dir}')
    return output_dir


def parse_args():
    parser = argparse.ArgumentParser(description='PhageHostLearn inference pipeline')
    parser.add_argument('--phages_fasta', type=str, required=True,
                        help='Multi-FASTA file containing all phage genomes')
    parser.add_argument('--bacteria_list', type=str, required=True,
                        help='Text file with one host genome FASTA path per line')
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

    # Create output directory
    os.makedirs(path, exist_ok=True)

    # Prepare input directories from user's data formats
    print('Preparing input data...')
    phages_dir = os.path.join(path, 'phage_genomes')
    bacteria_dir = os.path.join(path, 'bacteria_genomes')
    split_phage_fasta(args.phages_fasta, phages_dir)
    link_bacteria_genomes(args.bacteria_list, bacteria_dir)

    # ---- Step 1: Data processing ----
    phage_genes_file = path + '/phage_genes' + suffix + '.csv'
    gene_embeddings_file = path + '/phage_protein_embeddings' + suffix + '.csv'
    rbpbase_file = path + '/RBPbase' + suffix + '.csv'
    locibase_file = path + '/Locibase' + suffix + '.json'
    rbp_embeddings_path = path + '/esm2_embeddings_rbp' + suffix + '.csv'
    loci_embeddings_path = path + '/esm2_embeddings_loci' + suffix + '.csv'

    if os.path.exists(phage_genes_file):
        print('Step 1/4: PHANOTATE output found, skipping...')
    else:
        print('Step 1/4: Running PHANOTATE for phage gene calling...')
        phlp.phanotate_processing(path, phages_dir, args.phanotate_path, data_suffix=suffix)

    if os.path.exists(gene_embeddings_file):
        print('Step 2/4: Protein embeddings found, skipping...')
    else:
        print('Step 2/4: Computing protein embeddings for RBP detection...')
        phlp.compute_protein_embeddings(path, data_suffix=suffix)

    if os.path.exists(rbpbase_file):
        print('Step 3/4: RBPbase found, skipping...')
    else:
        print('Step 3/4: Detecting phage RBPs...')
        phlp.phageRBPdetect(path, args.pfam_path, args.hmmer_path, args.xgb_rbpdetect_path,
                            gene_embeddings_file, data_suffix=suffix)

    if os.path.exists(locibase_file):
        print('Step 4/4: Locibase found, skipping...')
    else:
        print('Step 4/4: Processing bacterial genomes with Kaptive...')
        phlp.process_bacterial_genomes(path, bacteria_dir, args.kaptive_db, data_suffix=suffix)

    # ---- Step 2: Feature construction ----
    if os.path.exists(rbp_embeddings_path):
        print('ESM-2 RBP embeddings found, skipping...')
    else:
        print('Computing ESM-2 embeddings for RBPs...')
        phlf.compute_esm2_embeddings_rbp(path, data_suffix=suffix)

    if os.path.exists(loci_embeddings_path):
        print('ESM-2 loci embeddings found, skipping...')
    else:
        print('Computing ESM-2 embeddings for K-loci...')
        phlf.compute_esm2_embeddings_loci(path, data_suffix=suffix)

    print('Constructing feature matrices...')
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
