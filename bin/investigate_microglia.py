#!/usr/bin/python3

import warnings
warnings.filterwarnings("ignore")

from pathlib import Path
import os
import sys
import random
import argparse
import json

import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import cellxgene_census
import cellxgene_census.experimental
import scvi
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import yaml

import utils
from utils import *
from sklearn.model_selection import train_test_split

# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
    parser.add_argument('--organism', type=str, default='mus_musculus', help='Organism name (e.g., homo_sapiens)')
    parser.add_argument('--census_version', type=str, default='2024-07-01', help='Census version (e.g., 2024-07-01)')
    parser.add_argument('--subsample_ref', type=int, default=50)
    parser.add_argument('--h5ad', type=str, default="/space/grp/rschwartz/rschwartz/eval-references/2024-07-01/mus_musculus/dataset_id/SCT/gap_false/ref_500/refs/scvi_integrated/whole_cortex.h5ad", help='Path to the h5ad file')
    parser.add_argument('--resolution', type=float, default=1.0, help='Resolution for Leiden clustering')
    parser.add_argument('--seed', type=int, default=42, help='Random seed for UMAP')
    parser.add_argument('--ref_keys', nargs='+', default=["subclass","class","family","global"], help='Reference keys for classification')
    parser.add_argument('--outdir', type=str, default="/space/grp/rschwartz/rschwartz/eval-references/microglia-investigate", help='Output directory for results')
    parser.add_argument('--mapping_file', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_mouse_author.tsv")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args
      

def rank_and_plot(ref_subset, reference, outdir, embedding_key="scvi", n_genes=25):
    """
    Run sc.tl.rank_genes_groups and plot/save the results for a given reference group.
    """
    sc.tl.rank_genes_groups(ref_subset, groupby="subclass", reference=reference, method="wilcoxon")
    plt.figure(figsize=(10, 6))
    sc.pl.rank_genes_groups(ref_subset, show=False, n_genes=n_genes, save=False)
    if reference == "rest":
        fname = "rest_vs_microglia_macrophage_de_genes.png"
    else:
        fname = f"{reference.lower()}_vs_{'microglia' if reference=='Macrophage' else 'macrophage'}_de_genes.png"
    plt.savefig(os.path.join(outdir, fname))
    plt.close()

def main():
    args = parse_arguments()
    h5ad_path = args.h5ad
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    ref_subset = sc.read_h5ad(h5ad_path)
    seed = args.seed
    resolution = args.resolution
    sc.pp.normalize_total(ref_subset, target_sum=1e4) 
    sc.pp.log1p(ref_subset)
    os.makedirs(outdir, exist_ok=True)
    ref_keys = args.ref_keys
    mapping_df = pd.read_csv(args.mapping_file, sep="\t")
 

    # restrict to .obs["family"]=="CNS macrophage"
    ref_subset = ref_subset[ref_subset.obs["family"] == "CNS macrophage"]
    # value counts of author cell type, subclass, class, family, global, dataset title, collection name
    ref_subset.obs[["author_cell_type", "subclass", "class", "family", "global", "dataset_title", "collection_name"]].value_counts().reset_index().to_csv(
        os.path.join(outdir, "microglia_cell_types.tsv"), sep="\t", index=False
    )
    embedding_key = "scvi"
    sc.pp.neighbors(ref_subset, use_rep=embedding_key)
    leiden_key = f"leiden_{embedding_key}"
    sc.tl.leiden(ref_subset, resolution=resolution, key_added=leiden_key, flavor="igraph", n_iterations=2)
    sc.tl.umap(ref_subset, random_state=seed)
    
    sc.pl.umap(
        ref_subset,
        color=["author_cell_type", "subclass", "collection_name"],
        ncols=1,
        save="_microglia_macrophage.png",
    )
    # Differential expression and plotting for each reference
    rank_and_plot(ref_subset, "Macrophage", outdir, embedding_key)
    rank_and_plot(ref_subset, "Microglia", outdir, embedding_key)
    rank_and_plot(ref_subset, "rest", outdir, embedding_key)
    
    
    # train test split on ref_subset
    train, test = train_test_split(ref_subset.obs, test_size=0.2, random_state=seed)
    train_subset = ref_subset[ref_subset.obs_names.isin(train.index)]
    test_subset = ref_subset[ref_subset.obs_names.isin(test.index)]
    probs = rfc_pred(train_subset, test_subset, ref_keys, seed)
    
    probabilities = probs[ref_keys[0]]['probabilities']
    class_labels = probs[ref_keys[0]]['class_labels']
    prob_df = pd.DataFrame(probabilities, columns=class_labels)
 
   
    obs = test_subset.obs.copy()
    obs = classify_cells(query=obs, ref_keys=ref_keys, cutoff=0, probabilities=prob_df, mapping_df=mapping_df, use_gap=False)
    
if __name__ == "__main__":
    main()