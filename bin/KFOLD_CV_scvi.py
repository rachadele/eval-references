#!/user/bin/python3

from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import cellxgene_census
import scvi
from scipy.sparse import csr_matrix
import warnings
import cellxgene_census
import cellxgene_census.experimental
import scvi
from sklearn.ensemble import RandomForestClassifier
import utils
from utils import *
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_curve, auc
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
import os
import json
from types import SimpleNamespace

# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Cross-validate classifier using k-fold cross-validation on reference data")
    parser.add_argument('--ref_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/refs/whole_cortex.h5ad")
    parser.add_argument('--ref_keys', type=str, nargs='+', default=["subclass", "class", "family"])
    parser.add_argument('--n_folds', type=int, default=5, help='Number of folds for k-fold cross-validation')

    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

    
    
def main():
    SEED = 42
    random.seed(SEED)         # For `random`
    np.random.seed(SEED)      # For `numpy`
    # For `torch`'
    scvi.settings.seed = SEED # For `scvi`
    # Parse command line arguments
    args = parse_arguments()

    # Set variables from arguments
   # organism = args.organism
   # census_version = args.census_version
    ref_path = args.ref_path
    ref_keys = args.ref_keys
    n_folds = args.n_folds

    # Load reference dataset
    ref = sc.read_h5ad(ref_path)
    ref_name = os.path.basename(ref_path).replace(".h5ad", "")

    # Use only the most granular key
    granular_key = ref_keys[0]

    # Remove cells with missing labels for the granular key
    ref = ref[~pd.isnull(ref.obs[granular_key])].copy()

    X = ref.obsm["scvi"]
    y = ref.obs[granular_key].values
    idx_all = np.arange(ref.n_obs)

    from sklearn.model_selection import StratifiedKFold
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=SEED)

    all_prob_dfs = []
    all_obs_dfs = []

    for fold, (train_idx, test_idx) in enumerate(skf.split(X, y)):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        ref_train = ad.AnnData(X_train, obs=ref.obs.iloc[train_idx].copy())
        ref_train.obsm["scvi"] = X_train
        ref_test = ad.AnnData(X_test, obs=ref.obs.iloc[test_idx].copy())
        ref_test.obsm["scvi"] = X_test

        probs = rfc_pred(ref=ref_train, query=ref_test, ref_keys=[granular_key], seed=SEED)
        probabilities = probs[granular_key]['probabilities']
        class_labels = probs[granular_key]['class_labels']
        prob_df = pd.DataFrame(probabilities, columns=class_labels, index=ref_test.obs.index)
        all_prob_dfs.append(prob_df)
        obs_df = ref_test.obs.copy()
        all_obs_dfs.append(obs_df)

    # Concatenate all folds
    prob_df_all = pd.concat(all_prob_dfs)
    obs_df_all = pd.concat(all_obs_dfs)

    outdir = "probs"
    os.makedirs(outdir, exist_ok=True)
    prob_df_all.to_csv(os.path.join(outdir, f"{ref_name}_kfold.prob.df.tsv"), sep="\t", index=False)
    obs_df_all.to_csv(f"{ref_name}_kfold.obs.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
    
