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
    parser.add_argument('--ref_path', type=str, default="/space/grp/rschwartz/rschwartz/eval-references/work/8f/3b1c4802884f8b1cd5f56121d59a15/whole_cortex.h5ad")
    parser.add_argument('--ref_keys', type=str, nargs='+', default=["subclass", "class", "family"])
    parser.add_argument('--n_folds', type=int, default=5, help='Number of folds for k-fold cross-validation')

    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

    
    
def main():
  SEED = 42
  random.seed(SEED)
  np.random.seed(SEED)
  scvi.settings.seed = SEED
  args = parse_arguments()

  ref_path = args.ref_path
  ref_keys = args.ref_keys

  ref = sc.read_h5ad(ref_path)
  ref_name = os.path.basename(ref_path).replace(".h5ad", "")
  granular_key = ref_keys[0]

  # Remove cells with missing labels for the granular key
  ref = ref[~pd.isnull(ref.obs[granular_key])].copy()

  # Get all unique dataset titles
  dataset_titles = ref.obs["dataset_title"].unique()

  all_prob_dfs = []
  all_obs_dfs = []

  for held_out_title in dataset_titles:
      # Split reference into train and test by dataset_title
      test_mask = ref.obs["dataset_title"] == held_out_title
      train_mask = ~test_mask

      if train_mask.sum() == 0 or test_mask.sum() == 0:
          continue  # Skip if no train or test samples

      X_train = ref.obsm["scvi"][train_mask]
      y_train = ref.obs[granular_key][train_mask].values
      X_test = ref.obsm["scvi"][test_mask]
      y_test = ref.obs[granular_key][test_mask].values

      ref_train_obs = ref.obs.iloc[np.array(train_mask)].copy()
      ref_test_obs = ref.obs.iloc[np.array(test_mask)].copy()
      
      
      ref_train = ad.AnnData(X_train, obs=ref_train_obs)
      ref_train.obsm["scvi"] = X_train
      ref_test = ad.AnnData(X_test, obs=ref_test_obs)
      ref_test.obsm["scvi"] = X_test
      probs = rfc_pred(ref=ref_train, query=ref_test, ref_keys=[granular_key], seed=SEED)
      probabilities = probs[granular_key]['probabilities']
      class_labels = probs[granular_key]['class_labels']
      prob_df = pd.DataFrame(probabilities, columns=class_labels, index=ref_test.obs.index)
      obs_df = ref_test.obs.copy()
      obs_df["held_out_dataset_title"] = held_out_title
      
      held_out_title = (held_out_title.replace(" ", "_")
            .replace("/", "_")
            .replace("'", "_")
            .replace(",", "_")
            .replace("(", "")
            .replace(")", "")
            .replace("&", "")

      )
      
      prob_df.to_csv(f"{held_out_title}_loocv.prob.df.tsv", sep="\t", index=False)
      obs_df.to_csv(f"{held_out_title}_loocv.obs.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()

