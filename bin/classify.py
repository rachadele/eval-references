
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
import yaml
from sklearn.metrics import precision_recall_curve, average_precision_score, PrecisionRecallDisplay
from collections import defaultdict
# silence warnings 
warnings.filterwarnings("ignore")


# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Classify cells given 1 ref and 1 query")
    parser.add_argument("--query_name", type=str, default = None)
    parser.add_argument('--obs', type=str, default = "")
    parser.add_argument('--ref_name', type=str, default = "Adult_mouse_cortical_cell_taxonomy_revealed_by_single_cell_transcriptomics")
    parser.add_argument('--ref_keys', type=str, nargs='+', default = ["subclass", "class", "family", "global"])
    parser.add_argument('--cutoff', type=float, help = "Cutoff threshold for positive classification")
    parser.add_argument('--probs', type=str, default="")
    parser.add_argument('--mapping_file', type=str, default="/space/grp/rschwartz/rschwartz/eval-references/meta/census_map_mouse_author.tsv")
    parser.add_argument('--ref_region_mapping', type=str)
    parser.add_argument('--use_gap', action='store_true', help="Use gap analysis for classification")
    parser.add_argument('--method', type=str, default="scvi") 
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
    ref_name = args.ref_name
    query_name = args.query_name
    ref_keys = args.ref_keys
    cutoff = args.cutoff
    ref_region_mapping = args.ref_region_mapping
    if args.use_gap:
        use_gap = True
    else:
        use_gap = False 
    method = args.method


    probs = pd.read_csv(args.probs, sep="\t")
    mapping_df = pd.read_csv(args.mapping_file, sep="\t")
    obs = pd.read_csv(args.obs, sep="\t")


    # Classify cells and evaluate
    obs = classify_cells(query=obs, ref_keys=ref_keys, cutoff=cutoff, probabilities=probs, mapping_df=mapping_df, use_gap=use_gap)

    outdir = os.path.join("predicted_meta")
    os.makedirs(outdir, exist_ok=True)

    # map valid labels for given query granularity and evaluate
    # remove this for cross validation
    obs = map_valid_labels(obs, ref_keys, mapping_df)
    class_metrics = evaluate_sample_predictions(obs, ref_keys, mapping_df)

    obs.to_csv(os.path.join(outdir,f"{ref_name}.predictions.{cutoff}.tsv"), index=False, sep="\t")

    #class_metrics = update_classification_report(class_metrics, ref_keys)

    # Plot confusion matrices
    for key in ref_keys:
        outdir = os.path.join("confusion")
        os.makedirs(outdir, exist_ok=True)
        plot_confusion_matrix(ref_name, ref_name, key, class_metrics[key]["confusion"], output_dir=outdir)

    # Collect F1 scores
    performance_metrics = []
    for key in ref_keys:
        label_metrics = class_metrics[key]["label_metrics"]
        weighted_metrics = class_metrics[key]["weighted_metrics"]
        macro_metrics = class_metrics[key]["macro_metrics"]
        overall_accuracy = class_metrics[key]["overall_accuracy"]
        for label, metrics in label_metrics.items():
            performance_metrics.append({
                'reference': ref_name,
                'label': label,
                'f1_score': metrics['f1_score'],
                'accuracy': metrics['accuracy'],
                'precision': metrics['precision'],
                'recall': metrics['recall'],
                'support': metrics['support'],
                'weighted_f1': weighted_metrics.get('f1_score', None),
                'weighted_precision': weighted_metrics.get('precision', None),
                'weighted_recall': weighted_metrics.get('recall', None),
                'macro_f1': macro_metrics.get('f1_score', None),
                'macro_precision': macro_metrics.get('precision', None),
                'macro_recall': macro_metrics.get('recall', None),
                'overall_accuracy': overall_accuracy,
                'key': key,
                'cutoff': cutoff,
                'method': method
                })

    # Save F1 scores to a file
    df = pd.DataFrame(performance_metrics)
    
    outdir = "label_transfer_metrics"
    os.makedirs(outdir, exist_ok=True)
    df.to_csv(os.path.join(outdir, f"{ref_name}_{method}_{cutoff}.summary.scores.tsv"), sep="\t", index=False)
    
if __name__ == "__main__":
    main()
    

