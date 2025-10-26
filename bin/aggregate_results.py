#!/user/bin/python3

from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import warnings
#import adata_functions
#from adata_functions import *
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
import os
import json
import ast
import sys
import matplotlib.lines as mlines
from utils import *


# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--metrics', type=str, nargs = "+", 
                        help="files containing all summary scores",)

    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args


def main():
    # Parse command line arguments
    args = parse_arguments()
    metrics_files = args.metrics
    
    # Initialize an empty DataFrame to store all metrics
    all_metrics = pd.DataFrame()
    # Loop through each file and read the metrics
    for file in metrics_files:
        # Read the metrics file
        metrics_df = pd.read_csv(file, sep="\t", low_memory=False)
        
        # Check if the DataFrame is empty
        if metrics_df.empty:
            print(f"Warning: {file} is empty. Skipping.")
            continue
        
        # Extract the prefix from the file name
        prefix = os.path.basename(file).split("_metrics.tsv")[0]
        
        # Add a column for the prefix
        metrics_df['prefix'] = prefix
        
        # Append to the all_metrics DataFrame
        all_metrics = pd.concat([all_metrics, metrics_df], ignore_index=True)


    # if support is 0, drop the row
    all_metrics.drop(all_metrics[all_metrics['support'] == 0].index, inplace=True)
    # write to tabular format
    all_metrics.to_csv("all_metrics.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
    
