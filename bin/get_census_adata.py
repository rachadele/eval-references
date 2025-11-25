#!/user/bin/python3

import warnings
from weakref import ref
warnings.filterwarnings("ignore")
from pathlib import Path
import os
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import warnings
import utils
from utils import *
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import yaml


# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
    parser.add_argument('--organism', type=str, default='mus_musculus', help='Organism name (e.g., homo_sapiens)')
    parser.add_argument('--census_version', type=str, default='2024-07-01', help='Census version (e.g., 2024-07-01)')
    parser.add_argument('--subsample_ref', type=int, default=50)
    parser.add_argument('--relabel_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_mouse_author.tsv")
    parser.add_argument('--ref_collections', type=str, nargs = '+', default = [
        "A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation",
        "An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types",
        "Tabula Muris Senis"
    ]) 
    parser.add_argument('--split_column', type=str, default="dataset_id")
    parser.add_argument('--ref_keys', type=str, nargs="+", default=["subclass","class","family","global"])
    parser.add_argument('--seed', type=int, default=42)
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args
        
     
def create_ref_region_yaml(refs, outdir):
    ref_names = list(refs.keys())
    ref_region_yaml = {}
    for ref_name in ref_names:
        unique_regions = refs[ref_name].obs["tissue"].unique().tolist()
        
        # If there are multiple region types, handle them by including them as a list
        if len(unique_regions) == 1:
            ref_region_yaml[ref_name.replace(" ", "_").replace("\\/", "_") \
                .replace("(","").replace(")","") \
                .replace("\\", "") \
                .replace("'", "") \
                .replace(":", "") \
                .replace(";", "") \
                .replace("&", "") ] = unique_regions[0]
        else:
            ref_region_yaml[ref_name.replace(" ", "_").replace("\\/", "_") \
                .replace("(","").replace(")","") \
                .replace("\\", "") \
                .replace("'", "") \
                .replace(":", "") \
                .replace(";", "") \
                .replace("&", "") ] = "multiple regions"
    
    with open("ref_region.yaml", 'w') as file:
        yaml.dump(ref_region_yaml, file)

         
def main():
    # Parse command line arguments
    args = parse_arguments()

    # Set organism and census_version from arguments
    organism = args.organism
    census_version = args.census_version
    subsample_ref = args.subsample_ref
    relabel_path = args.relabel_path
    ref_collections = args.ref_collections
    split_column = args.split_column
    ref_keys = args.ref_keys
    SEED = args.seed

    if organism == "mus_musculus":
        original_celltypes = get_original_celltypes(columns_file=f"/space/grp/rschwartz/rschwartz/eval-references/meta/author_cell_annotations/{census_version}/original_celltype_columns.tsv",
                                                    author_annotations_path=f"/space/grp/rschwartz/rschwartz/eval-references/meta/author_cell_annotations/{census_version}") 
    else:
        original_celltypes = None

    refs = utils.get_census(
        organism=organism,
        organ="brain",
        subsample=subsample_ref,
        split_column=split_column,
        census_version=census_version,
        relabel_path=relabel_path,
        ref_collections=ref_collections,
        seed=SEED,
        ref_keys=ref_keys,
        original_celltypes = original_celltypes
    )
    
    refs.pop('All - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - Smart-seq2', None)
    refs.pop('Microglia - 24 months old wild-type and Rag1-KO', None)
    refs.pop('Single-cell of aged oligodendrocytes', None)

    print("finished fetching anndata")
    outdir = "refs"
    os.makedirs(outdir, exist_ok=True)

    for ref_name, ref in refs.items():
        if ref.shape[0] < 50:
            continue
        
        new_dataset_title = (
            ref_name
            .replace(" ", "_")
            .replace("/", "_")
            .replace("'", "_")
            .replace(",", "_")
            .replace("(", "")
            .replace(")", "")
            .replace("&", "")
        )

        ref.write(f"{new_dataset_title}.h5ad")
        ref.obs.to_csv(f"{new_dataset_title}.obs.tsv", sep="\t")



    create_ref_region_yaml(refs, outdir)

      
if __name__ == "__main__":
    main()