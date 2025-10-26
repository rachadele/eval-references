import argparse
import os
import random
from weakref import ref
import numpy as np
import pandas as pd
import scvi
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
import scanpy as sc

# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Classify cells given 1 ref and 1 query")
    parser.add_argument('--ref_name', type=str, default = "whole_cortex")
    parser.add_argument('--scvi_h5ad', type=str, help = "Path to h5ad reference", default="/space/grp/rschwartz/rschwartz/eval-references/2024-07-01/mus_musculus/dataset_id/SCT/gap_false/ref_50/refs/scvi_integrated/whole_cortex.h5ad")
    parser.add_argument('--sct_h5ad', type=str, help = "Path to h5ad reference with sctransform integrated PCA embedding", default="/space/grp/rschwartz/rschwartz/eval-references/work/60/5ea3209f06182b9c9d2729d7c3dcae/whole_cortex_sct_integrated.h5ad")
    parser.add_argument('--batch_fields', nargs="+", help="Fields to use for batch key", default=["dataset_title", "tissue", "assay"])
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args
    

def explore_umap(adata, embedding_key, ref_name, seed=42, resolution=1.0, n_iterations=2, ncols=2):
    """
    Compute neighbors, clusters, and UMAP for a given embedding and plot the result.
    """
    sc.pp.neighbors(adata, use_rep=embedding_key)
    leiden_key = f"leiden_{embedding_key}"
    sc.tl.leiden(adata, resolution=resolution, key_added=leiden_key, flavor="igraph", n_iterations=n_iterations)
    sc.tl.umap(adata, random_state=seed)
    sc.pl.umap(adata, color=["subclass", leiden_key, "batch"], ncols=ncols, save=f"_{ref_name}_{embedding_key}_umap.png")

def get_scib_benchmarks(adata, batch_key="batch", label_key="cell_type", embedding_keys=["scvi"], pre_integrated_embedding="Unintegrated"):
  bm = Benchmarker(adata,
    batch_key=batch_key,
    label_key=label_key,
    bio_conservation_metrics=BioConservation(isolated_labels=True, 
                                            nmi_ari_cluster_labels_leiden=True,
                                            nmi_ari_cluster_labels_kmeans=True),
    batch_correction_metrics=BatchCorrection(bras=False),
    pre_integrated_embedding_obsm_key = pre_integrated_embedding,
    embedding_obsm_keys=embedding_keys,
    n_jobs=6
  )
  bm.benchmark()
  bm.plot_results_table(save_dir=".") # save figure
  results_df = bm.get_results()
  return results_df
  
  


def main():
  args = parse_arguments()
  SEED = 42
  random.seed(SEED)         # For `random`
  np.random.seed(SEED)      # For `numpy`
  # For `torch`'
  scvi.settings.seed = SEED # For `scvi`
  batch_fields = args.batch_fields
  # get scib metrics for the given reference
  ref_name = args.ref_name
  scvi_adata = sc.read_h5ad(args.scvi_h5ad) # load scvi integrated adata
  sct_adata = sc.read_h5ad(args.sct_h5ad) # load sct integrated adata
  
  # take intersection of cells
  intersection = set(scvi_adata.obs_names).intersection(set(sct_adata.obs_names))
  
  order = list(intersection)
  scvi_adata = scvi_adata[order].copy()
  sct_adata = sct_adata[order].copy()
  
  # add sct integrated pca to scvi adata
  # this wasn't working
  new_adata = scvi_adata.copy()
  new_adata.obsm["X_pca_sct"] = sct_adata.obsm["X_pca"]

 # create a new batch key that combines dataset_title, tissue, and assay
  
  for cat in batch_fields:
      new_adata.obs[cat] = new_adata.obs[cat].astype(str)
       
  new_adata.obs["batch"] = new_adata.obs[batch_fields].agg("_".join, axis=1)
  #new_adata.obs["batch"] = new_adata.obs["dataset_title"] + "_" + new_adata.obs["tissue"] + "_" + new_adata.obs["assay"]



  explore_umap(
      new_adata,
      embedding_key="X_pca_sct",
      ref_name=ref_name,
      seed=SEED,
      ncols=1
  )
  
  explore_umap(
    new_adata,
    embedding_key="scvi",
    ref_name=ref_name,
    seed=SEED,
    ncols=1
  )

 
# get scib benchmarks for sct PCA and scvi latent embeddings
# not sure if i actually need to normalize the data and compute HGV first? it seems like the benchmarker class does that automatically

  #sc.pp.normalize_total(new_adata, target_sum=1e4) 
  #sc.pp.log1p(new_adata)
  sc.pp.highly_variable_genes(
    new_adata,
    n_top_genes=2000,
    flavor="cell_ranger",
    batch_key="batch",
    subset=False,
    inplace=True
)

  # compute unintegrated pca for pre integrated embedding
  sc.pp.pca(new_adata, n_comps=30, use_highly_variable=True)
  new_adata = new_adata[:, new_adata.var.highly_variable].copy() # not sure why this is necessary
  new_adata.obsm["Unintegrated"] = new_adata.obsm["X_pca"]

  explore_umap(
      new_adata,
      embedding_key="Unintegrated",
      ref_name=ref_name,
      seed=SEED,
      ncols=1
  )

 
  results_df = get_scib_benchmarks(
    new_adata,
    batch_key="batch",
    label_key="subclass",
    embedding_keys=["scvi", "X_pca_sct", "Unintegrated"],
    pre_integrated_embedding="X_pca" # can be replaced with "Unintegrated" ? unsure how this changes anything
    
  )
  results_df.to_csv(f"scib_benchmark_results_{ref_name}.tsv", sep="\t")


if __name__ == "__main__":
    main()