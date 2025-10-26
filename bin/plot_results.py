#!/user/bin/python3

import argparse 
from pathlib import Path
import os
import sys
import numpy as np
import pandas as pd
import warnings
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

def parse_arguments():
    parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
    parser.add_argument('--aggregated_results', type=str,
                        help="files containing cell_type aggregate f1 results with params", default="/space/grp/rschwartz/rschwartz/eval-references/work/a8/554eda2d1cec78350123f3caa24bbd/all_metrics.tsv")
    parser.add_argument('--color_mapping_file', type=str, help="Mapping file", default = "/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/color_mapping.tsv")
    parser.add_argument('--mapping_file', type=str, help="Mapping file", default = "/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_human.tsv")
    parser.add_argument('--ref_keys', type=str, help="Reference keys to plot", default = ['subclass', 'class', 'family'], nargs="+")
    known_args, _ = parser.parse_known_args()
    return known_args



def make_stable_colors(color_mapping_df):
    
    all_subclasses = sorted(color_mapping_df["subclass"])
    # i need to hardcode a separate color palette based on the mmus mapping file
    # Generate unique colors for each subclass
    color_palette = sns.color_palette("husl", n_colors=len(all_subclasses))
    subclass_colors = dict(zip(all_subclasses, color_palette))
    return subclass_colors
    
    
def plot_line(df, x, y, hue, col, style, title, xcell_type, ycell_type, save_path):
    # set global fontsize
    plt.rcParams.update({'font.size': 14})  # Set default font size for all plot elements
    
    all_levels = df[hue].unique()
    
    colors = {subclass: color for subclass, color in zip(all_levels, sns.hls_palette(len(all_levels)))}

    # change figure size
    plt.figure(figsize=(22, 10))
    g = sns.relplot(
        data=df, x=x, y=y,
        col=col, hue=hue, style=style, palette=colors,
        kind="line", height=4, aspect=0.75, legend="full"  # Adjust figure size
    )
   # title = ""
    #title.replace("_", " ").title()  # Capitalize and substitute "_" with " " 
    g.figure.suptitle("", y=1, x = 0.5)  # Add title above plots
    g.set_axis_cell_types(xcell_type, ycell_type)  # Set axis cell_types

    g.set(xticks=[0,0.25,0.5,0.75])
    # set xtick fontsize
    # Rotate x-tick cell_types for better visibility
    for ax in g.axes.flat:
        plt.setp(ax.get_xtickcell_types(), rotation=45, ha="right", fontsize =12)  # Rotate 45 degrees and align
    
    plt.xcell_type(xcell_type)
    plt.ycell_type(ycell_type)
    # make legend fontzie smaller
    plt.setp(g._legend.get_texts(), fontsize=12)  # Adjust legend font size
    
   # plt.tight_layout()
    plt.savefig(save_path, bbox_inches="tight")


def plot_score_by_celltype(aggregated_results, levels, color_mapping_df, mapping_df, 
                           outdir="cutoff_plots", level="family", score_col="f1_score", subclass_col = "subclass"):
    os.makedirs(outdir, exist_ok=True)
    new_outdir = os.path.join(outdir, subclass_col)
    os.makedirs(new_outdir, exist_ok=True)
    plt.rcParams.update({'font.size': 20})

    all_subclasses = sorted(levels[subclass_col])
    subclass_colors = make_stable_colors(color_mapping_df)

    celltypes = levels[level]
    # remove "unknown" if remove_unknown is set
    if "unknown" in celltypes:
        celltypes = [ct for ct in celltypes if ct != "unknown"]
    methods = sorted(aggregated_results["method"].unique())

    rows, cols = len(celltypes), len(methods)
    fig, axes = plt.subplots(rows, cols, 
                             figsize=(cols * 6, rows * 5), squeeze=False, sharex=True, sharey=True)

    for i, celltype in enumerate(celltypes):
        group_subclasses = mapping_df[mapping_df[level] == celltype][subclass_col].unique()
        subclasses_to_plot = [subclass for subclass in all_subclasses if subclass in group_subclasses]
        if len(subclasses_to_plot) == 0:
           subclasses_to_plot = [celltype]

        filtered_df = aggregated_results[(aggregated_results["label"].isin(subclasses_to_plot)) & (aggregated_results["key"] == subclass_col)]
        
        if filtered_df.empty:
            print(f"No data available for cell type '{celltype}' with subclasses {subclasses_to_plot}. Skipping this cell type.")
            continue
        
        
        for j, method in enumerate(methods):
            ax = axes[i, j]
            method_df = filtered_df[filtered_df["method"] == method]
            
            sns.lineplot(data=method_df, x="cutoff", y=score_col, hue="label", 
                         marker="o", 
                         palette={cell_type: subclass_colors[cell_type] for cell_type in subclasses_to_plot}, ax=ax)
            if i == 0:
                ax.set_title(method)

            # Add legend only for the first subplot in each row
            if j == len(methods) - 1:
                ax.legend(title="label", bbox_to_anchor=(1, 0.5), loc="center left", fontsize=14)
            else:
                ax.legend_.remove()  # Remove legend from other subplots in the row
    for ax in axes.flat:
        ax.set_xticks([0, 0.25, 0.5, 0.75])
        ax.set_xticklabels([0, 0.25, 0.5, 0.75])
        ax.set_ylim(-0.05, 1.05)
        ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    plt.suptitle(f"{score_col.replace('_', ' ').title()} vs. Cutoff Across Cell Types", y=1)
    fig.tight_layout(rect=[0, 0, 1, 1])  # Adjust layout to fit legends

    # Save single figure
    save_path = os.path.join(new_outdir, f"{level}_{score_col}.png")
    plt.savefig(save_path, bbox_inches="tight")
    plt.close()



def main():
  # Parse command line arguments
  args = parse_arguments()
  aggregated_results_path = args.aggregated_results
  aggregated_results = pd.read_csv(aggregated_results_path, sep="\t", low_memory=False)
  color_mapping_df = pd.read_csv(args.color_mapping_file, sep="\t", low_memory=False)
  mapping_df = pd.read_csv(args.mapping_file, sep="\t", low_memory=False)
  ref_keys = args.ref_keys
  print(f"Reference keys to plot: {ref_keys}")
  # Split by reference and plot individually
  for ref_name, ref_df in aggregated_results.groupby("reference"):
      print(f"Plotting for reference: {ref_name}")
      subclasses = ref_df[ref_df["key"] == "subclass"]["label"].unique()
      classes = ref_df[ref_df["key"] == "class"]["label"].unique()
      families = ref_df[ref_df["key"] == "family"]["label"].unique()
      globalss = ref_df[ref_df["key"] == "global"]["label"].unique()

      levels = {
          "subclass": subclasses,
          "class": classes,
          "family": families,
          "global": globalss
      }

      for metric in ["f1_score", "precision", "recall"]:
          for level in ["subclass", "class", "family"]:
              plot_score_by_celltype(
                  ref_df, levels, color_mapping_df, mapping_df,
                  outdir=f"cutoff_plots/{ref_name}",
                  level=ref_keys[-1], score_col=metric, subclass_col=level
              )

if __name__ == "__main__":
    main()
