source("/space/grp/rschwartz/rschwartz/eval-references/bin/seurat_functions.R")
library(reticulate)
use_condaenv("~/anaconda3/envs/r4.3/")
library(Seurat)
library(dplyr)
library(gplots)
library(argparse)
library(cellxgene.census)
library(data.tree)
set.seed(123)
library(caret)
options(future.globals.maxSize = 10 * 1024^3)  # 10 GB




# set max
parser <- ArgumentParser(description = "Process Seurat objects and transfer labels.")
parser$add_argument("--integration_method", type="character", help="Integration method of query and reference", default="pcaproject")
parser$add_argument("--ref_keys", type="character", nargs="*", help="List of reference keys to pass to query_transfer", default=c("subclass", "class","family","global"))
parser$add_argument("--dims", type="integer", help="Number of dimensions", default=50)
parser$add_argument("--max.features", type="integer", help="Maximum number of features", default=200)
parser$add_argument("--k.anchor", type="integer", help="Number of anchors", default=10)
parser$add_argument("--k.score", type="integer", help="?", default=30)
parser$add_argument("--cutoff", type="numeric", help="Cutoff threshold for label transfer prediction scores", default=0)
parser$add_argument("--ref_path", type="character", help="path to references", default="/space/grp/rschwartz/rschwartz/eval-references/work/08/e27fc2041e59458284cfd632a29f17/whole_cortex.rds")
parser$add_argument("--k.weight", type="integer", help="k.weight", default=50)
parser$add_argument("--n_folds", type="integer", help="Number of folds for cross-validation", default=5)
parser$add_argument("--normalization_method", type="character", help="Normalization method", default="SCT")
parser$add_argument("--nfeatures", type="integer", help="Number of variable features to use for dim reduction", default=2000)
parser$add_argument("--seed", type="integer", help="Random seed", default=42)
args <- parser$parse_args()
# Extract arguments from the parsed list    
#batch_key <- args$batch_key
integration_method <- args$integration_method
ref_keys <- args$ref_keys
dims <- args$dims
max.features <- args$max.features
k.anchor <- args$k.anchor
k.score <- args$k.score
project.query <- NULL
k.weight <- args$k.weight
ref_path <- args$ref_path
seed <- args$seed
normalization_method <- args$normalization_method
n_features <- args$nfeatures

ref = readRDS(ref_path)
ref_key <- ref_keys[1]

# Get all unique dataset_title values
# remove "whole cortex"
dataset_titles <- unique(ref@meta.data$dataset_title)
dataset_titles <- dataset_titles[dataset_titles != "whole_cortex"]

all_pred_scores <- list()
all_test_meta <- list()
ref_name <- basename(ref_path) %>% gsub(".rds", "", .)

for (held_out_title in dataset_titles) {
    test_idx <- which(ref@meta.data$dataset_title == held_out_title)
    train_idx <- setdiff(seq_len(ncol(ref)), test_idx)
    if (length(train_idx) == 0 || length(test_idx) == 0) next


    # subsetting is causing issues with integrated assays
    # ok i see the problem, the SCT assay in the integrated reference is from the dataset "Single-cell RNA-seq for all cortical & hippocampal regions (10x)"
    # so when we subset to remove that dataset, the SCT assay is empty
    # try recomputing SCT assay in both train and test
    train_obj <- subset(ref, cells = colnames(ref)[train_idx])
    test_obj <- subset(ref, cells = colnames(ref)[test_idx])

    if (normalization_method == "SCT") {
        train_obj <- SCTransform(train_obj, verbose = FALSE, variable.features.n = n_features)
        test_obj <- SCTransform(test_obj, verbose = FALSE, variable.features.n = n_features)
        train_obj <- RunPCA(train_obj, verbose = FALSE, npcs=dims)
        test_obj <- RunPCA(test_obj, verbose = FALSE, npcs=dims)
    }
    
    pred_scores <- transfer_multiple_labels(
        query = test_obj,
        reference = train_obj,
        reduction = integration_method,
        normalization_method = normalization_method,
        ref_keys = ref_keys,
        dims = dims,
        max.features = max.features,
        k.anchor = k.anchor,
        k.score = k.score,
        k.weight = k.weight
    )
    pred_scores <- as.data.frame(pred_scores)
    #pred_scores$held_out_dataset_title <- held_out_title
    #all_pred_scores[[held_out_title]] <- pred_scores
    test_meta <- test_obj@meta.data
    #test_meta$held_out_dataset_title <- held_out_title
    #all_test_meta[[held_out_title]] <- test_meta



    # Clean up column names as before
    pred_scores <- pred_scores %>%
      as.data.frame() %>%
      select(-any_of(c("predicted.id", "prediction.score.max"))) %>%
      rename_with(~ gsub("prediction.score.", "", .)) %>%
      rename_with(~ gsub("\\.", " ", .)) %>%
      rename_with(~ gsub("L([0-9]+) ([0-9]+) (.*)", "L\\1/\\2 \\3", .)) %>%
      rename_with(~ gsub("L2/3 6 IT", "L2/3-6 IT", .)) %>%
      rename_with(~ gsub("CA1 ProS", "CA1-ProS", .)) %>%
      rename_with(~ gsub("L4 RSP ACA", "L4 RSP-ACA", .)) %>%
      rename_with(~ gsub("CA2 IG FC", "CA2-IG-FC", .)) %>%
      rename_with(~ gsub("SUB ProS", "SUB-ProS", .)) %>%
      rename_with(~ gsub("L4 IT ET", "L4 IT/ET", .)) %>%
      rename_with(~ gsub("Cajal Retzius", "Cajal-Retzius", .))



    # write to file with name of hold out dataset

    # replace weird characters in held out title
    held_out_title = held_out_title %>% 
                        gsub(" ", "_", .) %>%
                        gsub("/", "_", .) %>%
                        gsub("'", "_", .) %>%
                        gsub(",", "_", .) %>%
                        gsub("\\(", "", .) %>%
                        gsub("\\)", "", .) %>%
                        gsub("&", "", .) 

    write.table(pred_scores, file=paste0(held_out_title,"_loocv.prediction.scores.tsv"), sep="\t", row.names=FALSE)
    write.table(test_meta, file=paste0(held_out_title,"_loocv.obs.tsv"), row.names=FALSE, sep= "\t")



}



