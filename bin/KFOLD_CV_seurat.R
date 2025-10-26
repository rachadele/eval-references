
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

parser <- ArgumentParser(description = "Process Seurat objects and transfer labels.")
parser$add_argument("--integration_method", type="character", help="Integration method of query and reference", default="pcaproject")
parser$add_argument("--ref_keys", type="character", nargs="*", help="List of reference keys to pass to query_transfer", default=c("subclass", "class","family"))
parser$add_argument("--dims", type="integer", help="Number of dimensions", default=50)
parser$add_argument("--max.features", type="integer", help="Maximum number of features", default=200)
parser$add_argument("--k.anchor", type="integer", help="Number of anchors", default=10)
parser$add_argument("--k.score", type="integer", help="?", default=30)
parser$add_argument("--cutoff", type="numeric", help="Cutoff threshold for label transfer prediction scores", default=0)
parser$add_argument("--ref_path", type="character", help="path to references", default="")
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

# K-fold cross-validation using transfer_multiple_labels (moved to bottom)
n_folds <- args$n_folds
ref_key <- ref_keys[1]
folds <- createFolds(ref@meta.data[[ref_key]], k = n_folds, list = TRUE, returnTrain = FALSE)

all_pred_scores <- list()
all_test_meta <- list()
ref_name <- basename(ref_path) %>% gsub(".rds", "", .)

for (i in seq_along(folds)) {
        test_idx <- folds[[i]]
        train_idx <- setdiff(seq_len(ncol(ref)), test_idx)
        train_obj <- subset(ref, cells = colnames(ref)[train_idx])
        test_obj <- subset(ref, cells = colnames(ref)[test_idx])
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
        all_pred_scores[[i]] <- pred_scores %>% as.data.frame()
        all_test_meta[[i]] <- test_obj@meta.data
}

# Combine all folds
prediction_scores <- do.call(rbind, all_pred_scores)
test_meta <- do.call(rbind, all_test_meta)


# Only use ref_name for output in k-fold CV
ref_name = basename(ref_path) %>% gsub(".rds", "", .)

prediction_scores <- prediction_scores %>%
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


write.table(prediction_scores, file=paste0(ref_name,"_kfold_prediction_scores.tsv"), sep="\t", row.names=FALSE)
write.table(test_meta, file=paste0(ref_name,"_kfold.obs.tsv"), row.names=FALSE, sep= "\t")

