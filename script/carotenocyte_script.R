---
  title: "Identification of carotenocyte cluster"
author: 'Noé Wambreuse'
date: "2025-05-20"
output: html_document
---
  
#Import data

##For this analysis, we advise loading directly the Rdata "Seurat_object_UMAP"

load("data/Seurat_object_UMAP.RDATA")

#Load necessary packages
library(patchwork)
library(readxl)
library(Seurat)
library(tibble)
library(tidyr)    
library(ggplot2)  
library(cowplot)   
library(dplyr)
library(scales)
library(tximport)
library(tximeta)
library(DoubletFinder)
library(ggplot2)
library(SingleCellExperiment)
library(scran)
library(SingleR)

# 1. Signature of carotenocytes

library(readxl)
gene_list <- read_excel("data/annotation/gene_carotenocyte.xlsx", col_names = TRUE)
genes <- gene_list$gene 

obj$signature_score <- Matrix::colMeans(GetAssayData(obj, slot = "data")[genes, , drop = FALSE])

## UMAP with the signature
FeaturePlot(obj, features = "signature_score", pt.size = 0.3) +
  scale_color_gradient(low = "grey90", high = "#FF5353") +
  labs(title = "signature expression", color = "Score (%)") +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

## Violin plot of the signature by cluster
VlnPlot(obj, features = "signature_score", group.by = "seurat_clusters", pt.size = 0.5) +
  scale_fill_manual(values = c("grey90", "grey90", "grey90", "grey90", "grey90", "grey90", "#FF5353", "grey90", "grey90", "grey90")) +
  labs(title = "signature expression", y = "Score (%)", x = "Clusters") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

annotate_with_SingleR <- function(obj, ref_file, ref_labels_vec = c("PF", "HF"), label_col = "label.main")
  
  
#2. Use SingleR to identify carotenocyte cluster
  
## Prepare the single cell object
sce <- as.SingleCellExperiment(obj)
sce <- logNormCounts(sce)

## load the differential RNA-seq matrix
rna_seq_deg <- read_excel("data/annotation/bRNAseq_diff_exp.xlsx")

## Ajust the matrix
rownames_mat <- rna_seq_deg[[1]]
counts_mat <- as.matrix(rna_seq_deg[, -1])
rownames(counts_mat) <- rownames_mat

## Create the SingleCellExperiment object
sce_ref <- SingleCellExperiment(assays = list(counts = counts_mat))
ref_labels_vec <- c("PF", "HF")
ref_labels <- rep(ref_labels_vec, length.out = ncol(counts_mat))
colData(sce_ref)[["label_main"]] <- factor(ref_labels)
sce_ref <- logNormCounts(sce_ref)

## Annotation with SingleR
pred <- SingleR(test = sce, ref = sce_ref, labels = sce_ref[["label_main"]])
sce$SingleR_labels <- pred$labels

## Convert to Seurat object
seurat_annot <- as.Seurat(sce)

## Visualisation
custom_colors <- c(
  "PF" = "#D9D9D9",
  "HF" = "#FF5353"
)

p1 <- DimPlot(seurat_annot, group.by = "SingleR_labels", label = FALSE, pt.size = 0.5) +
  scale_color_manual(values = custom_colors) +  
  ggtitle("Annotation SingleR") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),       
    panel.background = element_blank(),      
  )

p1

## Calculate the proportion by cluster
cell_types <- sce$SingleR_labels
clusters <- as.factor(sce$seurat_clusters)
df <- data.frame(cluster = clusters, cell_type = cell_types)

proportion_df <- df %>%
  group_by(cluster, cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(proportion = count / sum(count))

## Visualisation

proportion_df$cell_type <- factor(proportion_df$cell_type, levels = c("HF", "PF"))
custom_colors <- c(
  "HF" = "#FF5353", 
  "PF" = "#D9D9D9"   
)

p2 <- ggplot(proportion_df, aes(x = cluster, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +  # reverse l'empilement
  scale_fill_manual(values = custom_colors) +
  labs(
    x = "Cluster",
    y = "Cell proportion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )

p2

#General info
sessionInfo()
