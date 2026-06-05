---
  title: "Functional enrichement based on marker genes"
author: 'Noé Wambreuse'
date: "2025-05-20"
output: html_document
---
  
#Load necessary packages

library(tibble)
library(tidyverse)
library(GO.db)
library(clusterProfiler)
library(enrichplot) 
library(openxlsx)
library(readxl)
library(stringr)
library(dplyr)

# 1. KEGG functional enrichment
#load the KEGG annotation file

library(readxl)
ko_genes <- read_excel("data/annotation/ko_genes.xlsx")

#Filter the unannotated genes
ko_genes_filtered<- ko_genes %>%
  filter(!is.na(KO))

# Create a list of KEGG terms by cluster
ko_lists_by_cluster <- ko_genes_filtered %>%
  group_by(cluster) %>%
  summarise(KO = list(unique(KO))) %>%
  deframe()

str(ko_lists_by_cluster)

## Perform the enrichiment analysis
kegg_enrich <- compareCluster(
  gene = ko_lists_by_cluster,
  fun = "enrichKEGG",
  organism = "ko",         
  keyType = "kegg",        
  pAdjustMethod = "BH",
  pvalueCutoff = 1,  # No filtering
  qvalueCutoff = 1   # No filtering
)

## Create the functional enrichment table
kegg_table <- as.data.frame(kegg_enrich)

## Export the table with the unigene list per KO
file_path <- "data/table"
file_path <- file.path(file_path, "kegg_table.xlsx")

ko_to_gene_cluster <- ko_genes_filtered %>%
  select(cluster, KO, gene) %>%
  distinct()

map_ko_to_genes_cluster <- function(ko_string, cluster_name) {
  ko_list <- str_split(ko_string, "/")[[1]]
  gene_list <- ko_to_gene_cluster %>%
    filter(cluster == cluster_name, KO %in% ko_list) %>%
    pull(gene) %>%
    unique()
  paste(gene_list, collapse = "/")
}

kegg_table$geneNames <- mapply(
  map_ko_to_genes_cluster,
  kegg_table$geneID,
  kegg_table$Cluster
)

write.xlsx(kegg_table, file = file_path)



## Built the graph using ggplot
### Keep the five best terms per cluster based on the Fold enrichment value
top5_kegg <- kegg_table %>%
  group_by(Cluster) %>%
  slice_max(FoldEnrichment, n = 5) %>%
  ungroup()

top5_kegg$Description <- str_trunc(top5_kegg$Description, 50)

ggplot(top5_kegg, aes(x = factor(Cluster), 
                      y = reorder(Description, FoldEnrichment), 
                      size = Count, 
                      color = -log10(qvalue))) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(3, 8)) + 
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal(base_size = 12) +
  labs(
    title = "",
    x = "Cluster",
    y = "KEGG pathways",
    color = expression(-log[10](qvalue)),
    size = "Gene number"
  ) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.y = element_text(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.title = element_blank(),
    panel.border = element_rect(color = "black", fill = "NA", size = 0.5)  
  )



# 2. Biological process (BP) functional enrichment
## load the GO annotation file
go_genes <- read_excel("data/annotation/go_genes.xlsx")

# Filter the unannotated genes
df_clean <- go_genes %>%
  filter(!is.na(GO) & GO != "NA" & GO != "")

# Keep only GO terms corresponding to "biological_process"
extract_bp_go <- function(go_string) {
  go_ids <- unique(unlist(regmatches(go_string, gregexpr("GO:[0-9]+", go_string))))
  types <- Ontology(go_ids)
  bp_ids <- go_ids[types == "BP"]
  return(bp_ids)
}

# Build the BP terms matrix
df_long <- df_clean %>%
  rowwise() %>%
  mutate(GO_BP = list(extract_bp_go(GO))) %>%
  tidyr::unnest(GO_BP) %>%
  dplyr::select(gene, cluster, GO_BP)

# Add the description of each BP
df_long <- df_long %>%
  mutate(Description = Term(GO_BP)) %>%
  filter(!is.na(Description) & Description != "")

# Prepare data for the functional enrichment analysis
genes_by_cluster <- split(df_long$gene, df_long$cluster)

term2gene <- df_long %>%
  select(GO_BP, gene) %>%
  distinct()

term2name <- df_long %>%
  select(GO_BP, Description) %>%
  distinct()

# BP functional enrichment anaysis
enrich_results <- lapply(names(genes_by_cluster), function(clust) {
  enricher(
    gene = unique(genes_by_cluster[[clust]]),
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pvalueCutoff = 1,   
    qvalueCutoff = 1
  )
})
names(enrich_results) <- names(genes_by_cluster)

# Combine result in one data frame
enrich_BP_table <- do.call(rbind, lapply(names(enrich_results), function(clust) {
  res <- enrich_results[[clust]]
  
  if (!is.null(res) && inherits(res, "enrichResult") && nrow(res@result) > 0) {
    res_df <- res@result
    res_df$cluster <- clust
    return(res_df)
  } else {
    return(NULL)
  }
}))

# Exporte the data
file_path <- "data/table"
file_path <- file.path(file_path, "enrich_BP_table.xlsx")
write.xlsx(enrich_BP_table, file = file_path)


## Built the graph using ggplot
### Keep the five best terms per cluster based on the Fold enrichment value
top5_enrich <- enrich_BP_table %>%
  group_by(cluster) %>%
  slice_max(FoldEnrichment, n = 5) %>%
  ungroup()
top5_enrich$Description <- str_trunc(top5_enrich$Description, 50)

ggplot(top5_enrich, aes(x = factor(cluster), 
                        y = reorder(Description, FoldEnrichment), 
                        size = Count, 
                        color = -log10(qvalue))) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(3, 8)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal(base_size = 12) +
  labs(
    title = "",
    x = "Cluster",
    y = "GO terms (BP)",
    color = expression(-log[10](qvalue)),
    size = "gene number"
  ) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.y = element_text(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )



# 3. Cellular (CC) functional enrichment
## Keep only GO terms corresponding to "Cellular component"

extract_cc_go <- function(go_string) {
  go_ids <- unique(unlist(regmatches(go_string, gregexpr("GO:[0-9]+", go_string))))
  types <- Ontology(go_ids)
  cc_ids <- go_ids[types == "CC"]
  return(cc_ids)
}

# Build the CC terms matrix
df_cc_long <- df_clean %>%
  rowwise() %>%
  mutate(GO_CC = list(extract_cc_go(GO))) %>%
  tidyr::unnest(GO_CC) %>%
  dplyr::select(gene, cluster, GO_CC)

# Add the description of each CC
df_cc_long <- df_cc_long %>%
  mutate(Description = Term(GO_CC)) %>%
  filter(!is.na(Description) & Description != "")

# Prepare data for the functional enrichment analysis
genes_by_cluster <- split(df_cc_long$gene, df_cc_long$cluster)

term2gene <- df_cc_long %>%
  dplyr::select(GO_CC, gene) %>%
  distinct()

term2name <- df_cc_long %>%
  dplyr::select(GO_CC, Description) %>%
  distinct()

# CC functional enrichment analysis
enrich_results_cc <- lapply(names(genes_by_cluster), function(clust) {
  enricher(
    gene = unique(genes_by_cluster[[clust]]),
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pvalueCutoff = 1,   # pas de filtre pour explorer tout
    qvalueCutoff = 1
  )
})
names(enrich_results_cc) <- names(genes_by_cluster)

# Combine the results in one data.frame
enrich_CC_table <- do.call(rbind, lapply(names(enrich_results_cc), function(clust) {
  res <- enrich_results_cc[[clust]]
  
  if (!is.null(res) && inherits(res, "enrichResult") && nrow(res@result) > 0) {
    res_df <- res@result
    res_df$cluster <- clust  
    return(res_df)
  } else {
    return(NULL)
  }
}))


# Export the table
file_path <- "data/table"
file_path <- file.path(file_path, "enrich_CC_table.xlsx")
write.xlsx(enrich_CC_table, file = file_path)


## Built the graph using ggplot
### Keep the five best terms per cluster based on the Fold enrichment value
top5_enrich_cc <- enrich_CC_table %>%
  group_by(cluster) %>%
  slice_max(FoldEnrichment, n = 5) %>%
  ungroup()
top5_enrich_cc$Description <- str_trunc(top5_enrich_cc$Description, 50)

ggplot(top5_enrich_cc, aes(x = factor(cluster), 
                           y = reorder(Description, FoldEnrichment), 
                           size = Count, 
                           color = -log10(qvalue))) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(3, 8)) +  
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal(base_size = 12) +
  labs(
    title = "",
    x = "Cluster",
    y = "GO terms (CC)",
    color = expression(-log[10](qvalue)),
    size = "gene number"
  ) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.y = element_text(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )


# 4. Molecular function (MF) functional enrichment
## Keep only GO terms corresponding to "Molecular function"

extract_mf_go <- function(go_string) {
  go_ids <- unique(unlist(regmatches(go_string, gregexpr("GO:[0-9]+", go_string))))
  types <- Ontology(go_ids)
  mf_ids <- go_ids[types == "MF"]
  return(mf_ids)
}

# Build the BP terms matrix
df_mf_long <- df_clean %>%
  rowwise() %>%
  mutate(GO_MF = list(extract_mf_go(GO))) %>%
  tidyr::unnest(GO_MF) %>%
  dplyr::select(gene, cluster, GO_MF)

# Add the MF description
df_mf_long <- df_mf_long %>%
  mutate(Description = Term(GO_MF)) %>%
  filter(!is.na(Description) & Description != "")


# Prepare the object for the functional enrichment
genes_by_cluster <- split(df_mf_long$gene, df_mf_long$cluster)

library(dplyr)
term2gene <- df_mf_long %>% dplyr::select(GO_MF, gene) %>% distinct()
term2name <- df_mf_long %>% dplyr::select(GO_MF, Description) %>% distinct()

# MF functional enrichment
enrich_results_mf <- lapply(names(genes_by_cluster), function(clust) {
  enricher(
    gene = unique(genes_by_cluster[[clust]]),
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
})

names(enrich_results_mf) <- names(genes_by_cluster)

# Combine the results in one data.frame
enrich_MF_table <- do.call(rbind, lapply(names(enrich_results_mf), function(clust) {
  res <- enrich_results_mf[[clust]]
  
  if (!is.null(res) && inherits(res, "enrichResult") && nrow(res@result) > 0) {
    res_df <- res@result
    res_df$cluster <- clust  # Ajouter le numéro ou nom du cluster
    return(res_df)
  } else {
    return(NULL)
  }
}))

#Export the table
file_path <- "data/table"
file_path <- file.path(file_path, "enrich_MF_table.xlsx")
write.xlsx(enrich_MF_table, file = file_path)

## Built the graph using ggplot
### Keep the five best terms per cluster based on the Fold enrichment value
top5_enrich_mf <- enrich_MF_table %>%
  group_by(cluster) %>%
  slice_max(FoldEnrichment, n = 5) %>%
  ungroup()

top5_enrich_mf$Description <- str_trunc(top5_enrich_mf$Description, 50)

ggplot(top5_enrich_mf, aes(x = factor(cluster), 
                           y = reorder(Description, FoldEnrichment), 
                           size = Count, 
                           color = -log10(qvalue))) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(3, 8)) +  # <- size of the point (min 3, max 8)
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal(base_size = 12) +
  labs(
    title = "",
    x = "Cluster",
    y = "GO terms (MF)",
    color = expression(-log[10](qvalue)),
    size = "gene number"
  ) +
  theme(
    text = element_text(family = "Arial"), 
    axis.text.y = element_text(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

#General information
sessionInfo()
