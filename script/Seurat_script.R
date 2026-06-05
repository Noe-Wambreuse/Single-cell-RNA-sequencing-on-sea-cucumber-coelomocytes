---
  title: "single cell on sea cucumber coelomocytes"
author: 'Noé'
date: "2025-05-20"
output: html_document
---

#Load the necessary packages

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


# 1. Import data and create seurat object

##Note: because this step require an older version of R, we recommend to upload the R data provided

files <- file.path("data/alevin/")
txi <- tximport(files, type="alevin")

##Creating the Seurat object and filtering gene that are expressed in less than 3 cells and gene with less than 100 counts

obj <- CreateSeuratObject(counts = txi$counts,
                          min.cells = 3, 
                          min.features = 100)


##Display the number of cells and gene in the seurat object
ncol(obj)
nrow(obj) 




# 2. Clustering analysis

## Normalisation 
obj <- NormalizeData(obj)  #(default => scale.factor = 10000)
obj <- FindVariableFeatures(obj) #Identification of the most variable genes => default = 2000)
obj <- ScaleData(obj, features = VariableFeatures(obj), 
                 vars.to.regress = c("nCount_RNA")) #ajusting variables reduce UMIs redundancy

## Visualisation of the ten most varaible genes
top10 <- head(VariableFeatures(obj), 10)
# plot variable features with and without labels
plot_variable <- VariableFeaturePlot(obj)
plot_variable_labels <- LabelPoints(plot = plot_variable, points = top10, repel = TRUE)
plot_variable_labels

## Run PCA
obj <- RunPCA(obj, verbose = T) 

### Examine and visualize PCA results
pca_stdev <- obj[["pca"]]@stdev
pct_var <- round(100 * (pca_stdev^2) / sum(pca_stdev^2), 1)
xlab <- paste0("PC1 (", pct_var[1], "%)")
ylab <- paste0("PC2 (", pct_var[2], "%)")
DimPlot(obj, reduction = "pca", pt.size = 2) +
  labs(x = xlab, y = ylab, color = "Cluster") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )

## Run UMAP
### Building the UMAP and Visualise UMAP
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.25) 
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, min.dist = 0.35)
D1<-DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.8, label.size = 6)+ NoLegend()
D1


## Calculate the number of cell per cluster
cluster_counts <- table(obj$seurat_clusters)
cluster_df <- as.data.frame(cluster_counts)
colnames(cluster_df) <- c("Cluster", "cell_number")

### Calculate their proportion
cluster_df$Proportion <- (cluster_df$cell_number / sum(cluster_df$cell_number)) * 100

### Build the graph
ggplot(cluster_df, aes(x = Cluster, y = cell_number, fill = Cluster)) +
  geom_bar(stat = "identity") +  # PAS de contour
  geom_text(aes(label = paste0(round(Proportion, 1), "%")), vjust = -0.5, size = 6) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)),
    breaks = scales::pretty_breaks(n = 5)
  ) +
  theme_minimal(base_size = 14) +  
  labs(
    x = "Cell cluster",
    y = "Number of cells"
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20),    
    axis.title = element_text(size = 20, face = "bold", color = "black"), 
    panel.grid.major.y = element_line(color = "gray80"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1) 
  )




# 3. Quality control
## Quality control of distribution of feature and RNA

### Extract the data
df_count <- FetchData(obj, vars = c("nCount_RNA", "orig.ident")) %>%
  filter(nCount_RNA > 0) %>%  # avoid log(0)
  mutate(log_nCount = log10(nCount_RNA))

### Build density graphs
ggplot(df_count, aes(x = log_nCount, fill = orig.ident)) +
  geom_density(alpha = 1) +
  geom_vline(xintercept = log10(300), linetype = "dashed", color = "black", size = 0.5) +  
  scale_x_continuous(
    name = "nGene",
    limits = c(2, 4),
    breaks = log10(c(100, 1000, 10000)),
    labels = c("100", "1k", "10k")
  ) +
  scale_fill_manual(values = c("#F8766D")) +
  labs(y = "Cell density") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )

### Extract the data
df_feature <- FetchData(obj, vars = c("nFeature_RNA", "orig.ident")) %>%
  filter(nFeature_RNA > 0) %>%  # avoid log(0)
  mutate(log_nFetaure = log10(nFeature_RNA))

### Build density graphs
ggplot(df_feature, aes(x = log_nFetaure, fill = orig.ident)) +
  geom_density(alpha = 1) +
  geom_vline(xintercept = log10(300), linetype = "dashed", color = "black", size = 0.5) +  
  scale_x_continuous(
    name = "nUMI",
    limits = c(2, 4),
    breaks = log10(c(100, 1000, 10000)),
    labels = c("100", "1k", "10k")
  ) +
  scale_fill_manual(values = c("#F8766D")) +
  labs(y = "Cell density") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  )


## Features versus counts for each cluster
plot_features_scatter <- FeatureScatter(
  obj,
  feature1 = "nFeature_RNA",
  feature2 = "nCount_RNA"
)

plot_features_scatter <- plot_features_scatter +
  labs(x = "nGene",
       y = "nUMI") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    panel.grid = element_blank(),  
  )

plot_features_scatter


## Proportion of mitochondrial genes

mt.genes <- read_excel("data/annotation/mit_gene.xlsx", 
                       col_types = c("text", "text", "numeric", 
                                     "text"))

### Keep only mitochondrial genes with e-value < 1e-20
mt.genes <- mt.genes$gene[mt.genes$e_value < 1e-20]

### Calculate percent mitochondrial reads
obj[["percent.mt"]] <- PercentageFeatureSet(obj, features = mt.genes)

### Visualise mitochondrial gene expression
FeaturePlot(obj, features = "percent.mt", reduction = "umap", cols = c("grey80", "firebrick"), pt.size = 0.8) + labs(title = NULL)

###Calculating average mt gene expression per cluster
mt_by_cluster <- obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(mean_percent_mt = mean(percent.mt, na.rm = TRUE))

### Comparison of average mt gene expression per cluster on a bar plot
library(scales)
Idents(obj) <- "seurat_clusters"
clusters <- levels(obj$seurat_clusters)
cluster_colors <- hue_pal()(length(clusters))
names(cluster_colors) <- clusters

ggplot(mt_by_cluster,
       aes(x = seurat_clusters,
           y = mean_percent_mt,
           fill = seurat_clusters)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  ylab("Mintochondrial genes (%)") +
  xlab("Cluster") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, color = "black"),
    panel.grid.major.y = element_line(color = "gray80"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "black", fill = NA, linewidth = 1
    )
  )

### Comparison of the proportion of cells with a mt gene expression > 20% per cluster on a bar plot
mt20_by_cluster <- obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    prop_mt20 = mean(percent.mt > 20, na.rm = TRUE) * 100
  )
mt20_by_cluster


Idents(obj) <- "seurat_clusters"
clusters <- levels(obj$seurat_clusters)
cluster_colors <- hue_pal()(length(clusters))
names(cluster_colors) <- clusters

ggplot(mt20_by_cluster,
       aes(x = seurat_clusters,
           y = prop_mt20,
           fill = seurat_clusters)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  ylab("> 20% of mitochondrial genes (%)") +
  xlab("Cluster") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, color = "black"),
    panel.grid.major.y = element_line(color = "gray80"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "black", fill = NA, linewidth = 1
    )
  )

### UMAP mapping cells with a mt gene expression > 20%
obj$high_mt20 <- ifelse(obj$percent.mt > 20, "mt > 20%", "mt ≤ 20%")

DimPlot(
  obj,
  reduction = "umap",
  group.by = "high_mt20",
  cols = c("firebrick","grey80"),
  pt.size = 0.8
) +
  labs(title = "Cells with mitochondrial proportion > 20%") + NoLegend()+ ggplot2::theme(plot.title = ggplot2::element_blank())


###Extract some metrics and filtering
sum(obj$percent.mt < 20)
round(100 * sum(obj$percent.mt > 20) / ncol(obj), 2)
mean(obj$percent.mt)
sd(obj$percent.mt)
median(obj$percent.mt)

mt.percent.mean <- obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    prop_mt20 = sd(percent.mt, na.rm = TRUE)
  )
mt.percent.mean

### filter the seurat object by cell expressing more than 20% of mitochondrial genes
obj$old_clusters <- Idents(obj) #define cluster for next visualisations
obj_mt.filtered <- subset(obj, subset = percent.mt < 20)

### UMAP after filtering MT gene expression > 20% and visulation
obj_mt.filtered <- FindNeighbors(obj_mt.filtered, dims = 1:30)
obj_mt.filtered <- FindClusters(obj_mt.filtered, resolution = 0.25) 
obj_mt.filtered <- RunUMAP(obj_mt.filtered, reduction = "pca", dims = 1:30, min.dist = 0.35)
D1_MT_filtered <-DimPlot(obj_mt.filtered, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.8, label.size = 6)+ NoLegend()
D1_MT_filtered

### UMAP with older clusters on the new mt filtered UMAP

DimPlot(obj_mt.filtered, group.by = "old_clusters",label = TRUE, repel = TRUE, pt.size = 0.8, label.size = 6)+ NoLegend()+ ggplot2::theme(plot.title = ggplot2::element_blank())

## Run Doublet finder as a quality control
## Find optimal pK

library(DoubletFinder)
sweep.res <- paramSweep(obj, PCs = 1:20)
sweep.stats <- summarizeSweep(sweep.res)
pK <- find.pK(sweep.stats)$pK[which.max(find.pK(sweep.stats)$BCmetric)]

pK.table <- find.pK(sweep.stats)
pK <- as.numeric(as.character(
  pK.table$pK[which.max(pK.table$BCmetric)]
))

print(pK)
class(pK)

## Estimate expected doublets
nExp <- round(0.1 * ncol(obj))  # 10% typical for 10X
print(nExp)

## Run DoubletFinder
obj_doublet <- doubletFinder(
  obj,
  PCs = 1:20,
  pN = 0.25,
  pK = pK,
  nExp = nExp
)

colnames(obj_doublet@meta.data)

##Visualise doublet on the UMAP
obj_doublet$DF.classifications_0.25_0.3_328 <- factor(
  obj_doublet$DF.classifications_0.25_0.3_328,
  levels = c("Singlet", "Doublet")
)
D1_doublet <- DimPlot(obj_doublet, group.by = "DF.classifications_0.25_0.3_328", cols = c("grey80", "firebrick"), pt.size = 0.8)+ NoLegend()+ ggplot2::theme(plot.title = ggplot2::element_blank())
D1_doublet

##Visualise pANN score acorss the UMAP using a color scale
FeaturePlot(
  obj_doublet,
  features = "pANN_0.25_0.3_328",
  reduction = "umap",
  cols = c("grey90", "firebrick"),
  pt.size = 1.5
) + ggplot2::theme(plot.title = ggplot2::element_blank())

## filter the doublet
single_filtered_obj <- subset(obj_doublet, subset = DF.classifications_0.25_0.3_328 == "Singlet")

## Running Umap after filtering doublets
single_filtered_obj <- FindNeighbors(single_filtered_obj, dims = 1:30)
single_filtered_obj <- FindClusters(single_filtered_obj, resolution = 0.25) 
single_filtered_obj <- RunUMAP(single_filtered_obj, reduction = "pca", dims = 1:30, min.dist = 0.35)
D1_doublet.filtered <-DimPlot(single_filtered_obj, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.8, label.size = 6)+ NoLegend()
D1_doublet.filtered

###Map older cluster on the new UMAP
DimPlot(single_filtered_obj, group.by = "old_clusters",label = TRUE, repel = TRUE, pt.size = 0.8, label.size = 6)+ NoLegend()+ ggplot2::theme(plot.title = ggplot2::element_blank())

###Identify the proportion of doublet and the paNN score per cluster
###proportion of doublet
df_doublet <- obj_doublet@meta.data

proportion_doublets <- df_doublet %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_cells = n(),
    n_doublets = sum(DF.classifications_0.25_0.3_328 == "Doublet"),
    prop_doublets = 100 * n_doublets / n_cells
  )
proportion_doublets

###pANN score
pANN_by_cluster <- df_doublet  %>%
  group_by(seurat_clusters) %>%
  summarise(mean_pANN = mean(pANN_0.25_0.3_328))
pANN_by_cluster

###Plot the results
Idents(obj) <- "seurat_clusters"
clusters <- levels(obj$seurat_clusters)
cluster_colors <- hue_pal()(length(clusters))
names(cluster_colors) <- clusters

ggplot(proportion_doublets, aes(x = seurat_clusters, y = prop_doublets, fill = seurat_clusters)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  ylab("Doublet proportion (%)") +
  xlab("Cluster") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),    
        axis.title = element_text(size = 16, color = "black"), 
        panel.grid.major.y = element_line(color = "gray80"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

#Average paNN score
ggplot(pANN_by_cluster, aes(x = seurat_clusters, y = mean_pANN, fill = seurat_clusters)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  ylab("Mean pANN score") +
  xlab("Cluster") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),    
        axis.title = element_text(size = 16, color = "black"), 
        panel.grid.major.y = element_line(color = "gray80"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))




# 3. Identification of markers genes and visualisation
Markers <- FindAllMarkers(obj, return.thresh = 0.0001, only.pos = T) 
write.table(Markers, file = "Markers.csv", sep = "\t", dec = ".", col.names = T, row.names = F)

#Filetring by avg_log2FC and Padj
Marker <- Markers %>%
  filter(avg_log2FC >= 0.5, p_val_adj <= 0.01)

#Creation of gene marker graph
Dotplot <-function(x,y,z){
  Data.to.plot <- FetchData(x, vars = y)
  Data.to.plot$cell <- rownames(Data.to.plot)
  Data.to.plot$cluster <- x@meta.data$seurat_clusters
  
  Data.to.plot <- Data.to.plot %>% gather(key = genes.plot, 
                                          value = expression, -c(cell, cluster))
  Data.to.plot <- Data.to.plot %>% group_by(cluster, genes.plot) %>% 
    summarize(avg.exp = ExpMean(x = expression), pct.exp = PercentAbove(x = expression, threshold = 0), Expr = "0")
  ZSCORE<-Data.to.plot %>% group_by(genes.plot) %>% 
    summarize(M = ExpMean(x=avg.exp), SD = ExpSD(x=avg.exp))
  Data.to.plot$zscore = (Data.to.plot$avg.exp - ZSCORE$M) / ZSCORE$SD
  test <- Data.to.plot %>% group_by(genes.plot) %>% summarize(gene.expr = mean(x = avg.exp), gene.max = max(x=avg.exp))
  Data.to.plot$Ratio = Data.to.plot$avg.exp / test$gene.max
  
  G1<-ggplot(Data.to.plot, mapping = aes(x = cluster, y = genes.plot)) +
    geom_point(mapping = aes(size = pct.exp, color = zscore)) + theme(axis.text.x = element_text(angle = 0, hjust = 50)) + 
    scale_color_gradientn(colours = c("Blue", "Grey", "Yellow", "Orange", "Red"), guide = "colourbar", values = c(0,0.25,0.5,0.75,1)) + 
    scale_y_discrete(limits=y) + ggtitle(z) + theme()
  
  G2<-ggplot(Data.to.plot, mapping = aes(x = cluster, y = genes.plot)) +
    geom_point(mapping = aes(size = pct.exp, color = Ratio)) + theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 16)) + 
    scale_color_gradientn(colours = c("Blue", "Cyan", "Yellow", "Orange", "Red"), guide = "colourbar") + 
    scale_y_discrete(limits=y) + ggtitle(z) + theme()
  G3 <-plot_grid(G1,G2,ncol=2)
  return(G2)
}


### Adding annotation information of marker genes
annotation <- read_excel("~/R shared/single_cell_new/data/alevin/annotation.xlsx")
Markers_and_annotation <- left_join(Marker, annotation, by = "gene")

###Filtering for only having marker genes annotated
Markers_filtered_annotation <- Markers_and_annotation %>% filter(!is.na(annotation), annotation != "NA")

###Isoltate TOP5 marker gene for each cluster
Top5_filtered <- Markers_filtered_annotation %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
Top5 <- Markers_and_annotation %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

### Build dotplots
P1 <-Dotplot(obj, unique(Top5$gene), "De novo markers")
P1_filtered <-Dotplot(obj, unique(Top5_filtered$gene), "De novo markers")

### Visualisation
P1
P1_filtered

### Visualising the expression of one gene of interest
F1<-FeaturePlot(obj, features = "CL2228.Contig1", order=T, pt.size = 1, cols = c("grey", "blue"))
F1

### Create a signature based on marker genes
signature_list <- Top5 %>% 
  group_by(cluster) %>%
  summarise(genes = list(unique(gene))) %>%
  deframe() 

signature_list_NA_filtered <- Top5_filtered %>% 
  group_by(cluster) %>%
  summarise(genes = list(unique(gene))) %>%
  deframe() 

### Calculate the percentage score
for (clust in names(signature_list_NA_filtered)) {
  gene_set <- signature_list[[clust]]
  obj[[paste0("signature_", clust)]] <- PercentageFeatureSet(obj, features = gene_set)
}

for (clust in names(signature_list)) {
  gene_set <- signature_list[[clust]]
  obj[[paste0("signature_", clust)]] <- PercentageFeatureSet(obj, features = gene_set)
}

### Visualise signature on the umap
FeaturePlot(
  obj,
  features = paste0("signature_", names(signature_list)),
  ncol = 3,
  pt.size = 0.3,      
  order = TRUE,      
  min.cutoff = "q10",  #Improve contrast
  max.cutoff = "q90"   #Avoid that the scale is dominated by outliners
)

FeaturePlot(
  obj,
  features = paste0("signature_", names(signature_list_NA_filtered)),
  ncol = 3,
  pt.size = 0.3,      
  order = TRUE,      
  min.cutoff = "q10",  #Improve contrast
  max.cutoff = "q90"   #Avoid that the scale is dominated by outliners
)




#8. Get average expression per cluster

###Isolate avg expression
avg_exp <- AverageExpression(obj, assays = "RNA", slot = "data", return.seurat = FALSE)
avg_exp_matrix <- avg_exp$RNA  # rows = genes, columns = clusters

clusters <- Idents(obj)
genes <- rownames(obj)
    
pct_mat <- sapply(levels(clusters), function(clust) {
      cells <- WhichCells(obj, idents = clust)
      expr <- GetAssayData(obj, slot = "data")[, cells]
      rowSums(expr > 0) / length(cells) * 100
    })
    
rownames(pct_mat) <- genes
    
### Convert sparse matrices to regular data frames
avg_exp_df <- as.data.frame(as.matrix(avg_exp_matrix))
pct_df <- as.data.frame(as.matrix(pct_mat))

### Add suffixes to column names
colnames(avg_exp_df) <- paste0("avg_", colnames(avg_exp_df))
colnames(pct_df) <- paste0("pct_", colnames(pct_df))

### Combine
combined_wide <- cbind(avg_exp_df, pct_df)

### Add gene names as a column
combined_wide <- cbind(gene = rownames(combined_wide), combined_wide)

### Save
write.csv(combined_wide, "data/gene_avg_exp.csv", row.names = FALSE)


#. General informcation
sessionInfo()
