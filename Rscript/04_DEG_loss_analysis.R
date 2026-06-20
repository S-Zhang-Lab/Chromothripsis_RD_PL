##========================================##
##          DEG analysis                  ##
##========================================##
#Load Packages
library(dplyr)
library(Seurat)
library(ggplot2)
library(sctransform)
library(scales)
library(plotly)

# Set seed
set.seed(42)

#set working directory for loading in raw data
#Note: opened 03_3p_loss.RData for following analysis
setwd("/path/to/Chromothripsis_2")  # <-- set to your local path to the project root
head(obj@meta.data)
Idents(obj) <- "seurat_clusters"
DimPlot(obj, label = TRUE)

### Show 3p loss cells on UMAP/FeaturePlot
library(Seurat)
library(ggplot2)

# Create a new column in metadata for binary labeling (3p loss vs other)
obj$highlight_3p_loss <- ifelse(obj$recoded_loss_chr3_2_bin == "02_chr3_loss", "3p_loss", "other")

# Set colors: gray for "other", and a bright color (e.g., purple) for "3p_loss"
color_mapping <- c("other" = "gray", "3p_loss" = "purple")

# Plot the UMAP using DimPlot to color cells based on the new metadata column
DimPlot(
  object = obj,
  group.by = "highlight_3p_loss",  # Use the new binary column
  reduction = "umap",  # Assuming you're using UMAP for dimensional reduction
  cols = color_mapping  # Apply the custom color mapping
) +
  theme_minimal() +  # Clean up the background
  labs(title = "UMAP Highlighting 3p Loss Cells") + 
  theme(legend.position = "none")  # Remove the legend if not needed


### DEG from cells with 3p loss on Day5 vs Day10 
### goal: to assess if there are any differential cellular responses to 3p loss over time
# Lane 3 cells are cells from day 10 and lane 1,2 cells are from day 5

# Assuming your Seurat object is called 'seurat_object'

# Create a new column 'timepoint' based on HTO_classification
obj@meta.data$timepoint <- ifelse(grepl("^L1|^L2", obj@meta.data$HTO_classification), 
                                            "day5", 
                                            ifelse(grepl("^L3", obj@meta.data$HTO_classification), 
                                                   "day10", 
                                                   NA))  # Assign NA if it doesn't match any pattern

# Check if the 'timepoint' column was added correctly
table(obj@meta.data$timepoint)
# day10  day5 
# 11669 13451

# Subset day5 cells
day5_cells <- subset(obj, subset = timepoint == "day5")

# Subset day10 cells
day10_cells <- subset(obj, subset = timepoint == "day10")

# DEG analysis for day5 cells: Chromosome 3 loss vs No_loss
Idents(day5_cells) <- "recoded_loss_chr3_2_bin"
day5_cells <- PrepSCTFindMarkers(day5_cells)
day5_DEGs <- FindMarkers(day5_cells, ident.1 = "02_chr3_loss", ident.2 = "01_No_loss", verbose = TRUE, logfc.threshold = 0.05)
# View top results
head(day5_DEGs)

# DEG analysis for day10 cells: Chromosome 3 loss vs No_loss
Idents(day10_cells) <- "recoded_loss_chr3_2_bin"
day10_cells <- PrepSCTFindMarkers(day10_cells)
day10_DEGs <- FindMarkers(day10_cells, ident.1 = "02_chr3_loss", ident.2 = "01_No_loss", verbose = TRUE, logfc.threshold = 0.05)

# View top results
head(day10_DEGs)

write.csv(day10_DEGs, file = "3pLoss_vs_NoLoss_day10_DEGs.csv", row.names = TRUE)
write.csv(day5_DEGs, file = "3pLoss_vs_NoLoss_day5_DEGs.csv", row.names = TRUE)

# For day5 DEGs
day5_DEGs$log2FC <- day5_DEGs$avg_log2FC
day5_DEGs$neg_log10_pval <- -log10(day5_DEGs$p_val_adj)

# For day10 DEGs
day10_DEGs$log2FC <- day10_DEGs$avg_log2FC
day10_DEGs$neg_log10_pval <- -log10(day10_DEGs$p_val_adj)

library(ggplot2)

# Volcano plot for day5 DEGs
ggplot(day5_DEGs, aes(x = log2FC, y = neg_log10_pval)) +
  geom_point(aes(color = (p_val_adj < 0.05 & abs(log2FC) > 0.25)), size = 2) +
  scale_color_manual(values = c("grey", "blue")) +  # Customize colors
  theme_minimal() +                                 # Use a minimal theme
  labs(title = "DEGs: Chromosome 3 Loss vs No_loss (Day 5)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.position = "none") +                 # Remove legend
  geom_text_repel(data = subset(day5_DEGs, p_val_adj < 0.05 & abs(log2FC) > 0.5), # Add labels for significant DEGs
                  aes(label = rownames(subset(day5_DEGs, p_val_adj < 0.05 & abs(log2FC) > 0.5))),
                  size = 3, max.overlaps = 10)      # Control overlap and label size

# Volcano plot for day10 DEGs
ggplot(day10_DEGs, aes(x = log2FC, y = neg_log10_pval)) +
  geom_point(aes(color = (p_val_adj < 0.05 & abs(log2FC) > 0.25)), size = 2) +
  scale_color_manual(values = c("grey", "red")) +  # Customize colors
  theme_minimal() +                                 # Use a minimal theme
  labs(title = "DEGs: Chromosome 3 Loss vs No_loss (Day 10)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.position = "none") +                 # Remove legend
  geom_text_repel(data = subset(day10_DEGs, p_val_adj < 0.05 & abs(log2FC) > 0.5), # Add labels for significant DEGs
                  aes(label = rownames(subset(day10_DEGs, p_val_adj < 0.05 & abs(log2FC) > 0.5))),
                  size = 3, max.overlaps = 10)      # Control overlap and label size


### Clustering information for 3p loss cells - split by timepoint 
# Set seed
set.seed(42)
#re-cluster subsetted cell type
day5_cells <- RunPCA(day5_cells, verbose = FALSE)
Idents(day5_cells) <- "seurat_clusters"
day5_cells <- RunUMAP(day5_cells, dims = 1:10, verbose = FALSE)
day5_cells <- FindNeighbors(day5_cells, dims = 1:30, verbose = FALSE)
day5_cells <- FindClusters(day5_cells, verbose = FALSE, resolution = 0.2)
DimPlot(day5_cells, label = TRUE)
DimPlot(day5_cells, split.by = "HTO_classification")

Idents(day5_cells) <- "seurat_clusters"
day5.markers <- FindAllMarkers(day5_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
day5.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
View(day5.markers)
top10_day5_DEG <- day5.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(day5_cells, features = top10_day5_DEG$gene) +scale_fill_gradientn(colors = c("blue", "white", "red")) 
write.csv(day5.markers,"day5_marker_genes_subsetted.csv")

#Compare Cluster Proportions Between day5 loss and no loss cells
table(day5_cells@meta.data$seurat_clusters,day5_cells@meta.data$recoded_loss_chr3_2_bin)
freq_table <- prop.table(x = table(day5_cells@meta.data$seurat_clusters,day5_cells@meta.data$recoded_loss_chr3_2_bin),margin = 2)
barplot(height = freq_table) 
coloridentities <- levels(day5_cells@meta.data$seurat_clusters) 
my_color_palette <- hue_pal()(length(coloridentities))
barplot(height = freq_table, col = my_color_palette)
# seems cluster 5 is expanded in day5 chr3 loss cells, find markers 

C5.markers <- FindMarkers(day5_cells, ident.1 ="5", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(C5.markers, file = "DEG_C5_compared_all.csv")

# Set seed
set.seed(42)
#re-cluster subsetted cell type
day10_cells <- RunPCA(day10_cells, verbose = FALSE)
Idents(day10_cells) <- "seurat_clusters"
day10_cells <- RunUMAP(day10_cells, dims = 1:10, verbose = FALSE)
day10_cells <- FindNeighbors(day10_cells, dims = 1:30, verbose = FALSE)
day10_cells <- FindClusters(day10_cells, verbose = FALSE, resolution = 0.2)
DimPlot(day10_cells, label = TRUE)
DimPlot(day10_cells, split.by = "HTO_classification")

Idents(day10_cells) <- "seurat_clusters"
day10.markers <- FindAllMarkers(day10_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
day10.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
View(day10.markers)
top10_day10_DEG <- day10.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(day10_cells, features = top10_day10_DEG$gene) +scale_fill_gradientn(colors = c("blue", "white", "red")) 
write.csv(day10.markers,"day10_marker_genes_subsetted.csv")

### Plot 3p loss cells on UMAP for each time point
# Set colors: gray for "other", and a bright color (e.g., purple) for "3p_loss"
color_mapping <- c("other" = "gray", "3p_loss" = "purple")

# Plot the UMAP using DimPlot to color cells based on the new metadata column
DimPlot(
  object = day5_cells,
  group.by = "highlight_3p_loss",  # Use the new binary column
  reduction = "umap",  # Assuming you're using UMAP for dimensional reduction
  cols = color_mapping  # Apply the custom color mapping
) +
  theme_minimal() +  # Clean up the background
  labs(title = "UMAP Highlighting 3p Loss Cells (day5)") + 
  theme(legend.position = "none")  # Remove the legend if not needed

DimPlot(
  object = day10_cells,
  group.by = "highlight_3p_loss",  # Use the new binary column
  reduction = "umap",  # Assuming you're using UMAP for dimensional reduction
  cols = color_mapping  # Apply the custom color mapping
) +
  theme_minimal() +  # Clean up the background
  labs(title = "UMAP Highlighting 3p Loss Cells (day10)") + 
  theme(legend.position = "none")  # Remove the legend if not neede

#Compare Cluster Proportions Between day10 loss and no loss cells
table(day10_cells@meta.data$seurat_clusters,day10_cells@meta.data$recoded_loss_chr3_2_bin)
freq_table <- prop.table(x = table(day10_cells@meta.data$seurat_clusters,day10_cells@meta.data$recoded_loss_chr3_2_bin),margin = 2)
barplot(height = freq_table) 
coloridentities <- levels(day10_cells@meta.data$seurat_clusters) 
my_color_palette <- hue_pal()(length(coloridentities))
barplot(height = freq_table, col = my_color_palette)
# seems cluster 1 is expanded in day10 chr3 loss cells, find markers 

C1.markers <- FindMarkers(day10_cells, ident.1 ="1", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(C1.markers, file = "DEG_C1_compared_all_day10.csv")

saveRDS(day5_cells, file = "day5_cells_EA.rds")
saveRDS(day10_cells, file = "day10_cells_EA.rds")

###DEGs from cells with yq loss and compare it with genes that have 3p loss to confirm
###3p specific cellular response (can use Venn diagram)
### First: I will generate the venn diagram for the no loss comparison

# Ensure the plotting area is cleared
dev.off()

# Set the significance threshold for p-value and log2 fold change
pval_threshold <- 0.01
logfc_threshold <- 0.5

# Filter upregulated genes for day 5 (log2FC > 0.25 and p_val_adj < 0.05)
day5_upregulated <- rownames(subset(day5_DEGs, avg_log2FC > logfc_threshold & p_val_adj < pval_threshold))

# Filter upregulated genes for day 10 (log2FC > 0.25 and p_val_adj < 0.05)
day10_upregulated <- rownames(subset(day10_DEGs, avg_log2FC > logfc_threshold & p_val_adj < pval_threshold))

# Use the VennDiagram package to create a Venn diagram
library(VennDiagram)

# Create a list of the upregulated genes for comparison
upregulated_list <- list("Day 5 Upregulated" = day5_upregulated, "Day 10 Upregulated" = day10_upregulated)

# Create the Venn diagram
venn.plot <- venn.diagram(
  x = upregulated_list,
  category.names = c("Day 5 Upregulated", "Day 10 Upregulated"),
  fill = c("blue", "red"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  filename = NULL
)

# Plot the Venn diagram
grid.draw(venn.plot)

# Step 4: Find unique and intersecting genes
# Unique to day 5
unique_day5 <- setdiff(day5_upregulated, day10_upregulated)

# Unique to day 10
unique_day10 <- setdiff(day10_upregulated, day5_upregulated)

# Shared between day 5 and day 10
shared_genes <- intersect(day5_upregulated, day10_upregulated)

# View the results
unique_day5
unique_day10
shared_genes

# If you want to save the gene lists as text files:
write.table(unique_day5, "unique_day5_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(unique_day10, "unique_day10_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(shared_genes, "shared_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

#-------------------------------------------------------------------------------
### Define cells with y chromosome loss

day5_cells <- readRDS("day5_cells_EA.rds")
day10_cells <- readRDS("day10_cells_EA.rds")
DefaultAssay(day5_cells) <- "SCT"
DefaultAssay(day10_cells) <- "SCT"

# Define cells with Y chromosome loss
# These are the genes of interest to identify Y chromosome loss
y_loss_genes <- c("EIF1AY", "KDM5D", "RPS4Y2", "TMSB4Y")
rps4y1_gene <- "RPS4Y1"

# Define loss threshold: cells with expression below this value for Y genes are considered "lost"
loss_threshold <- 0.1  # Adjust as needed

# Fetch the expression levels of the genes
expr_y_loss <- FetchData(day5_cells, vars = y_loss_genes)
expr_rps4y1 <- FetchData(day5_cells, vars = rps4y1_gene)

expr_y_loss.10 <- FetchData(day10_cells, vars = y_loss_genes)
expr_rps4y1.10 <- FetchData(day10_cells, vars = rps4y1_gene)

# Identify cells with loss of Y genes (below threshold for all Y loss genes and above threshold for RPS4Y1)
cells_y_loss <- rownames(expr_y_loss)[
  rowSums(expr_y_loss < loss_threshold) == length(y_loss_genes) &  # Loss of Y genes
    expr_rps4y1 > loss_threshold  # Retained RPS4Y1 expression
]

cells_y_loss.10 <- rownames(expr_y_loss.10)[
  rowSums(expr_y_loss.10 < loss_threshold) == length(y_loss_genes) &  # Loss of Y genes
    expr_rps4y1.10 > loss_threshold  # Retained RPS4Y1 expression
]


# Add Y chromosome loss to metadata
day5_cells$y_loss <- ifelse(rownames(day5_cells@meta.data) %in% cells_y_loss, "y_loss", "no_y_loss")

day10_cells$y_loss <- ifelse(rownames(day10_cells@meta.data) %in% cells_y_loss.10, "y_loss", "no_y_loss")

# Check the metadata
table(day5_cells$y_loss)
table(day10_cells$y_loss)

# Define the new metadata labels
day5_cells$new_metadata <- ifelse(
  day5_cells$recoded_loss_chr3_2_bin == "02_chr3_loss" & day5_cells$y_loss == "y_loss", "3p_and_y_loss",
  ifelse(day5_cells$recoded_loss_chr3_2_bin == "02_chr3_loss" & day5_cells$y_loss == "no_y_loss", "3p_loss_only",
         ifelse(day5_cells$recoded_loss_chr3_2_bin == "01_No_loss" & day5_cells$y_loss == "y_loss", "y_loss_only", 
                "no_loss"))
)

# Check the distribution of the new metadata labels
table(day5_cells$new_metadata)
#3 p_and_y_loss  3p_loss_only       no_loss   y_loss_only 
# 31             883               11875      662
# Similarly for day10_cells
day10_cells$new_metadata <- ifelse(
  day10_cells$recoded_loss_chr3_2_bin == "02_chr3_loss" & day10_cells$y_loss == "y_loss", "3p_and_y_loss",
  ifelse(day10_cells$recoded_loss_chr3_2_bin == "02_chr3_loss" & day10_cells$y_loss == "no_y_loss", "3p_loss_only",
         ifelse(day10_cells$recoded_loss_chr3_2_bin == "01_No_loss" & day10_cells$y_loss == "y_loss", "y_loss_only", 
                "no_loss"))
)

# Check the distribution for day 10
table(day10_cells$new_metadata)
# 3p_and_y_loss  3p_loss_only       no_loss   y_loss_only 
# 41             151                8907      2570 

# Validate cell groupings
Idents(day5_cells) <- "new_metadata"
VlnPlot(day5_cells, features = c("EIF1AY", "KDM5D", "RPS4Y2", "TMSB4Y", "RPS4Y1"))

Idents(day10_cells) <- "new_metadata"
VlnPlot(day10_cells, features = c("EIF1AY", "KDM5D", "RPS4Y2", "TMSB4Y", "RPS4Y1"))

# DEG analysis for day5 cells: Chromosome 3 loss vs Y loss
Idents(day5_cells) <- "new_metadata"
day5_cells <- PrepSCTFindMarkers(day5_cells)
day5_DEGs.y <- FindMarkers(day5_cells, ident.1 = "3p_loss_only", ident.2 = "y_loss_only", verbose = TRUE, logfc.threshold = 0.05)
# View top results
head(day5_DEGs.y)
write.csv(day5_DEGs.y, file = "DEG_day5_3p_vs_Yq.csv")

# DEG analysis for day10 cells: Chromosome 3 loss vs Y loss
Idents(day10_cells) <- "new_metadata"
day10_cells <- PrepSCTFindMarkers(day10_cells)
day10_DEGs.y <- FindMarkers(day10_cells, ident.1 = "3p_loss_only", ident.2 = "y_loss_only", verbose = TRUE, logfc.threshold = 0.05)
write.csv(day10_DEGs.y, file = "DEG_day10_3p_vsYq.csv")

# For day5 DEGs
day5_DEGs.y$log2FC <- day5_DEGs.y$avg_log2FC
day5_DEGs.y$neg_log10_pval <- -log10(day5_DEGs.y$p_val_adj)

# For day10 DEGs
day10_DEGs.y$log2FC <- day10_DEGs.y$avg_log2FC
day10_DEGs.y$neg_log10_pval <- -log10(day10_DEGs.y$p_val_adj)

# Volcano plot for day5 DEGs
ggplot(day5_DEGs.y, aes(x = log2FC, y = neg_log10_pval)) +
  geom_point(aes(color = (p_val_adj < 0.05 & abs(log2FC) > 0.25)), size = 2) +
  scale_color_manual(values = c("grey", "blue")) +  # Customize colors
  theme_minimal() +                                 # Use a minimal theme
  labs(title = "DEGs: Chromosome 3 Loss vs Y loss (Day 5)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.position = "none") +                 # Remove legend
  geom_text_repel(data = subset(day5_DEGs.y, p_val_adj < 0.05 & abs(log2FC) > 0.5), # Add labels for significant DEGs
                  aes(label = rownames(subset(day5_DEGs.y, p_val_adj < 0.05 & abs(log2FC) > 0.5))),
                  size = 3, max.overlaps = 15)      # Control overlap and label size

# Volcano plot for day10 DEGs
ggplot(day10_DEGs.y, aes(x = log2FC, y = neg_log10_pval)) +
  geom_point(aes(color = (p_val_adj < 0.05 & abs(log2FC) > 0.25)), size = 2) +
  scale_color_manual(values = c("grey", "red")) +  # Customize colors
  theme_minimal() +                                 # Use a minimal theme
  labs(title = "DEGs: Chromosome 3 Loss vs Y loss (Day 10)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.position = "none") +                 # Remove legend
  geom_text_repel(data = subset(day10_DEGs.y, p_val_adj < 0.05 & abs(log2FC) > 0.5), # Add labels for significant DEGs
                  aes(label = rownames(subset(day10_DEGs.y, p_val_adj < 0.05 & abs(log2FC) > 0.5))),
                  size = 3, max.overlaps = 20)      # Control overlap and label size


###Proportion of cells with 3p loss (out of total) in both timepoints
library(ggplot2)
library(dplyr)

# For day 5 cells: Calculate percentages
day5_population <- as.data.frame(table(day5_cells$new_metadata)) %>%
  mutate(percentage = Freq / sum(Freq) * 100,
         timepoint = "day5")

# For day 10 cells: Calculate percentages
day10_population <- as.data.frame(table(day10_cells$new_metadata)) %>%
  mutate(percentage = Freq / sum(Freq) * 100,
         timepoint = "day10")

# Combine both datasets
combined_population <- rbind(day5_population, day10_population)

# Rename columns for clarity
colnames(combined_population) <- c("Group", "Count", "Percentage", "Timepoint")

# Set the order of the 'Timepoint' factor so that day5 comes before day10
combined_population$Timepoint <- factor(combined_population$Timepoint, levels = c("day5", "day10"))

# Print the table
print(combined_population)
# Export combined_population as a CSV file
write.csv(combined_population, file = "combined_3p_yq_counts.csv", row.names = FALSE)

# Plot the stacked bar chart
ggplot(combined_population, aes(x = Timepoint, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Proportion of Cell Populations at Day 5 and Day 10",
       x = "Timepoint", y = "Percentage of Cells") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")  # Choose a color palette for the groups


### Venn diagram for comparing Yq and 3p loss cells
### 1. Perform DEG for each population against no loss cells individually


# Ensure the plotting area is cleared
dev.off()

# Set the significance threshold for p-value and log2 fold change
pval_threshold <- 0.1
logfc_threshold <- 0.5

# Filter upregulated genes for day 5 (log2FC > 0.25 and p_val_adj < 0.05)
day5_upregulated <- rownames(subset(day5_DEGs, avg_log2FC > logfc_threshold & p_val_adj < pval_threshold))

# Filter upregulated genes for day 10 (log2FC > 0.25 and p_val_adj < 0.05)
day10_upregulated <- rownames(subset(day10_DEGs, avg_log2FC > logfc_threshold & p_val_adj < pval_threshold))

# Use the VennDiagram package to create a Venn diagram
library(VennDiagram)

# Create a list of the upregulated genes for comparison
upregulated_list <- list("Day 5 Upregulated" = day5_upregulated, "Day 10 Upregulated" = day10_upregulated)

# Create the Venn diagram
venn.plot <- venn.diagram(
  x = upregulated_list,
  category.names = c("Day 5 Upregulated", "Day 10 Upregulated"),
  fill = c("blue", "red"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  filename = NULL
)

# Plot the Venn diagram
grid.draw(venn.plot)

# Step 4: Find unique and intersecting genes
# Unique to day 5
unique_day5 <- setdiff(day5_upregulated, day10_upregulated)

# Unique to day 10
unique_day10 <- setdiff(day10_upregulated, day5_upregulated)

# Shared between day 5 and day 10
shared_genes <- intersect(day5_upregulated, day10_upregulated)

# View the results
unique_day5
unique_day10
shared_genes

# If you want to save the gene lists as text files:
write.table(unique_day5, "unique_day5_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(unique_day10, "unique_day10_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(shared_genes, "shared_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

#### Identify Yq loss in the big integrated object + perform Venn diagram

# These are the genes of interest to identify Y chromosome loss
y_loss_genes <- c("EIF1AY", "KDM5D", "RPS4Y2", "TMSB4Y")
rps4y1_gene <- "RPS4Y1"

# Define loss threshold: cells with expression below this value for Y genes are considered "lost"
loss_threshold <- 0.1

# Fetch the expression levels of the Y loss genes and RPS4Y1 for the entire object (obj)
expr_y_loss <- FetchData(obj, vars = y_loss_genes)
expr_rps4y1 <- FetchData(obj, vars = rps4y1_gene)

# Identify cells with Y loss: below threshold for all Y loss genes and above threshold for RPS4Y1
cells_y_loss <- rownames(expr_y_loss)[
  rowSums(expr_y_loss < loss_threshold) == length(y_loss_genes) &  # Loss of Y genes
    expr_rps4y1 > loss_threshold  # Retained RPS4Y1 expression
]

# Add 3p loss information (already in the metadata slot 'recoded_loss_chr3_2_bin')
# Create a new metadata column that classifies based on both 3p and Y loss
obj$new_metadata <- ifelse(rownames(obj) %in% intersect(rownames(obj), cells_y_loss) & 
                             obj$recoded_loss_chr3_2_bin == "02_chr3_loss", "3p_and_y_loss",
                           ifelse(rownames(obj) %in% intersect(rownames(obj), cells_y_loss), "y_loss_only",
                                  ifelse(obj$recoded_loss_chr3_2_bin == "02_chr3_loss", "3p_loss_only", "no_loss")))

# Fetch the metadata related to 3p loss into a data frame
metadata <- obj@meta.data  # Extract metadata

# Create a new column for identifying 3p and Y loss statuses
metadata$new_metadata <- "no_loss"  # Default classification as no loss

# Now update based on Y loss and 3p loss information
metadata$new_metadata[rownames(metadata) %in% cells_y_loss & metadata$recoded_loss_chr3_2_bin == "02_chr3_loss"] <- "3p_and_y_loss"
metadata$new_metadata[rownames(metadata) %in% cells_y_loss & metadata$recoded_loss_chr3_2_bin != "02_chr3_loss"] <- "y_loss_only"
metadata$new_metadata[rownames(metadata) %in% rownames(metadata) & metadata$recoded_loss_chr3_2_bin == "02_chr3_loss" & !(rownames(metadata) %in% cells_y_loss)] <- "3p_loss_only"

# Assign the modified metadata back to the Seurat object
obj$new_metadata <- metadata$new_metadata
# Check the distribution of the new classifications
table(obj$new_metadata)
#3p_and_y_loss  3p_loss_only       no_loss   y_loss_only 
#137           969         20023          3991
# Set identity based on the new metadata
Idents(obj) <- "new_metadata"

# Differential expression for 3p loss cells compared to no loss
deg_3p_loss <- FindMarkers(obj, ident.1 = "3p_loss_only", ident.2 = "no_loss", logfc.threshold = 0.25, min.pct = 0.1)

# Differential expression for Y loss cells compared to no loss
deg_y_loss <- FindMarkers(obj, ident.1 = "y_loss_only", ident.2 = "no_loss", logfc.threshold = 0.25, min.pct = 0.1)

deg_3p_loss_vs_y <- FindMarkers(obj, ident.1 = "3p_loss_only", ident.2 = "y_loss_only", logfc.threshold = 0.25, min.pct = 0.1)

write.csv(deg_3p_loss, "deg_3p_loss_vs_no_loss.csv")
write.csv(deg_y_loss, "deg_y_loss_vs_no_loss.csv")
write.csv(deg_3p_loss_vs_y, "deg_3p_loss_vs_y_loss.csv")

# Filter genes for 3p loss with adj. p-value < 0.05 and log2FC > 0.25
significant_3p_loss_genes <- rownames(subset(deg_3p_loss, p_val_adj < 0.05 & avg_log2FC > 0.25))

# Filter genes for Y loss with adj. p-value < 0.05 and log2FC > 0.25
significant_y_loss_genes <- rownames(subset(deg_y_loss, p_val_adj < 0.05 & avg_log2FC > 0.25))

# Load the VennDiagram package if not already loaded
library(VennDiagram)

# Create a list of the significant genes for 3p loss and Y loss
gene_list <- list("3p Loss Genes" = significant_3p_loss_genes, "Y Loss Genes" = significant_y_loss_genes)

# Find intersected genes between 3p loss and Y loss
intersected_genes <- intersect(significant_3p_loss_genes, significant_y_loss_genes)

# Find unique genes for 3p loss and Y loss
unique_3p_loss_genes <- setdiff(significant_3p_loss_genes, significant_y_loss_genes)
unique_y_loss_genes <- setdiff(significant_y_loss_genes, significant_3p_loss_genes)


# Generate the Venn diagram

# Create the Venn diagram with labels for intersected genes
# Load the VennDiagram package
library(VennDiagram)

# Create the Venn diagram and customize the intersect region with gene names
# Load the VennDiagram package
library(VennDiagram)
library(grid)

# Create the Venn diagram with basic settings
venn.plot <- venn.diagram(
  x = gene_list,
  category.names = c("3p Loss Genes", "Y Loss Genes"),
  fill = c("purple", "gray"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 1.5,  # Adjust size of category labels
  filename = NULL,
  main = "Venn Diagram of 3p Loss vs Y Loss Genes",
  label.col = "black",  # Set color for all labels
  fontface = "bold",
  lwd = 2,  # Line width for circles
  cat.pos = c(-20, 20),  # Adjust position of labels
  cat.dist = c(0.05, 0.05)  # Adjust distance of labels from circles
)

# Draw the Venn diagram
grid.draw(venn.plot)

# Overlay the actual intersecting gene names in the intersect region
grid.text(
  paste(intersected_genes, collapse = "\n"),  # Join the gene names with line breaks
  x = 0.5,  # Horizontal position (adjust as needed)
  y = 0.5,  # Vertical position (adjust as needed)
  gp = gpar(fontsize = 12, col = "black", fontface = "bold")
)


# Annotate the intersected genes on the plot (replacing the count number)
grid::grid.text(paste(intersected_genes, collapse = "\n"), x = 0.5, y = 0.5, gp = grid::gpar(fontsize = 10, col = "black"))


# Export the unique and intersected genes to text files
write.table(unique_3p_loss_genes, file = "unique_3p_loss_genes_vs_y.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(unique_y_loss_genes, file = "unique_y_loss_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(intersected_genes, file = "intersected_genes_3p_vs_y.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

Idents(obj) <- "new_metadata"
VlnPlot(obj, features = c("EIF1AY", "KDM5D", "RPS4Y2", "TMSB4Y", "RPS4Y1"))

day5_DEGs <- read.csv("3pLoss_vs_NoLoss_day5_DEGs.csv", header = TRUE, stringsAsFactors = FALSE)

day10_DEGs <- read.csv("3pLoss_vs_NoLoss_day10_DEGs.csv", header = TRUE, stringsAsFactors = FALSE)

# Filter the day 5 data for upregulated genes (log2FC > 0.25 and adj p < 0.05)
day5_upregulated <- day5_DEGs %>%
  filter(avg_log2FC > 0.25 & p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))  # Order by log2 fold change for display purposes

# Similarly, filter the day 10 data
day10_upregulated <- day10_DEGs %>%
  filter(avg_log2FC > 0.25 & p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))



# Filter to keep the top 100 genes based on -log10(p_val_adj)
top_genes_5 <- day5_upregulated %>%
  arrange(desc(-log10(p_val_adj))) %>%  # Order by -log10(p_val_adj) in descending order
  head(100)  # Take the top 100 genes

top_genes_10 <- day10_upregulated %>%
  arrange(desc(-log10(p_val_adj))) %>%  # Order by -log10(p_val_adj) in descending order
  head(100)  # Take the top 100 genes

# Plot with bars and color gradient for log2 Fold Change
library(ggplot2)

ggplot(top_genes_5, aes(x = reorder(X, -log10(p_val_adj)), y = -log10(p_val_adj), fill = avg_log2FC)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_gradient(low = "lightgray", high = "blue", name = "Log2 Fold Change") +
  labs(title = "Day 5 Upregulated Genes (3p Loss vs No Loss)",
       x = "Genes", y = "-log10 Adjusted P-Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip()  # Flip coordinates to make it easier to read gene names

ggplot(top_genes_10, aes(x = reorder(X, -log10(p_val_adj)), y = -log10(p_val_adj), fill = avg_log2FC)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_gradient(low = "lightgray", high = "red", name = "Log2 Fold Change") +
  labs(title = "Day 10 Upregulated Genes (3p Loss vs No Loss)",
       x = "Genes", y = "-log10 Adjusted P-Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip()  # Flip coordinates to make it easier to read gene names

### remove labels from left side of DEG volcano plots

# Assuming day5_DEGs is the data frame containing your results
# Calculate negative log10 adjusted p-value
day5_DEGs$neg_log10_pval <- -log10(day5_DEGs$p_val_adj)
day10_DEGs$neg_log10_pval <- -log10(day10_DEGs$p_val_adj)

# Filter out rows with NA or zero p-values to prevent issues with log transformation
day5_DEGs <- day5_DEGs[!is.na(day5_DEGs$neg_log10_pval) & 
                         day5_DEGs$p_val_adj > 0 & 
                         is.finite(day5_DEGs$neg_log10_pval), ]
day10_DEGs <- day10_DEGs[!is.na(day10_DEGs$neg_log10_pval) & 
                         day10_DEGs$p_val_adj > 0 & 
                         is.finite(day10_DEGs$neg_log10_pval), ]

library(ggplot2)
library(ggrepel)


# Create the volcano plot using day5_DEGs
ggplot(day5_DEGs, aes(x = avg_log2FC, y = neg_log10_pval)) +
  geom_point(aes(color = ifelse(p_val_adj < 0.05 & avg_log2FC > 0.25, "upregulated", 
                                ifelse(p_val_adj < 0.05 & avg_log2FC < -0.25, "downregulated", "not_significant"))), 
             size = 2) +
  scale_color_manual(values = c("upregulated" = "blue", "downregulated" = "grey", "not_significant" = "lightgrey")) +  # Customize colors
  theme_minimal() +                                 # Use a minimal theme
  labs(title = "DEGs: Chromosome 3 Loss vs Y Loss (Day 5)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.position = "none") +                 # Remove legend
  geom_text_repel(data = subset(day5_DEGs, p_val_adj < 0.05 & avg_log2FC > 0.25),  # Subset for labels
                  aes(label = X),  # Use the 'X' column for gene names
                  size = 3, max.overlaps = 20) +
  coord_cartesian(ylim = c(0, max(day5_DEGs$neg_log10_pval, na.rm = TRUE) + 1))  # Adjust y-axis limits

## day10
ggplot(day10_DEGs, aes(x = avg_log2FC, y = neg_log10_pval)) +
  geom_point(aes(color = ifelse(p_val_adj < 0.05 & avg_log2FC > 0.25, "upregulated", 
                                ifelse(p_val_adj < 0.05 & avg_log2FC < -0.25, "downregulated", "not_significant"))), 
             size = 2) +
  scale_color_manual(values = c("upregulated" = "red", "downregulated" = "grey", "not_significant" = "lightgrey")) +  # Customize colors
  theme_minimal() +                                 # Use a minimal theme
  labs(title = "DEGs: Chromosome 3 Loss vs Y Loss (Day 10)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.position = "none") +                 # Remove legend
  geom_text_repel(data = subset(day10_DEGs, p_val_adj < 0.05 & avg_log2FC > 0.25),  # Subset for labels
                  aes(label = X),  # Use the 'X' column for gene names
                  size = 3, max.overlaps = 20) +
  coord_cartesian(ylim = c(0, max(day10_DEGs$neg_log10_pval, na.rm = TRUE) + 1))  # Adjust y-axis limits










