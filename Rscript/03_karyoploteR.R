##========================================##
##    DE and karyoploteR analysis         ##
##========================================##
#Load Packages

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(karyoploteR)
library(tidyverse)
library(future)
library(Matrix)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(readxl)
library(ggrepel)


plan()
plan(multicore, workers = 22) # make the computer has enough cores. Use carefully. 
plan()
options(future.globals.maxSize = 50 * 1024^3)  # allocate 50 GB RAM
options(scipen = 100) # set to avoid using scientific notation. 

base_dir <- "/path/to/Chromothripsis_2"  # <-- set to your local path to the project root
figures_dir <- file.path(base_dir, "Figures/integrated/KaryoploteR")

# load data
obj <- readRDS(file=file.path(base_dir, "/DATA/Integrated/InferCNV/Post_inferCNV_obj.rds"))
colnames(obj@meta.data)

# extract annnotation
HMM_anno <- colnames(obj@meta.data[, 18:343])
Chr3_changes <- c(colnames(obj@meta.data[, 36:44]))
Chr3_changes               

meta <- obj@meta.data
           
pdf(file.path(figures_dir, "Featureplot_Chr3_changes.pdf"),
    width=12,height=10,paper='special')
FeaturePlot(obj, features = Chr3_changes)
dev.off()


######################################################
#    select cells with highest chr3 changes          # 
######################################################
Idents(obj) <- "seurat_clusters"

pdf(file.path(figures_dir, "DimPlot.pdf"),width=7,height=6,paper='special')
DimPlot(obj)
dev.off()

Idents(obj) <- "proportion_loss_chr3"

pdf(file.path(figures_dir, "DimPlot_chr3loss_HTO.pdf"),width=12,height=8,paper='special')
DimPlot(obj, split.by = "HTO_classification", ncol = 4)
dev.off()

pdf(file.path(figures_dir, "Hist_chr3loss.pdf"),width=6,height=4,paper='special')
hist(meta$proportion_scaled_loss_chr3, breaks = 100)
dev.off()

table(meta$proportion_scaled_loss_chr3)

# Set identities to "proportion_loss_chr3"
Idents(obj) <- "proportion_loss_chr3"

# Retrieve all unique clusters
clusters <- levels(Idents(obj))

# Define custom colors: "0" = gray, "1" = red
custom_colors <- c("0" = "gray", "1" = "red")

# Number of remaining levels to color
num_other_levels <- 17

# Generate a colorblind-friendly palette for the other levels
# Use the 'Set3' palette from RColorBrewer, which has 12 colors, and then use additional colors from 'Paired'
additional_colors <- c(brewer.pal(12, "Set3"), brewer.pal(9, "Paired")[1:(num_other_levels - 12)])

# Combine custom colors with generated colors
all_colors <- c(custom_colors, additional_colors)

# Assign colors to levels
names(all_colors) <- clusters

# Open the PDF device
pdf(output_pdf, width = 12, height = 8, paper = 'special')

# Create the DimPlot with custom colors, jitter, and adjusted point size
DimPlot(
  obj,
  split.by = "HTO_classification",
  ncol = 4,
  cols = all_colors,
  pt.size = 5,  # Increase point size for better visibility
  raster = TRUE  # Ensure higher quality plotting
) +
  ggtitle("DimPlot of chr3 Loss by HTO Classification") +
  theme_minimal() +
  # Set transparency directly in the plot
  scale_alpha_manual(values = 0.4)  # Adjust alpha for all points

# Close the PDF device
dev.off()


# === recode : define chr3_loss range ===
# Define the breaks and labels for the new categories
# Define breakpoints for the four categories
breaks <- c(-Inf, 0, 0.33, 0.67, Inf)  # Adjust these values as needed
labels = c("No_loss", "Low_loss", "Medium-High_loss", "High_loss")

# Use the cut function to categorize the proportion_loss_chr3 values
meta$recoded_loss_chr3 <- cut(meta$proportion_loss_chr3,
                          breaks = breaks,
                          labels = labels,
                          include.lowest = TRUE)

# Check the categorization
table(meta$recoded_loss_chr3)

# Create a dataframe for shading the categories
shade_df <- data.frame(
  xmin = breaks[-length(breaks)],
  xmax = breaks[-1],
  category = labels
)

# Filter out the rows where proportion_loss_chr3 is 0
meta_filtered <- meta %>% filter(proportion_loss_chr3 != 0)

# Generate the histogram with the filtered data

output_pdf <- file.path(figures_dir, "Hist_chr3_loss_category.pdf")
pdf(output_pdf, width = 12, height = 8, paper = 'special')
ggplot(meta_filtered, aes(x = proportion_loss_chr3)) +
  # Shaded rectangles for categories
  geom_rect(data = shade_df, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = category),
            alpha = 0.2) +
  # Histogram bars
  geom_histogram(binwidth = 0.01, color = "black", fill = "white", boundary = 0) +
  # Vertical lines at breakpoints
  geom_vline(xintercept = breaks[-c(1, length(breaks))], linetype = "dashed", color = "red") +
  # Labels and title
  labs(
    title = "Histogram of Proportion Loss chr3 with Category Boundaries",
    x = "Proportion Loss chr3",
    y = "Frequency",
    fill = "Category"
  ) +
  # Custom fill colors (optional)
  scale_fill_brewer(palette = "Pastel1") +
  # Theme adjustments
  theme_minimal()
dev.off()

# write back to obj
obj@meta.data <- meta

# plot by category. 
output_pdf <- file.path(figures_dir, "DimPlot_Chr3_loss_category.pdf")
pdf(output_pdf, width = 10, height = 8, paper = 'special')
DimPlot(
  obj,
  split.by = "recoded_loss_chr3",
  ncol = 2,
  cols = all_colors,
  pt.size = 5,  # Increase point size for better visibility
  raster = TRUE  # Ensure higher quality plotting
  ) +
  ggtitle("DimPlot of chr3 Loss by Chr3_loss_category") +
  theme_minimal() +
  # Set transparency directly in the plot
  scale_alpha_manual(values = 0.4)  # Adjust alpha for all points
dev.off()

# Examine Chr3p genes
DefaultAssay(obj) <- "SCT"
Chr_3p_genes_list <- read_excel(file.path(base_dir, "Misc/Chr.3p_genes_list.xlsx"))
View(Chr_3p_genes_list)

Chr_3p_genes <- Chr_3p_genes_list$`Gene name`

pdf(file.path(figures_dir, "FeaturePlot_3ploss_markers.pdf"),width=12,height=12,paper='special')
FeaturePlot(obj, features = c("GBE1", "ROBO1", "ROBO2", "ZNF717"))
dev.off()




# ========== Global DE gene analysis between chr3 or ch5 CNV TRUE / FALSE group. =============
# set a volcano plot function 

library(ggrepel)
create_volcano_plot <- function(data, logFC_column, p_value_column, title, p_value_threshold) {
  # Add gene names from row names
  data$gene <- rownames(data)
  
  # Transform the p-value and determine significance
  data$neg_log10_p_val_adj <- -log10(data[[p_value_column]])
  data$significant <- data[[p_value_column]] < p_value_threshold & abs(data[[logFC_column]]) > 1
  
  # Create the plot
  plot <- ggplot(data, aes_string(x = logFC_column, y = "neg_log10_p_val_adj")) +
    geom_point(aes(color = significant), alpha = 0.5) +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
    geom_text_repel(data = subset(data, significant),
                    aes(label = gene),
                    max.overlaps = Inf,
                    box.padding = 0.35,
                    point.padding = 0.5,
                    segment.color = 'grey50') +
    ggtitle(title) +
    xlab("Average Log2 Fold Change") +
    ylab("-log10(Adjusted P-value)") +
    geom_hline(yintercept = -log10(p_value_threshold), color = "red", linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), color = "blue", linetype = "dashed") +
    theme_minimal() +
    guides(color = FALSE)  # Optionally hide the legend
  
  return(plot)
}

table(meta$recoded_loss_chr3)

# recoded 
# Recode 'recoded_loss_chr3' into 'recoded_loss_chr3_2_bin'
obj$recoded_loss_chr3_2_bin <- ifelse(
  obj$recoded_loss_chr3 == "No_loss",
  "01_No_loss",
  "02_chr3_loss"
)

# Convert the new column to a factor with specified levels
obj$recoded_loss_chr3_2_bin <- factor(
  obj$recoded_loss_chr3_2_bin,
  levels = c("01_No_loss", "02_chr3_loss")
)

## Identify global DE genes by comparing cell with/without Chr3.loss
Idents(obj) <- "recoded_loss_chr3_2_bin"
Idents(obj) <- factor(x = Idents(obj), levels = sort(levels(obj)))

pdf(file.path(figures_dir, "VlnPlot_obj_VHL.pdf"),
    width=6,height=4,paper='special')
VlnPlot(obj, feature = "VHL")
dev.off()

pdf(file.path(figures_dir, "VlnPlot_obj_GBE1.pdf"),
    width=6,height=4,paper='special')
VlnPlot(obj, feature = "GBE1")
dev.off()

pdf(file.path(figures_dir, "VlnPlot_obj_SQSTM1.pdf"),
    width=6,height=4,paper='special')
VlnPlot(obj, feature = "SQSTM1")
dev.off()

table(Idents(obj))
DefaultAssay(obj) <- "SCT"
obj <- PrepSCTFindMarkers(obj)

DEG_Chr3.loss <- FindMarkers(obj, ident.1 = c('02_chr3_loss'), ident.2 = c('01_No_loss'), min.pct = 0.25)
write.csv(DEG_Chr3.loss, file = file.path(figures_dir, "DEG_Chr3.loss_DEG.csv"))

# filter data before valcano plot to perserve the p-value = 0 genes on volcano
# Calculate min_p_value
min_p_value <- min(DEG_Chr3.loss$p_val_adj[DEG_Chr3.loss$p_val_adj > 0], na.rm = TRUE)

# Create new column
DEG_Chr3.loss$p_val_adj_non_zero <- ifelse(DEG_Chr3.loss$p_val_adj == 0, 
                                           min_p_value / 10, 
                                           DEG_Chr3.loss$p_val_adj)


pdf(file.path(figures_dir, "Volcano_obj.DEG.Chr3.loss_noZero_adjP.pdf"),
    width=10,height=10,paper='special')
create_volcano_plot(DEG_Chr3.loss, "avg_log2FC", "p_val_adj_non_zero", "DEG_Chr3.loss (right: increase in Chr3.loss)", 0.05)
dev.off()

pdf(file.path(figures_dir, "VlnPlot_obj_RBMS3.pdf"),
    width=4,height=4,paper='special')
VlnPlot(obj, feature = "RBMS3")
dev.off()

pdf(file.path(figures_dir, "VlnPlot_obj_OSBPL10.pdf"),
    width=4,height=4,paper='special')
VlnPlot(obj, feature = "OSBPL10")
dev.off()

pdf(file.path(figures_dir, "VlnPlot_obj_GBE1.pdf"),
    width=4,height=4,paper='special')
VlnPlot(obj, feature = "GBE1")
dev.off()

# save current working env
save.image("./DATA/Integrated/KaryoploteR/03_3p_loss.RData")
# To restore work environment later use: load("./DATA/Integrated/KaryoploteR/03_3p_loss.RData") 

## ======== plot DEG on the Chr3 and Chr5 using karayoploter package ===========
# check the tutorial: https://bernatgel.github.io/karyoploter_tutorial//Examples/GeneExpression/GeneExpression.html 

# Load the libraries
library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
hg.genes <- genes(txdb)

# map gene symbol to geneID
# Convert the gene IDs to a character vector
gene_ids <- as.character(mcols(hg.genes)$gene_id)

# Map gene IDs to symbols
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = gene_ids, 
                       column = "SYMBOL", 
                       keytype = "ENTREZID", 
                       multiVals = "first")

mcols(hg.genes)$symbol <- gene_symbols[as.character(mcols(hg.genes)$gene_id)]
head(hg.genes, n=50)

# ======= DEG list of  Chr3.loss T/F comparison  ====== 
DEG_Chr3.loss$symbol <- row.names(DEG_Chr3.loss)
head(DEG_Chr3.loss)

# Merge DEG data into GRanges by matching "symbol"
# Create a named list from the DEG table, indexed by the symbol
DEG_Chr3.loss_list <- split(DEG_Chr3.loss, DEG_Chr3.loss$symbol)

# Function to fetch DEG data by symbol
fetch_deg_data <- function(symbol) {
  if (symbol %in% names(DEG_Chr3.loss_list)) {
    return(DEG_Chr3.loss_list[[symbol]])
  } else {
    # Return NA values if symbol is not found (adjust the length according to the number of DEG columns)
    return(rep(NA, length(DEG_Chr3.loss_list[[1]]) - 1))
  }
}

# Add DEG data as metadata columns to the GRanges object
mcols(hg.genes) <- cbind(mcols(hg.genes), do.call(rbind, lapply(mcols(hg.genes)$symbol, fetch_deg_data)))

# The GRanges object now includes DEG data as additional metadata columns
head(hg.genes, n=50)

# establish karyotype plot
kp <- plotKaryotype(plot.type=2, objsomes = c("chr3"))
kpDataBackground(kp, data.panel = 1, col="#AACBFF")
kpDataBackground(kp, data.panel = 2, col="#FFAACB")

# sort gene by avg_log2FC and select top 50 genes. 
ordered <- hg.genes[order(hg.genes$avg_log2FC, na.last = TRUE),]

top.genes <- ordered[1:50]

# map names to symbol
symbol_vector <- mcols(hg.genes)$symbol

# Ensure names are set to allow direct indexing by symbol
names(symbol_vector) <- symbol_vector

# Directly extract the top 50 gene symbols from the ordered 'hg.genes'
ordered_symbols <- mcols(ordered)[1:50, "symbol"]
ordered_symbols

kp <- kpPlotMarkers(kp, ordered[1:50], labels = ordered_symbols, text.orientation = "horizontal")


filtered.hg.genes <- ordered[!is.na(ordered$p_val_adj)]
head(filtered.hg.genes)
log.pval <- -log10(filtered.hg.genes$p_val_adj)
mcols(filtered.hg.genes)$log.pval <- log.pval
filtered.hg.genes
head(filtered.hg.genes)

sign.genes <- filtered.hg.genes[filtered.hg.genes$p_val_adj < 0.01,]
head(sign.genes)
kp <- plotKaryotype(plot.type=2, objsomes = c("chr3"))
kpDataBackground(kp, data.panel = 1, col="#AACBFF")
kpPoints(kp, data=sign.genes, y=sign.genes$log.pval, ymax=max(sign.genes$log.pval))
kp <- kpPlotMarkers(kp, ordered[1:50], labels = ordered_symbols, text.orientation = "vertical")

kp <- plotKaryotype(plot.type=2, objsomes = c("chr3"))
kpDataBackground(kp, data.panel = 1, col="#AACBFF")
fc.ymax <- ceiling(max(abs(range(sign.genes$avg_log2FC))))
fc.ymin <- -fc.ymax
kpPoints(kp, data=sign.genes, y=sign.genes$avg_log2FC, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin)
kpAddLabels(kp, labels = "avg_log2FC", srt=90, pos=1, label.margin = 0.06, ymax=fc.ymax, ymin=fc.ymin)
kp <- kpPlotMarkers(kp, ordered[1:50], labels = ordered_symbols, text.orientation = "vertical")

# plot significant DE genes across genome

# Set fold change and p-value thresholds
log2FC_threshold <- 1.5  # Adjust this value as needed
pval_threshold <- 0.05  # Example threshold for adjusted p-value

# Filter genes based on these thresholds
filtered_genes <- subset(sign.genes, abs(avg_log2FC) > log2FC_threshold & p_val_adj < pval_threshold)

# Recalculate y-axis limits based on the filtered data
fc.ymax <- ceiling(max(abs(range(filtered_genes$avg_log2FC))))
fc.ymin <- -fc.ymax

# Plot the filtered genes
pdf(file.path(figures_dir, "Karyotype_obj.DEG.Chr3_filtered.pdf"),
    width = 30, height = 30, paper = 'special')
kp <- plotKaryotype(plot.type = 2)
kpDataBackground(kp, data.panel = 1, col = "#AACBFF")
kpPoints(kp, data = filtered_genes, y = filtered_genes$avg_log2FC, ymax = fc.ymax, ymin = fc.ymin)
kpAxis(kp, ymax = fc.ymax, ymin = fc.ymin)
kpAddLabels(kp, labels = "avg_log2FC", srt = 90, pos = 1, label.margin = 0.02, ymax = fc.ymax, ymin = fc.ymin)
dev.off()


#  add p value to the plot. The size of the dot represent the log p-value. Bigger the size = more significant
cex.val <- sqrt(sign.genes$log.pval)/3
kp <- plotKaryotype(genome="hg38")
kpPoints(kp, data=filtered_genes, y=sign.genes$avg_log2FC, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin)
kpAddLabels(kp, labels = "avg_log2FC", srt=90, pos=1, label.margin = 0.02, ymax=fc.ymax, ymin=fc.ymin)

# add top 50 gene symbol
pdf(file.path(figures_dir, "Karyotype_obj.DEG.Chr3_top50.pdf"),
    width=30,height=30,paper='special')

kp <- plotKaryotype(genome="hg38")
kpPoints(kp, data=filtered_genes, y=sign.genes$avg_log2FC, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin)
kpAddLabels(kp, labels = "avg_log2FC", srt=90, pos=1, label.margin = 0.02, ymax=fc.ymax, ymin=fc.ymin)
kpPlotMarkers(kp, top.genes, labels = ordered_symbols, text.orientation = "vertical")

dev.off()

# Prep for dual data panel plot
points.top <- 0.8

col.over <- "#FFBD07AA"
col.under <- "#00A6EDAA"
sign.col <- rep(col.over, length(sign.genes))
sign.col[sign.genes$avg_log2FC<0] <- col.under

# Data panel 1
kp <- plotKaryotype(genome="hg38", plot.type=2)
kpPoints(kp, data=sign.genes, y=sign.genes$avg_log2FC, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpAddLabels(kp, labels = "avg_log2FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, 
           chr=as.character(seqnames(top.genes)), 
           x0=gene.mean, 
           x1=gene.mean, 
           y0=top.genes$avg_log2FC,  # Starting y-values for each segment
           y1=fc.ymax,  # Assuming this is a constant or a vector of ending y-values for the segments
           ymax=fc.ymax, 
           ymin=fc.ymin, 
           r1=points.top, 
           col="#777777")

kpPlotMarkers(kp, top.genes, labels = ordered_symbols, text.orientation = "vertical", r0=points.top, cex=0.5)

# Data panel 2
kp <- kpPlotDensity(kp, data=hg.genes, window.size = 10e4, data.panel = 2)



##  ===  Final plot (all chr) ==== 
pp <- getDefaultPlotParams(plot.type = 2)
pp$data2height <- 100
pp$ideogramheight <- 100

pdf(file.path(figures_dir, "Karyotype_obj.DEG.Chr3_top50_Allobj.pdf"),
    width=30,height=30,paper='special')

kp <- plotKaryotype(genome="hg38", plot.type=2, plot.params = pp)
kpAddMainTitle(kp, main = "Gene expression - Chr3 loss True vs False")
## Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$avg_log2FC, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2

kpSegments(kp, 
           chr=as.character(seqnames(top.genes)), 
           x0=gene.mean, 
           x1=gene.mean, 
           y0=top.genes$avg_log2FC,  # Starting y-values for each segment
           y1=fc.ymax,  # Assuming this is a constant or a vector of ending y-values for the segments
           ymax=fc.ymax, 
           ymin=fc.ymin, 
           r1=points.top, 
           col="#777777")

kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin)
kpAddLabels(kp, labels = "avg_log2FC", srt=90, pos=1, label.margin = 0.02, ymax=fc.ymax, ymin=fc.ymin)
kpPlotMarkers(kp, top.genes, labels = ordered_symbols, text.orientation = "vertical", 
              r0=points.top, label.dist = 0.008, 
              label.color="#444444", 
              line.color = "#777777",
              cex=0.6)

## Data panel 2
kp <- kpPlotDensity(kp, data=hg.genes, window.size = 10e4, data.panel = 2)

dev.off()

# === Setup and Plot Parameters for Chr 3 plot ===
# Ensure plotting parameters are set

# === Setup and Plot Parameters ===
# Ensure plotting parameters are set
pp <- getDefaultPlotParams(plot.type = 2)
pp$data2height <- 100
pp$ideogramheight <- 100

# Specify the output file for the PDF
pdf(file.path(figures_dir, "Karyotype_obj.DEG.Chr3_only_top50_sig_genes.pdf"), width = 12, height = 8, paper = 'special')

# === Plot Karyotype for Chromosome 3 ===
# Adjust genome to "hg38" or "hg19" as needed, ensure all data is consistent with this genome
kp <- plotKaryotype(plot.type = 2, plot.params = pp, chromosomes = c("chr3"), genome = "hg38")

# Add main title to the plot
kpAddMainTitle(kp, main = "Gene expression - Chr3 loss True vs False")

# === Plot Data Panel 1 ===
# Handle adj-p values that are 0 by replacing them with a very small number
sign.genes$p_val_adj[sign.genes$p_val_adj == 0] <- 1e-300

# Derive bubble size (cex) based on adj-p values (smaller adj-p values = larger bubbles)
# Normalize adj-p values to a suitable range for bubble sizes (e.g., between 0.5 and 2)
cex.val <- -log10(sign.genes$p_val_adj)  # Use negative log10 of adj-p values to reflect significance
cex.val <- scales::rescale(cex.val, to = c(0.5, 2))  # Rescale bubble sizes between 0.5 and 2

# Plot gene expression points where y = log2FC and bubble size = adj-p values
kpPoints(kp, data = sign.genes, y = sign.genes$avg_log2FC, cex = cex.val, ymax = fc.ymax, ymin = fc.ymin, r1 = points.top, col = sign.col)

# Plot gene segments
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes)) / 2
kpSegments(kp, 
           chr = as.character(seqnames(top.genes)), 
           x0 = gene.mean, 
           x1 = gene.mean, 
           y0 = top.genes$avg_log2FC,  
           y1 = fc.ymax,  
           ymax = fc.ymax, 
           ymin = fc.ymin, 
           r1 = points.top, 
           col = "#777777")

# Add axis and labels (without ymax and ymin in kpAddLabels)
kpAxis(kp, ymax = fc.ymax, ymin = fc.ymin)
kpAddLabels(kp, labels = "avg_log2FC", srt = 90, pos = 1, label.margin = 0.06)

# Plot gene markers
kpPlotMarkers(kp, top.genes, labels = ordered_symbols, text.orientation = "vertical", 
              r0 = points.top, label.dist = 0.008, 
              label.color = "#444444", 
              line.color = "#777777", cex = 0.6)

# === Mark the Cas9 Cut Site Region (chr3:79,901,054 - chr3:79,901,074) ===
# Draw a short vertical line at the Cas9 cut site on top of the Chr3 location (adjust the height using r0 and r1)
kpSegments(kp, 
           chr = "chr3", 
           x0 = 79901054, 
           x1 = 79901074, 
           y0 = 0.8, y1 = 1,   # Adjust r0 and r1 to make the line short and on top of the ideogram
           col = "red", lwd = 2, r0 = 0.8, r1 = 1)  # Adjusting the top part of the plot

# Highlight genes beyond the Cas9 cut site
affected_genes <- top.genes[start(top.genes) > 79901054]
kpPoints(kp, data = affected_genes, y = affected_genes$avg_log2FC, cex = cex.val, ymax = fc.ymax, ymin = fc.ymin, r1 = points.top, col = "blue")

# === Plot Data Panel 2 (Density Plot) ===
# Ensure hg.genes is compatible with genome assembly (hg38 in this case)
kp <- kpPlotDensity(kp, data = hg.genes, window.size = 10e4, data.panel = 2)

# === Add Figure Legend (with real values for adj-p bubble size) ===
# Define bubble sizes for the legend corresponding to actual adj-p values
legend_sizes <- c(0.5, 1, 1.5, 2)  # Bubble sizes (for example, corresponding to different p-value ranges)
legend_values <- c(0.05, 0.01, 0.001, 0.0001)  # Corresponding adjusted p-values (lower p-values = larger bubbles)

# Manually plot the legend showing the relation between bubble size and adj-p values
legend("topright", 
       legend = paste(legend_values, "adj-p"),  # Display actual adj-p values
       title = "Bubble Size (adj-p)", 
       pch = 16,                                 # Use filled circles for the bubbles
       pt.cex = legend_sizes,                    # Set the sizes of the points (bubble sizes)
       col = "black",                            # Color of the bubbles in the legend
       bty = "n",                                # No box around the legend
       cex = 0.8)                                # Size of the text in the legend

# Close the PDF output
dev.off()

# ========= plot for Chr Y =====

# Specify the output file for the PDF
pdf(file.path(figures_dir, "Karyotype_obj.DEG.ChrY_only_top50_sig_genes.pdf"), width = 12, height = 8, paper = 'special')

# === Plot Karyotype for Chromosome 3 ===
# Adjust genome to "hg38" or "hg19" as needed, ensure all data is consistent with this genome
kp <- plotKaryotype(plot.type = 2, plot.params = pp, chromosomes = c("chrY"), genome = "hg38")

# Add main title to the plot
kpAddMainTitle(kp, main = "Gene expression - ChrY loss True vs False")

# === Plot Data Panel 1 ===
# Handle adj-p values that are 0 by replacing them with a very small number
sign.genes$p_val_adj[sign.genes$p_val_adj == 0] <- 1e-300

# Derive bubble size (cex) based on adj-p values (smaller adj-p values = larger bubbles)
# Normalize adj-p values to a suitable range for bubble sizes (e.g., between 0.5 and 2)
cex.val <- -log10(sign.genes$p_val_adj)  # Use negative log10 of adj-p values to reflect significance
cex.val <- scales::rescale(cex.val, to = c(0.5, 2))  # Rescale bubble sizes between 0.5 and 2

# Plot gene expression points where y = log2FC and bubble size = adj-p values
kpPoints(kp, data = sign.genes, y = sign.genes$avg_log2FC, cex = cex.val, ymax = fc.ymax, ymin = fc.ymin, r1 = points.top, col = sign.col)

# Plot gene segments
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes)) / 2
kpSegments(kp, 
           chr = as.character(seqnames(top.genes)), 
           x0 = gene.mean, 
           x1 = gene.mean, 
           y0 = top.genes$avg_log2FC,  
           y1 = fc.ymax,  
           ymax = fc.ymax, 
           ymin = fc.ymin, 
           r1 = points.top, 
           col = "#777777")

# Add axis and labels (without ymax and ymin in kpAddLabels)
kpAxis(kp, ymax = fc.ymax, ymin = fc.ymin)
kpAddLabels(kp, labels = "avg_log2FC", srt = 90, pos = 1, label.margin = 0.06)

# Plot gene markers
kpPlotMarkers(kp, top.genes, labels = ordered_symbols, text.orientation = "vertical", 
              r0 = points.top, label.dist = 0.008, 
              label.color = "#444444", 
              line.color = "#777777", cex = 0.6)

# === Mark the Cas9 Cut Site Region (chrY:13,042,797-chrY:13,042,817) ===
# Draw a short vertical line at the Cas9 cut site on top of the ChrY location (adjust the height using r0 and r1)
kpSegments(kp, 
           chr = "chrY", 
           x0 = 13042797, 
           x1 = 13042817, 
           y0 = 0.8, y1 = 1,   # Adjust r0 and r1 to make the line short and on top of the ideogram
           col = "red", lwd = 2, r0 = 0.8, r1 = 1)  # Adjusting the top part of the plot

# Highlight genes beyond the Cas9 cut site
affected_genes <- top.genes[start(top.genes) > 79901054]
kpPoints(kp, data = affected_genes, y = affected_genes$avg_log2FC, cex = cex.val, ymax = fc.ymax, ymin = fc.ymin, r1 = points.top, col = "blue")

# === Plot Data Panel 2 (Density Plot) ===
# Ensure hg.genes is compatible with genome assembly (hg38 in this case)
kp <- kpPlotDensity(kp, data = hg.genes, window.size = 10e4, data.panel = 2)

# === Add Figure Legend (with real values for adj-p bubble size) ===
# Define bubble sizes for the legend corresponding to actual adj-p values
legend_sizes <- c(0.5, 1, 1.5, 2)  # Bubble sizes (for example, corresponding to different p-value ranges)
legend_values <- c(0.05, 0.01, 0.001, 0.0001)  # Corresponding adjusted p-values (lower p-values = larger bubbles)

# Manually plot the legend showing the relation between bubble size and adj-p values
legend("topright", 
       legend = paste(legend_values, "adj-p"),  # Display actual adj-p values
       title = "Bubble Size (adj-p)", 
       pch = 16,                                 # Use filled circles for the bubbles
       pt.cex = legend_sizes,                    # Set the sizes of the points (bubble sizes)
       col = "black",                            # Color of the bubbles in the legend
       bty = "n",                                # No box around the legend
       cex = 0.8)                                # Size of the text in the legend

# Close the PDF output
dev.off()


# save current working env
save.image("./DATA/Integrated/KaryoploteR/03_3p_loss.RData")
# To restore work environment later use: load("./DATA/Integrated/KaryoploteR/03_3p_loss.RData") 

