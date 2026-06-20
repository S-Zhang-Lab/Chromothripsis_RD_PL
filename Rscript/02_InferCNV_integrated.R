##========================================##
##             inferCVN                   ##
##========================================##
#Load Packages
library(dplyr)
library(Seurat)
library(ggplot2)
library(sctransform)
library(infercnv)
library(patchwork)
library(karyoploteR)
library(tidyverse)
library(future)
library(Matrix)

plan()
plan(multicore, workers = 22) # make the computer has enough cores. Use carefully. 
plan()
options(future.globals.maxSize = 50 * 1024^3)  # 50 GiB
options(scipen = 100) # set to avoid using scientific notation. 

base_dir <- "/path/to/Chromothripsis_2"  # <-- set to your local path to the project root
figures_dir <- file.path(base_dir, "Figures/integrated/inferCNV")

# load data
obj <- readRDS(file=file.path(base_dir, "/DATA/Integrated/integrated_obj.rds"))

# Set up counts matrix
# counts_matrix <- as.data.frame(obj@assays$SCT$counts)
counts_matrix <- obj@assays$SCT$counts
meta <- obj@meta.data

annot <- as.data.frame(obj@meta.data$HTO_classification)
colnames(annot) <- c('V1')
annot$V1 <- as.character(annot$V1)
rownames(annot) <- as.character(colnames(counts_matrix))
annot$cell <- recode(
  annot$V1, 
  "L1-H-HTO-5" = "sgNTC-D5", 
  "L1-H-HTO-6" = "sgNTC-D5", 
  "L1-H-HTO-2" = "sgYq-D5", 
  "L1-H-HTO-3" = "sgYq-D5",
  "L2-H-HTO-2" = "sg3p-D5",
  "L2-H-HTO-3" = "sg3p-D5",
  "L3-H-HTO-1" = "sgNTC-D10", 
  "L3-H-HTO-6" = "sgNTC-D10", 
  "L3-H-HTO-4" = "sgYq-D10", 
  "L3-H-HTO-5" = "sgYq-D10",
  "L3-H-HTO-2" = "sg3p-D10",
  "L3-H-HTO-3" = "sg3p-D10"
)
annot$V1 <- NULL

# Load gene order file
gene_order <- read.table('./Misc/gene_order_GRCh28_2020_V2.txt', header = FALSE)

# Create unique row names using the first column (gene symbols)
gene_order$V1 <- make.unique(as.character(gene_order$V1))

# Set the first column as row names
rownames(gene_order) <- gene_order$V1

# Optionally, remove the first column (if it's now used as row names)
gene_order <- gene_order[ , -1]

# Filter for genes on ChrX and ChrY in the gene_order file
genes_on_chrX <- rownames(gene_order[gene_order$V2 == "chrX", ])
genes_on_chrY <- rownames(gene_order[gene_order$V2 == "chrY", ])

# Combine the two lists
genes_on_sex_chromosomes <- c(genes_on_chrX, genes_on_chrY)

# Display the number of genes found on ChrX and ChrY
length(genes_on_chrX)  # Number of genes on ChrX
length(genes_on_chrY)  # Number of genes on ChrY

# Check for matching genes on ChrX and ChrY in counts_matrix
matched_genes_on_sex_chromosomes <- rownames(counts_matrix)[rownames(counts_matrix) %in% genes_on_sex_chromosomes]

# Display the matched genes
matched_genes_on_sex_chromosomes

# Look up the chromosome information for the matched genes
gene_order[matched_genes_on_sex_chromosomes, ]

# Check the expression levels of the matched genes in the counts matrix
expression_data_chrX <- counts_matrix[matched_genes_on_sex_chromosomes, ]

# View the expression data for these genes
expression_data_chrX[1:10,1:10]

# Before creating the infercnv_obj, check the number of ChrY genes
genes_on_chrY <- rownames(gene_order[gene_order$V2 == "chrY", ])

# Check if these genes exist in the counts_matrix before creating the object
chrY_genes_in_matrix <- rownames(counts_matrix)[rownames(counts_matrix) %in% genes_on_chrY]
length(chrY_genes_in_matrix)  # Number of ChrY genes in the counts_matrix

# Before creating the infercnv_obj, check the number of ChrY genes
genes_on_chrX <- rownames(gene_order[gene_order$V2 == "chrX", ])

# Check if these genes exist in the counts_matrix before creating the object
chrX_genes_in_matrix <- rownames(counts_matrix)[rownames(counts_matrix) %in% genes_on_chrX]
length(chrX_genes_in_matrix)  # Number of ChrX genes in the counts_matrix

# Create inferCNV group

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = counts_matrix,
                                     annotations_file = annot,
                                     gene_order_file = gene_order,
                                     chr_exclude = NULL,
                                     ref_group_names = 'sgNTC-D5')


# Run inferCNV
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=file.path(base_dir, "DATA/Integrated/InferCNV/CNV_output"),
                             output_format="pdf",
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=FALSE,
                             num_threads=22,
                             no_plot=FALSE)


# HMM mode
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=file.path(base_dir, "/DATA/Integrated/InferCNV/CNV_HMM_output"),
                             output_format="pdf",
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE,
                             num_threads=22,
)



# write CNV results
cnvRes <- infercnv_obj@expr.data
cnvRes[1:4,1:4]


# write CNV back to Seurat meta data
obj = infercnv::add_to_seurat(
  infercnv_output_path=file.path(base_dir, "/DATA/Integrated/InferCNV/CNV_HMM_output"),
  seurat_obj=obj, # optional
  top_n=50
)

meta <- obj@meta.data
colnames(meta)

DimPlot(obj, group.by = "HTO_classification")
FeaturePlot(obj , features="has_loss_chr3") + 
  ggplot2::scale_colour_gradient(low="lightgrey", high="blue", limits=c(0,1))

FeaturePlot(obj , features="has_dupli_chr3") + 
  ggplot2::scale_colour_gradient(low="lightgrey", high="blue", limits=c(0,1))

FeaturePlot(obj , features="top_loss_3") + 
  ggplot2::scale_colour_gradient(low="lightgrey", high="blue", limits=c(0,1))

FeaturePlot(obj , features="proportion_scaled_cnv_chr3") + 
  ggplot2::scale_colour_gradient(low="lightgrey", high="blue", limits=c(0,1))

FeaturePlot(obj , features="has_loss_chrY") + 
  ggplot2::scale_colour_gradient(low="lightgrey", high="blue", limits=c(0,1))



save.image(file = "obj.RData")
saveRDS(obj, file = file.path(base_dir, "/DATA/Integrated/InferCNV/Post_inferCNV_obj.rds"))



