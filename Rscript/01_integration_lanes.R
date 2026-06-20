# ==== integration of L1,L2,L3 ======= #
#Load Packages
library(dplyr)
library(Seurat)
library(ggplot2)
library(sctransform)
library(infercnv)
library(patchwork)
library(karyoploteR)
library(tidyverse)

# Parallelization in Seurat
library(future)
plan()
plan(multicore, workers = 22) # make the computer has enough cores. Use carefully. 
plan()
options(future.globals.maxSize = 50 * 1024^3)  # 50 GiB

# Integration of L1-3 from cellranger output.  Reference: https://satijalab.org/seurat/articles/integration_introduction
# load Cellranger 7.2 output

base_dir <- "/path/to/Chromothripsis_2"  # <-- set to your local path to the project root
figures_dir <- file.path(base_dir, "Figures/integrated")

# load 10X data (NOTE: "filtered_feature_bc_matrix" is under 
# "outs" folder from Cellranger V7.2)
L1.data <- Read10X(
  data.dir = file.path(base_dir, "/DATA/UTSW20_CITE-seq_L1/filtered_feature_bc_matrix"))
L2.data <- Read10X(
  data.dir = file.path(base_dir, "/DATA/UTSW20_CITE-seq_L2/filtered_feature_bc_matrix"))
L3.data <- Read10X(
  data.dir = file.path(base_dir, "/DATA/UTSW20_CITE-seq_L3/filtered_feature_bc_matrix"))

L1.data

# Initialize the Seurat object with the raw (non-normalized data).
L1 <- CreateSeuratObject(counts = L1.data$`Gene Expression`, project = "L1")

# Modify HTO names for L1
hto_L1 <- L1.data$`Antibody Capture`
rownames(hto_L1) <- paste0("L1_", rownames(hto_L1))

# Create CITE assay for L1 with modified HTO names
L1[["Abs"]] <- CreateAssayObject(counts = hto_L1)

L1.Abs <- t(as.data.frame(L1$"Abs"@counts))

# Subset HTO groups for cell group demultiplexing
L1.HTO <- t(L1.Abs[,c(3,4,6,7)])

# Create HTO assay with only R2 reads
L1[["HTO"]] <- CreateAssayObject(counts = L1.HTO) 

rm(L1.Abs)

# View the assays presented
L1 

L1 <- NormalizeData(L1, assay = "HTO", normalization.method = "CLR")
L1 <- HTODemux(L1, assay = "HTO", positive.quantile = 0.99)

# Visualize demultiplexing results
table(L1$HTO_classification.global)

# Visualize enrichment for selected HTOs with ridge plots
# Group cells based on the max HTO signal
Idents(L1) <- "HTO_maxID"

pdf(file.path(figures_dir, "Ridgeplot_HTO.pdf"),width=8,height=12,paper='special')
RidgePlot(L1, assay = "HTO", features = rownames(L1[["HTO"]])[1:6], ncol = 2)
dev.off()

# pull metadata and examine the distributions table
meta.L1 <- L1@meta.data
table(meta.L1$HTO_maxID)
table(meta.L1$hash.ID)
table(meta.L1$HTO_classification)
table(meta.L1$HTO_classification.global)

# change Idents 
Idents(L1) = "hash.ID"

# Get the cell barcodes that are NOT "Negative" or "Doublet"
valid_barcodes <- rownames(meta.L1)[meta.L1$hash.ID != "Negative" & meta.L1$hash.ID != "Doublet"]

# Subset your Seurat object
L1 <- subset(L1, cells = valid_barcodes)

# L2 
# Initialize the Seurat object with the raw (non-normalized data).
L2 <- CreateSeuratObject(counts = L2.data$`Gene Expression`, project = "L2")

# Modify HTO names for L2
hto_L2 <- L2.data$`Antibody Capture`
rownames(hto_L2) <- paste0("L2_", rownames(hto_L2))

# Create CITE assay for L2 with modified HTO names
L2[["Abs"]] <- CreateAssayObject(counts = hto_L2)

L2.Abs <- t(as.data.frame(L2$"Abs"@counts))

# Subset HTO groups for cell group demultiplexing
L2.HTO <- t(L2.Abs[,c(3,4)])

# Create HTO assay with only R2 reads
L2[["HTO"]] <- CreateAssayObject(counts = L2.HTO) 

rm(L2.Abs)

# View the assays presented
L2 

L2 <- NormalizeData(L2, assay = "HTO", normalization.method = "CLR")
L2 <- HTODemux(L2, assay = "HTO", positive.quantile = 0.99)

# Visualize demultiplexing results
table(L2$HTO_classification.global)

# Visualize enrichment for selected HTOs with ridge plots
# Group cells based on the max HTO signal
Idents(L2) <- "HTO_maxID"

pdf(file.path(figures_dir, "Ridgeplot_L2_HTO.pdf"),width=8,height=12,paper='special')
RidgePlot(L2, assay = "HTO", features = rownames(L2[["HTO"]])[1:6], ncol = 2)
dev.off()

# pull metadata and examine the distributions table
meta.L2 <- L2@meta.data
table(meta.L2$HTO_maxID)
table(meta.L2$hash.ID)
table(meta.L2$HTO_classification)
table(meta.L2$HTO_classification.global)

# change Idents 
Idents(L2) = "hash.ID"

# Get the cell barcodes that are NOT "Negative" or "Doublet"
valid_barcodes <- rownames(meta.L2)[meta.L2$hash.ID != "Negative" & meta.L2$hash.ID != "Doublet"]

# Subset your Seurat object
L2 <- subset(L2, cells = valid_barcodes)

# L3
# Initialize the Seurat object with the raw (non-normalized data).
L3 <- CreateSeuratObject(counts = L3.data$`Gene Expression`, project = "L3")

# Modify HTO names for L3
hto_L3 <- L3.data$`Antibody Capture`
rownames(hto_L3) <- paste0("L3_", rownames(hto_L3))

# Create CITE assay for L3 with modified HTO names
L3[["Abs"]] <- CreateAssayObject(counts = hto_L3)

L3.Abs <- t(as.data.frame(L3$"Abs"@counts))

# Subset HTO groups for cell group demultiplexing
L3.HTO <- t(L3.Abs[,c(1,3:7)])
rownames(L3.HTO)
rownames(L3.HTO)[rownames(L3.HTO) == "L3-H-HTO-1.1"] <- "L3-H-HTO-1"

# Create HTO assay with only R2 reads
L3[["HTO"]] <- CreateAssayObject(counts = L3.HTO) 

rm(L3.Abs)

# View the assays presented
L3 

L3 <- NormalizeData(L3, assay = "HTO", normalization.method = "CLR")
L3 <- HTODemux(L3, assay = "HTO", positive.quantile = 0.99)

# Visualize demultiplexing results
table(L3$HTO_classification.global)

# Visualize enrichment for selected HTOs with ridge plots
# Group cells based on the max HTO signal
Idents(L3) <- "HTO_maxID"

pdf(file.path(figures_dir, "L3_Ridgeplot_HTO.pdf"),width=8,height=12,paper='special')
RidgePlot(L3, assay = "HTO", features = rownames(L3[["HTO"]])[1:6], ncol = 2)
dev.off()

# pull metadata and examine the distributions table
meta.L3 <- L3@meta.data
table(meta.L3$HTO_maxID)
table(meta.L3$hash.ID)
table(meta.L3$HTO_classification)
table(meta.L3$HTO_classification.global)

# change Idents 
Idents(L3) = "hash.ID"

# Get the cell barcodes that are NOT "Negative" or "Doublet"
valid_barcodes <- rownames(meta.L3)[meta.L3$hash.ID != "Negative" & meta.L3$hash.ID != "Doublet"]

# Subset your Seurat object
L3 <- subset(L3, cells = valid_barcodes)



# ===== merge different lanes to get a merged single seurat object ====
obj <- merge(x = L1, y = list(L2, L3)) 
Lane.backup <- obj

# remove the other raw objects to reduce project size

rm(L1.data)
rm(L2.data)
rm(L3.data)
rm(L1)
rm(L2)
rm(L3)


##### Test Objects integration in Seurat5 https://satijalab.org/seurat/articles/seurat5_integration ######### 
DefaultAssay(obj)
set.seed(123)

obj <- SCTransform(obj)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(
  object = obj,
  method = RPCAIntegration,
  normalization.method = "SCT",
  verbose = F
)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "integrated.dr")
obj <- FindClusters(obj, resolution = 2)
obj <- RunUMAP(obj, dims = 1:30, reduction = "integrated.dr")

# pull metadata
meta <- obj@meta.data

DimPlot(
  obj,
  reduction = "integrated.dr",
  group.by = c("orig.ident"),
  combine = FALSE, label.size = 2
)

DimPlot(
  obj,
  reduction = "integrated.dr",
  group.by = c("orig.ident"),
  combine = FALSE, label.size = 2, 
  split.by = "HTO_classification"
)


saveRDS(obj, file=file.path(base_dir, "/DATA/Integrated/integrated_obj.rds"))



sessionInfo()