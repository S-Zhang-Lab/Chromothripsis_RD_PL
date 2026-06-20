#Load Packages
library(dplyr)
library(Seurat)
library(ggplot2)
library(sctransform)
library(infercnv)
library(patchwork)
library(karyoploteR)

base_dir <- "/path/to/Chromothripsis_2"  # <-- set to your local path to the project root
figures_dir <- file.path(base_dir, "Figures/L1")

# load 10X data (NOTE: "filtered_feature_bc_matrix" is under 
# "outs" folder from Cellranger V7.2)
L1.data <- Read10X(
  data.dir = file.path(base_dir, "/DATA/UTSW20_CITE-seq_L1/filtered_feature_bc_matrix"))


# Examine the structure of the data. It should be a list of 2 objects
str(L1.data) 

# Initialize the Seurat object with the raw (non-normalized data).
L1 <- CreateSeuratObject(counts = L1.data$"Gene Expression", 
                             project = "L1")

# Add HTO objects
L1$Abs <- CreateAssayObject(counts = L1.data$"Antibody Capture")
L1.Abs <- t(as.data.frame(L1$"Abs"@counts))

# Subset HTO groups for cell group demultiplexing
L1.Abs <- t(L1.Abs[,c(3,4,6,7)])

# Create HTO assay with only R2 reads
L1[["HTO"]] <- CreateAssayObject(counts = L1.Abs) 
rm(L1.Abs)

# View the assays presented
L1 

# The [[ operator can add columns to object metadata. 
L1[["percent.mt"]] <- PercentageFeatureSet(L1, pattern = "^MT-")
meta.L1 <- L1@meta.data

# Visualize QC metrics as a violin plot
pdf(file.path(figures_dir, "QC.pdf"),width=6,height=10)
VlnPlot(L1, features = c(
  "nFeature_RNA", "nCount_RNA", "nCount_HTO","percent.mt"), ncol = 2)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, 
# but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf(file.path(figures_dir, "Featureplot.pdf"),width=10,height=4,paper='special')
plot1 <- FeatureScatter(L1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(L1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

pdf(file.path(figures_dir, "Featureplot_HTO.pdf"),width=10,height=4,paper='special')
plot3 <- FeatureScatter(L1, feature1 = "nFeature_HTO", feature2 = "nCount_HTO")
plot4 <- FeatureScatter(L1, feature1 = "nFeature_HTO", feature2 = "nCount_HTO")
CombinePlots(plots = list(plot3, plot4))
dev.off()

# Removing unwanted cells from the data based on above QC plots
L1 <- subset(L1, subset = nFeature_RNA > 200 & nFeature_RNA < 15000 & percent.mt < 25)


###################################################
#           HTO demultiplexing                    #
###################################################

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

###############################################
#       Run normalization and downstream      #
###############################################

options(future.globals.maxSize = 60 * 1024^3)  # set to ram use to 60 GB

# run sctransform on RNA assay
L1 <- SCTransform(L1, vars.to.regress = "percent.mt", verbose = TRUE)

# finding variable genes
L1 <- FindVariableFeatures(L1, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(L1), 10)

# plot variable features with and without labels
pdf(file.path(figures_dir, "Variable.Genes.pdf"),width=10,height=4,paper='special')
plot1 <- VariableFeaturePlot(L1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()


# Perform linear dimensional reduction
set.seed(1)
L1 <- RunPCA(L1, verbose = FALSE)
pdf(file.path(figures_dir, "Elbowplot.pdf"),width=7,height=5,paper='special')
ElbowPlot(L1)
dev.off()

print(L1[["pca"]], dims = 1:5, nfeatures = 5)


pdf(file.path(figures_dir, "Dim.plot.pdf"),width=5,height=4,paper='special')
DimPlot(L1, reduction = "pca")
dev.off()

# Cluster the cells

L1 <- FindNeighbors(L1, dims = 1:10)
L1 <- FindClusters(L1, resolution = 1)
head(Idents(L1), 5)

# Run non-linear dimensional reduction (UMAP/tSNE)
L1 <- RunUMAP(L1, dims = 1:10, verbose = TRUE)

# note that you can set `label = TRUE` or use the LabelClusters function 
# to help label individual clusters
pdf(file.path(figures_dir, "DimPlot.UMAP.pdf"),width=15,height=4,paper='special')
DimPlot(L1, reduction = "umap", label = TRUE, split.by = "hash.ID")
dev.off()

# find markers for every cluster compared to all remaining cells, report only the positive ones
L1.markers <- FindAllMarkers(
  L1, only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25)
L1.markers %>% 
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# do heatmap
top10 <- L1.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

pdf(file.path(figures_dir, "Heatmap.top10.pdf"),width=20,height=20,paper='special')
DoHeatmap(L1, features = top10$gene) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) 
dev.off()

head(L1@meta.data)
meta.L1 <- L1@meta.data

# stacked bar chart
pdf(file.path(figures_dir, "Barchart_Freq_Stack_fill.pdf"),width=8,height=5,paper='special')
ggplot(L1@meta.data, aes(x=hash.ID, fill=seurat_clusters)) + 
  geom_bar(position = "fill")
dev.off()

pdf(file.path(figures_dir, "Barchart_Freq_Stack_hash.ID.pdf"),width=8,height=5,paper='special')
ggplot(L1@meta.data, aes(x = hash.ID, fill=seurat_clusters)) +  
  geom_bar(aes(y = (..count..)/sum(..count..)))
dev.off()

pdf(file.path(figures_dir, "Barchart_Freq_byhash.ID.pdf"),width=10,height=5,paper='special')
ggplot(L1@meta.data, aes(seurat_clusters, group = hash.ID)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") + 
  scale_y_continuous(labels=scales::percent) +
  ylab("relative frequencies") +
  facet_grid(~hash.ID)
dev.off()

# plot marker genes for cluster
pdf(file.path(figures_dir, "Featureplot_C0_IQCD.pdf"),width=5,height=5,paper='special')
FeaturePlot(L1, features = "IQCD")
dev.off()

pdf(file.path(figures_dir, "Featureplot_TOP2A.pdf"),width=5,height=5,paper='special')
FeaturePlot(L1, features = "TOP2A")
dev.off()

DefaultAssay(L1) <- "SCT"
pdf(file.path(figures_dir, "Featureplot_VHL_split.pdf"),width=6,height=6,paper='special')
FeaturePlot(L1, features = "VHL") + 
  ggplot2::scale_colour_gradient(low="lightgrey", high="red", limits=c(0,1))
dev.off()


##========================================##
##             inferCVN                   ##
##========================================##

library(infercnv)

# Set up counts matrix
counts_matrix <- as.data.frame(L1@assays$RNA$counts)
meta <- L1@meta.data

# Annotate immune cells based on SingleR output
annot <- as.data.frame(L1@meta.data$hash.ID)
colnames(annot) <- c('V1')
annot$V1 <- as.character(annot$V1)
rownames(annot) <- as.character(colnames(counts_matrix))
control <- c("NTC")
annot$cell <- recode(
  annot$V1, 
  "H-HTO-5" = "NTC", 
  "H-HTO-6" = "NTC", 
  "H-HTO-2" = "Yq", 
  "H-HTO-3" = "Yq"
  )
annot$V1 <- NULL

# Load gene order file
gene_order <- read.table('./Misc/hg38_gencode_v27.txt', header = F, row.names = 1)

# Create inferCNV object
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = counts_matrix,
                                     annotations_file = annot,
                                     gene_order_file = gene_order,
                                     ref_group_names = 'NTC')

# Run inferCNV
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=file.path(base_dir, "/DATA/InfCNV/L1/CNV_output"),
                             output_format="pdf",
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=FALSE,
                             num_threads=22,
                             no_plot=FALSE)


# HMM mode
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=file.path(base_dir, "/DATA/InfCNV/L1/CNV_HMM_output"),
                             output_format="pdf",
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE,
                             num_threads=22,
                             no_plot=FALSE)



# write CNV results
cnvRes <- infercnv_obj@expr.data
write.table(
  cnvRes,
  file = ('cnvRes.csv'),
  sep = '\t',
  row.names = T,
  col.names = T,
  quote = F)

cnvRes[1:4,1:4]

# write CNV back to Seurat meta data
L1 = infercnv::add_to_seurat(
  infercnv_output_path=file.path(base_dir, "/DATA/InfCNV/L1/CNV_HMM_output"),
  seurat_obj=L1, # optional
  top_n=20
  )

meta <- L1@meta.data
colnames(meta)

FeaturePlot(L1 , features="has_loss_chr3") + 
  ggplot2::scale_colour_gradient(low="lightgrey", high="blue", limits=c(0,1))
FeaturePlot(L1 , features="has_dupli_chr3") + 
  ggplot2::scale_colour_gradient(low="lightgrey", high="blue", limits=c(0,1))


## extract annnotation
HMM_anno <- colnames(L1@meta.data[, 19:257])
Chr3_changes <- c(colnames(L1@meta.data[, 37:45]))
Chr3_changes               
Chr5_changes <- c(colnames(L1@meta.data[, 55:63]))
Chr5_changes               

                  
pdf(file.path(figures_dir, "Featureplot_CNV_3p.pdf"),
    width=12,height=10,paper='special')
FeaturePlot(L1, features = Chr3_changes)
dev.off()

pdf(file.path(figures_dir, "Featureplot_CNV_5q.pdf"),
    width=12,height=10,paper='special')
FeaturePlot(L1, features = Chr5_changes)
dev.off()

pdf(file.path(figures_dir, "Featureplot_CNV_3p_5q.pdf"),
    width=16,height=16,paper='special')
FeaturePlot(L1, features = c(Chr3_changes, Chr5_changes))
dev.off()

######################
# Test MULTIseqDemux #
######################
# NOTE: Did not performed better than Seurat default Demux method #############
# Initialize the Seurat object with the raw (non-normalized data).
df <- CreateSeuratObject(counts = L1.data$"Gene Expression", project = "MULTI-Chr")

# add ADT/HTO/BC objects
df$HTO <- CreateAssayObject(counts = L1.data$"Antibody Capture")

L1.HTOs <- t(as.data.frame(df$"HTO"@counts))
L1.HTOs <- t(L1.HTOs[,7:12]) #subset only R2 reads for cell group demultiplexing
df[["HTO"]] <- CreateAssayObject(counts = (L1.HTOs+0.1)*10) # create HTO assay with only R2 reads counts 
df[["MULTI"]] <- CreateAssayObject(counts = (L1.HTOs)) #create MULTI-seq assay based on raw HTO counts with only R2 reads
rm(L1.HTOs)

str(df)

# do MULTIseqDemux on L1 df HTO assay (contains both R1 and R2)
L1 <- MULTIseqDemux(
  L1,
  assay = "HTO",
  quantile = 0.9,
  autoThresh = TRUE,
  maxiter = 5,
  qrange = seq(from = 0.1, to = 0.99, by = 0.01),
  verbose = TRUE
)

meta <- L1@meta.data

table(L1$MULTI_ID)
table(L1$MULTI_classification)

# Re-code ID by combining R1 classification and R2 classification. 
L1@meta.data$MULTIdemuxR1R2 = recode(
  L1@meta.data$MULTI_ID, 
  "R1-BC-ID1" = "BC-ID1", 
  "R2-BC-ID1" = "BC-ID1", 
  "R1-BC-ID2" = "BC-ID2", 
  "R2-BC-ID2" = "BC-ID2",
  "R1-BC-ID3" = "BC-ID3", 
  "R2-BC-ID3" = "BC-ID3",
  "R1-BC-ID4" = "BC-ID4",
  "R2-BC-ID4" = "BC-ID4",
  "R1-BC-ID5" = "BC-ID5", 
  "R2-BC-ID5" = "BC-ID5",
  "R1-BC-ID6" = "BC-ID6",
  "R2-BC-ID6" = "BC-ID6"
)

table(L1$MULTIdemuxR1R2)


######################################################
#    select Chr3 cells highest scaled_cnv            # 
######################################################
Idents(L1) <- "seurat_clusters"
DimPlot(L1)

Idents(L1) <- "proportion_scaled_cnv_chr3"
DimPlot(L1)

# Select cells with the highest 3p loss with very stringent cut off (>0.15)
L1.3ploss = subset(L1, idents = "0", invert = TRUE)
L1.3ploss@meta.data$Chr3p.loss <- ifelse(
  L1.3ploss@meta.data$proportion_scaled_cnv_chr3 > 0.15, "TRUE", "FALSE"
  )
meta.3ploss = L1.3ploss@meta.data
table(meta.3ploss$Chr3p.loss)


Idents(L1.3ploss) <- "Chr3p.loss"
DimPlot(L1.3ploss, split.by = "HTO_maxID")

pdf(file.path(figures_dir, "DimPlot_3ploss_BC_ID.pdf"),
    width=16,height=4,paper='special')
DimPlot(L1.3ploss, split.by = "HTO_maxID")
dev.off()


table(meta.3ploss$HTO_classification)
table(meta.3ploss$HTO_maxID)
table(meta.3ploss$HTO_classification.global)
table(meta.3ploss$HTO_classification)

######################################################
#    select Chr5 cells highest scaled_cnv            # 
######################################################
Idents(L1) <- "seurat_clusters"
DimPlot(L1)

Idents(L1) <- "proportion_scaled_cnv_chr5"
DimPlot(L1)

# select cells with highest chr5 cnv with very stringent cut off (>0.2)
L1.5q.gain = subset(L1, idents = "0", invert = TRUE)
L1.5q.gain@meta.data$Chr5q.gain <- ifelse(
  L1.5q.gain@meta.data$proportion_scaled_cnv_chr5 > 0.2, "TRUE", "FALSE"
  )
meta.5q.gain = L1.5q.gain@meta.data
meta.5q.gain$Chr5q.gain

Idents(L1.5q.gain) <- "Chr5q.gain"
DimPlot(L1.5q.gain, split.by = "HTO_maxID")


pdf(file.path(figures_dir, "DimPlot_5q.gain_BC_ID.pdf"),
    width=16,height=4,paper='special')
DimPlot(L1.5q.gain, split.by = "HTO_maxID")
dev.off()


DefaultAssay(L1) <- "SCT"
DefaultAssay(L1.3ploss) <- "SCT"

# Examine Chr3p genes
library(readxl)
Chr_3p_genes_list <- read_excel(file.path(base_dir, "Genelist/Chr.3p_genes_list.xlsx"))
View(Chr_3p_genes_list)

Chr_3p_genes <- Chr_3p_genes_list$`Gene name`

pdf(file.path(figures_dir, "FeaturePlot_3ploss_markers.pdf"),width=12,height=12,paper='special')
FeaturePlot(L1.3ploss, features = c("GBE1", "ROBO1", "ROBO2", "ZNF717"))
dev.off()


# Examine Chr5q genes

Chr_5q_genes_list <- read_excel(file.path(base_dir, "Genelist/human_chr5_genes.xlsx"))
View(Chr5q)
Chr_5q_genes <- Chr_5q_genes_list$Gene_name
Chr5q15_21_gene <- c("CAST", "CHD1", "PCSK1")
DefaultAssay(L1.5q.gain) <- "SCT"
pdf(file.path(figures_dir, "FeaturePlot_5q.gain_markers.pdf"),width=12,height=12,paper='special')
FeaturePlot(L1.3ploss, features = Chr5q15_21_gene)
dev.off()


# ========= remove cells with high CNV from control group and reanalyze ========= #
# calculate mean of CNV of each cells and add back to meta.data. 
cnv_cols <- grep('proportion_scaled_cnv_chr', names(L1@meta.data), value = T)
cnvs <- L1@meta.data[, cnv_cols]
L1@meta.data$cnv_avg <- rowMeans(cnvs)

pdf(file.path(figures_dir, "cnv_avg_fullscale.pdf"),width=8,height=4,paper='special')
ggplot(L1@meta.data, aes(x = cnv_avg)) + 
  geom_histogram(colour = 2, fill = "white", binwidth = 0.0005) 
dev.off()

pdf(file.path(figures_dir, "cnv_avg_upto005.pdf"),width=8,height=4,paper='special')
ggplot(L1@meta.data, aes(x = cnv_avg)) + 
  geom_histogram(colour = 2, fill = "white", binwidth = 0.0005) + xlim(c(0, 0.05))
dev.off()

FeaturePlot(L1, features = "cnv_avg")
L1.clean <- subset(L1, subset = cnv_avg < 0.05)

L1.clean$cnv.class <- "cnv_avg.normal"
L1.clean$cnv.class[L1.clean$cnv_avg > 0.03] <- "cnv_avg.high"
L1.clean$cnv.class[L1.clean$cnv_avg > 0.01 & L1.clean$cnv_avg <= 0.03] <- "cnv_avg.mid"
L1.clean$cnv.class[L1.clean$cnv_avg > 0.00 & L1.clean$cnv_avg <= 0.01] <- "cnv_avg.low"

L1.clean$HTO_recode <- recode(
  L1.clean$HTO_maxID, 
  "R2-BC-ID1" = "01_ID1", 
  "R2-BC-ID2" = "02_ID2", 
  "R2-BC-ID3" = "03_ID3", 
  "R2-BC-ID4" = "04_ID4", 
  "R2-BC-ID5" = "05_ID5", 
  "R2-BC-ID6" = "06_ID6"
  )

Idents(L1.clean) <- "HTO_recode"
Idents(L1.clean) <- factor(x = Idents(L1.clean), levels = sort(levels(L1.clean))) # sort idents

Idents(L1.clean) <- "HTO_maxID"
Idents(L1.clean) <- factor(x = Idents(L1.clean), levels = sort(levels(L1.clean))) # sort idents

table(L1.clean$cnv.class)
pdf(file.path(figures_dir, "VlnPlot_L1.clean.VHL_cnv.class.pdf"),width=8,height=4,paper='special')
VlnPlot(L1.clean, feature = "VHL", split.by = "cnv.class")
dev.off()

#==  re select Chr3 cells highest scaled_cnv  ==# 

## Identify cells with highest chr3 changes
Idents(L1.clean) <- "seurat_clusters"
DimPlot(L1.clean)

Idents(L1.clean) <- "proportion_scaled_cnv_chr3"
DimPlot(L1.clean)

# Select cells with the highest 3p loss with very stringent cut off (>0.2)
L1.clean.3ploss = subset(L1.clean, idents = "0", invert = TRUE)
Idents(L1.clean.3ploss) <- "proportion_scaled_cnv_chr3"
DimPlot(L1.clean.3ploss)
L1.clean.3ploss@meta.data$Chr3p.loss <- ifelse(
  L1.clean.3ploss@meta.data$proportion_scaled_cnv_chr3 > 0.2, 
  "TRUE", "FALSE"
  )
meta.3ploss = L1.clean.3ploss@meta.data
table(meta.3ploss$Chr3p.loss)
Idents(L1.clean.3ploss) <- "Chr3p.loss"
DimPlot(L1.clean.3ploss)
DimPlot(L1.clean.3ploss, split.by = "HTO_maxID")

pdf(file.path(figures_dir, "DimPlot_L1.clean.3ploss_BC_ID.pdf"),
    width=16,height=4,paper='special')
DimPlot(L1.clean.3ploss, split.by = "HTO_maxID")
dev.off()


#==  re select Chr5 cells highest scaled_cnv  ==# 

Idents(L1.clean) <- "seurat_clusters"
DimPlot(L1.clean)

Idents(L1.clean) <- "proportion_scaled_cnv_chr5"
DimPlot(L1.clean)

# select cells with highest chr5 cnv with very stringent cut off (>0.2)
L1.clean.5q.gain = subset(L1.clean, idents = "0", invert = TRUE)
L1.clean.5q.gain@meta.data$Chr5q.gain <- ifelse(
  L1.clean.5q.gain@meta.data$proportion_scaled_cnv_chr5 > 0.2, 
  "TRUE", "FALSE"
  )
meta.5q.gain = L1.clean.5q.gain@meta.data
meta.5q.gain$Chr5q.gain

Idents(L1.clean.5q.gain) <- "Chr5q.gain"
DimPlot(L1.clean.5q.gain, split.by = "HTO_maxID")

pdf(file.path(figures_dir, "DimPlot_L1.clean.5q.gain_BC_ID.pdf"),
    width=16,height=4,paper='special')
DimPlot(L1.clean.5q.gain, split.by = "HTO_maxID")
dev.off()


Idents(L1.clean.5q.gain) <- "proportion_scaled_cnv_chr5"
pdf(file.path(figures_dir, "DimPlot_L1.clean.5q.gain_proportion_scaled_cnv_chr5.pdf"),
    width=16,height=4,paper='special')
DimPlot(L1.clean.5q.gain, split.by = "HTO_maxID")
dev.off()

# 5q cut site coordinates: 98,852,808. The cutting should happen between these two genes
FeaturePlot(L1.clean.5q.gain, feature = c("CAST", "CHD1")) 

Idents(L1.clean.5q.gain) <- "HTO_maxID"
Idents(L1.clean.5q.gain) <- factor(x = Idents(L1.clean.5q.gain), 
                                       levels = sort(levels(L1.clean.5q.gain))
                                       )
Idents(L1.clean.3ploss) <- "HTO_maxID"
Idents(L1.clean.3ploss) <- factor(x = Idents(L1.clean.3ploss), 
                                      levels = sort(levels(L1.clean.3ploss))
                                      )

pdf(file.path(figures_dir, "VlnPlot_L1.clean.5q.gain_CAST.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean.5q.gain, feature = "CAST", split.by = "proportion_scaled_cnv_chr5")
dev.off()

pdf(file.path(figures_dir, "VlnPlot_L1.clean.5q.gain_CHD1.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean.5q.gain, feature = "CHD1", split.by = "proportion_scaled_cnv_chr5")
dev.off()


VlnPlot(L1.clean.5q.gain, feature = "VHL", split.by = "proportion_scaled_cnv_chr5")

pdf(file.path(figures_dir, "VlnPlot_L1.clean.5q.gain_VHL_by.proportion_scaled_cnv_chr5.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "VHL", split.by = "proportion_scaled_cnv_chr5")
dev.off()

pdf(file.path(figures_dir, "VlnPlot_L1.clean.5q.gain_VHL_by.cnv.class.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "VHL", split.by = "cnv.class")
dev.off()

Idents(L1.clean) <- "HTO_maxID"

pdf(file.path(figures_dir, "VlnPlot_L1.clean.proportion_scaled_cnv_chr5_SQSTM1.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "SQSTM1", split.by = "proportion_scaled_cnv_chr5")
dev.off()

pdf(file.path(figures_dir, "VlnPlot_L1.clean.proportion_scaled_cnv_chr3_SQSTM1.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "SQSTM1", split.by = "proportion_scaled_cnv_chr3")
dev.off()

pdf(file.path(figures_dir, "VlnPlot_L1.clean.proportion_scaled_cnv_chr5_VHL.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "VHL", split.by = "proportion_scaled_cnv_chr5")
dev.off()

pdf(file.path(figures_dir, "VlnPlot_L1.clean.proportion_scaled_cnv_chr3_VHL.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "VHL", split.by = "proportion_scaled_cnv_chr3")
dev.off()

pdf(file.path(figures_dir, "VlnPlot_L1.clean.has_dupli_chr5_VHL.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "VHL", split.by = "has_dupli_chr5")
dev.off()

pdf(file.path(figures_dir, "VlnPlot_L1.clean.has_loss_chr3_VHL.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "VHL", split.by = "has_loss_chr3")
dev.off()

pdf(file.path(figures_dir, "VlnPlot_L1.clean.has_dupli_chr5_SQSTM1.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "SQSTM1", split.by = "has_dupli_chr5")
dev.off()

pdf(file.path(figures_dir, "VlnPlot_L1.clean.has_loss_chr3_SQSTM1.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "SQSTM1", split.by = "has_loss_chr3")
dev.off()

### good plots ====
pdf(file.path(figures_dir, "VlnPlot_L1.clean.proportion_scaled_dupli_chr5_SQSTM1.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "SQSTM1", split.by = "proportion_scaled_dupli_chr5")
dev.off()

pdf(file.path(figures_dir, "VlnPlot_L1.clean.proportion_scaled_loss_chr3_SQSTM1.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "SQSTM1", split.by = "proportion_scaled_loss_chr3")
dev.off()

pdf(file.path(figures_dir, "VlnPlot_L1.clean.proportion_scaled_loss_chr5_SQSTM1.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "SQSTM1", split.by = "proportion_scaled_loss_chr5")
dev.off()

pdf(file.path(figures_dir, "VlnPlot_L1.clean.proportion_scaled_dupli_chr3_SQSTM1.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "SQSTM1", split.by = "proportion_scaled_dupli_chr3")
dev.off()


VlnPlot(L1.clean, feature = "BAP1", split.by = "proportion_scaled_dupli_chr5")
VlnPlot(L1.clean, feature = "FANCD2", split.by = "proportion_scaled_dupli_chr5")


# ====  recode chr3.loss and chr5.dupli to simply the plotting

L1.clean$Chr3.loss <- ifelse(
  L1.clean@meta.data$proportion_scaled_loss_chr3 > 0.2, 
  "Chr3.loss_TRUE", "Chr3.loss_FALSE"
  )
L1.clean$Chr5.dupli <- ifelse(
  L1.clean@meta.data$proportion_scaled_dupli_chr5 > 0.2, 
  "Chr5.dupli_TRUE", "Chr5.dupli_FALSE"
  )

Idents(L1.clean) <- "HTO_maxID"
Idents(L1.clean) <- factor(x = Idents(L1.clean), levels = sort(levels(L1.clean)))

pdf(file.path(figures_dir, "VlnPlot_L1.clean.Chr3.loss_SQSTM1.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "SQSTM1", split.by = "Chr3.loss")
dev.off()

pdf(file.path(figures_dir, "VlnPlot_L1.clean.Chr5.dupli_SQSTM1.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "SQSTM1", split.by = "Chr5.dupli")
dev.off()

# ========== Global DE gene analysis between chr3 or ch5 CNV TRUE / FALSE group. =============
# set a volcano plot function 
# 
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

## Identify global DE genes by comparing cell with/without Chr3.loss
Idents(L1.clean) <- "Chr3.loss"
Idents(L1.clean) <- factor(x = Idents(L1.clean), levels = sort(levels(L1.clean)))

table(Idents(L1.clean))
DEG_Chr3.loss <- FindMarkers(L1.clean, ident.1 = c('Chr3.loss_TRUE'), ident.2 = c('Chr3.loss_FALSE'), min.pct = 0.25)
write.csv(DEG_Chr3.loss, file = file.path(figures_dir, "DEG_Chr3.loss_DEG.csv"))

pdf(file.path(figures_dir, "Volcano_L1.clean.DEG.Chr3.loss.pdf"),
    width=6,height=4,paper='special')
create_volcano_plot(DEG_Chr3.loss, "avg_log2FC", "p_val_adj", "DEG_Chr3.loss (right: increase in Chr.loss.TRUE)", 0.05)
dev.off()

Idents(L1.clean) <- "HTO_maxID"
Idents(L1.clean) <- factor(x = Idents(L1.clean), levels = sort(levels(L1.clean)))

pdf(file.path(figures_dir, "Volcano_L1.clean.DEG.Chr3.loss_INHBA.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "INHBA", split.by = "Chr3.loss")
dev.off()

pdf(file.path(figures_dir, "Volcano_L1.clean.DEG.Chr3.loss_SOD3.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "SOD3", split.by = "Chr3.loss")
dev.off()

## Identify global DE genes by comparing cell with/without Chr5.dupli
Idents(L1.clean) <- "Chr5.dupli"
Idents(L1.clean) <- factor(x = Idents(L1.clean), levels = sort(levels(L1.clean)))

table(Idents(L1.clean))
DEG_Chr5.dupli <- FindMarkers(L1.clean, ident.1 = c('Chr5.dupli_TRUE'), ident.2 = c('Chr5.dupli_FALSE'), min.pct = 0.25)
write.csv(DEG_Chr5.dupli, file = file.path(figures_dir, "DEG_Chr5.dupli_DEG.csv"))

# draw volcano plot
pdf(file.path(figures_dir, "Volcano_L1.clean.DEG.Chr5.dupli.pdf"),
    width=6,height=4,paper='special')
create_volcano_plot(DEG_Chr5.dupli, "avg_log2FC", "p_val_adj", "DEG_Chr5.dupli (right: increase in Chr5.dupli.TRUE)", 0.05)
dev.off()

Idents(L1.clean) <- "HTO_maxID"
Idents(L1.clean) <- factor(x = Idents(L1.clean), levels = sort(levels(L1.clean)))

pdf(file.path(figures_dir, "Volcano_L1.clean.DEG.Chr5.dupli_ANKHD1.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "ANKHD1", split.by = "Chr5.dupli")
dev.off()

pdf(file.path(figures_dir, "Volcano_L1.clean.DEG.Chr5.dupli_TAF7.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "TAF7", split.by = "Chr5.dupli")
dev.off()

pdf(file.path(figures_dir, "Volcano_L1.clean.DEG.Chr5.dupli_TARBP2.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "TARBP2", split.by = "Chr5.dupli")
dev.off()

pdf(file.path(figures_dir, "Volcano_L1.clean.DEG.Chr5.dupli_G3BP.pdf"),
    width=6,height=4,paper='special')
VlnPlot(L1.clean, feature = "G3BP1", split.by = "Chr5.dupli")
dev.off()

# ======= identified cells with potential both 3p-loss & 5q.gain ======
Idents(L1.clean) <- "Chr3.loss"
DimPlot(L1.clean)

Idents(L1.clean) <- "Chr5.dupli"
DimPlot(L1.clean)

Idents(L1.clean) <- "seurat_clusters"
L1.clean.3p5qcnv = subset(L1.clean, idents = "0", invert = TRUE)
L1.clean.3p5qcnv@meta.data$cnv.3p5q <- ifelse(
  L1.clean.3p5qcnv$Chr3.loss == "TRUE" & 
    L1.clean.3p5qcnv$Chr5.dupli == "TRUE", "TRUE", "FALSE"
  )

table(L1.clean.3p5qcnv$Chr3.loss)
table(L1.clean.3p5qcnv$Chr5.dupli)

L1.clean.3p5qcnv.meta <- L1.clean.3p5qcnv@meta.data
filtered_cells <- cbind(L1.clean.3p5qcnv$Chr3.loss, L1.clean.3p5qcnv$Chr5.dupli)
filtered_cells.3p5qcnv.selected <- rownames(
  L1.clean.3p5qcnv.meta)[L1.clean.3p5qcnv.meta$Chr3.loss
                             == "TRUE" & L1.clean.3p5qcnv.meta$Chr5.dupli 
                             == "TRUE"]

# ==== using less stringent cut offs
L1.clean$Chr3loss.Chr5.dupli <- ifelse(
  L1.clean@meta.data$proportion_scaled_loss_chr3 > 0 &
  L1.clean@meta.data$proportion_scaled_dupli_chr5 > 0, 
  "Chr3loss.Chr5.dupli_TRUE", "Chr3loss.Chr5.dupli_FALSE"
  )

table(L1.clean$Chr3loss.Chr5.dupli)

save.image("3p5q.RData")

# To restore work environment later use: load("3p5q.RData") 

## ======== plot DEG on the Chr3 and Chr5 using karayoploter package ===========
# https://bernatgel.github.io/karyoploter_tutorial//Examples/GeneExpression/GeneExpression.html 

BiocManager::install("GenomicFeatures")
BiocManager::install("org.Hs.eg.db")

library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicFeatures)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
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
kp <- plotKaryotype(plot.type=2, L1somes = c("chr3"))
kpDataBackground(kp, data.panel = 1, col="#AACBFF")
kpDataBackground(kp, data.panel = 2, col="#FFAACB")

# sort gene by p_val_adj and select top 50 genes. 
ordered <- hg.genes[order(hg.genes$p_val_adj, na.last = TRUE),]

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

sign.genes <- filtered.hg.genes[filtered.hg.genes$p_val_adj < 0.05,]
head(sign.genes)
kp <- plotKaryotype(plot.type=2, L1somes = c("chr3"))
kpDataBackground(kp, data.panel = 1, col="#AACBFF")
kpPoints(kp, data=sign.genes, y=sign.genes$log.pval, ymax=max(sign.genes$log.pval))
kp <- kpPlotMarkers(kp, ordered[1:50], labels = ordered_symbols, text.orientation = "vertical")

kp <- plotKaryotype(plot.type=2, L1somes = c("chr3"))
kpDataBackground(kp, data.panel = 1, col="#AACBFF")
fc.ymax <- ceiling(max(abs(range(sign.genes$avg_log2FC))))
fc.ymin <- -fc.ymax
kpPoints(kp, data=sign.genes, y=sign.genes$avg_log2FC, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin)
kpAddLabels(kp, labels = "avg_log2FC", srt=90, pos=1, label.margin = 0.06, ymax=fc.ymax, ymin=fc.ymin)
kp <- kpPlotMarkers(kp, ordered[1:50], labels = ordered_symbols, text.orientation = "vertical")

# plot significant DE genes across genome

fc.ymax <- ceiling(max(abs(range(sign.genes$avg_log2FC))))
fc.ymin <- -fc.ymax

pdf(file.path(figures_dir, "Karyotype_L1.clean.DEG.Chr3.pdf"),
    width=30,height=30,paper='special')
  kp <- plotKaryotype(plot.type=2)
  kpDataBackground(kp, data.panel = 1, col="#AACBFF")
  kpPoints(kp, data=sign.genes, y=sign.genes$avg_log2FC, ymax=fc.ymax, ymin=fc.ymin)
  kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin)
  kpAddLabels(kp, labels = "avg_log2FC", srt=90, pos=1, label.margin = 0.02, ymax=fc.ymax, ymin=fc.ymin)
dev.off()

#  add p value to the plot. The size of the dot represent the log p-value. Bigger the size = more significant
cex.val <- sqrt(sign.genes$log.pval)/3
kp <- plotKaryotype(genome="hg19")
kpPoints(kp, data=sign.genes, y=sign.genes$avg_log2FC, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin)
kpAddLabels(kp, labels = "avg_log2FC", srt=90, pos=1, label.margin = 0.02, ymax=fc.ymax, ymin=fc.ymin)

# add top 50 gene symbol
pdf(file.path(figures_dir, "Karyotype_L1.clean.DEG.Chr3_top50.pdf"),
    width=30,height=30,paper='special')

kp <- plotKaryotype(genome="hg19")
kpPoints(kp, data=sign.genes, y=sign.genes$avg_log2FC, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin)
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

#Data panel 1
kp <- plotKaryotype(genome="hg19", plot.type=2)
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

pdf(file.path(figures_dir, "Karyotype_L1.clean.DEG.Chr3_top50_AllL1.pdf"),
    width=30,height=30,paper='special')

kp <- plotKaryotype(genome="hg19", plot.type=2, plot.params = pp)
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

# === Only plot Chr3  ===
pdf(file.path(figures_dir, "Karyotype_L1.clean.DEG.Chr3_only_top50.pdf"),
    width=10,height=8,paper='special')

kp <- plotKaryotype(plot.type=2, plot.params = pp, L1somes = c("chr3"))
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
kpAddLabels(kp, labels = "avg_log2FC", srt=90, pos=1, label.margin = 0.06, ymax=fc.ymax, ymin=fc.ymin)
kpPlotMarkers(kp, top.genes, labels = ordered_symbols, text.orientation = "vertical", 
              r0=points.top, label.dist = 0.008, 
              label.color="#444444", 
              line.color = "#777777",
              cex=0.6)

## Data panel 2
kp <- kpPlotDensity(kp, data=hg.genes, window.size = 10e4, data.panel = 2)
kp <- kpPoints(kp, data=GRanges(seqnames="chr3", ranges=IRanges(start=79901053, end=79901054)), 
               r0=0.5, r1=0.5, cex=1.5, col="red", y=0.05)
kp <- kpText(kp, data=GRanges(seqnames="chr3", ranges=IRanges(start=79901053, end=79901054)), 
             labels="gRNA", cex=0.8, vjust=-1, y=0.1)
dev.off()


# ======= DEG list of  DEG_Chr5.dupli True / False =======
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
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

DEG_Chr5.dupli$symbol <- row.names(DEG_Chr5.dupli)
head(DEG_Chr5.dupli)

# Merge DEG data into GRanges by matching "symbol"
# Create a named list from the DEG table, indexed by the symbol
DEG_Chr5.dupli_list <- split(DEG_Chr5.dupli, DEG_Chr5.dupli$symbol)

# Function to fetch DEG data by symbol
fetch_deg_data <- function(symbol) {
  if (symbol %in% names(DEG_Chr5.dupli_list)) {
    return(DEG_Chr5.dupli_list[[symbol]])
  } else {
    # Return NA values if symbol is not found (adjust the length according to the number of DEG columns)
    return(rep(NA, length(DEG_Chr5.dupli_list[[1]]) - 1))
  }
}

# Add DEG data as metadata columns to the GRanges object
mcols(hg.genes) <- cbind(mcols(hg.genes), do.call(rbind, lapply(mcols(hg.genes)$symbol, fetch_deg_data)))

# The GRanges object now includes DEG data as additional metadata columns
head(hg.genes, n=50)

# establish karyotype plot
kp <- plotKaryotype(plot.type=2, L1somes = c("chr5"))
kpDataBackground(kp, data.panel = 1, col="#AACBFF")
kpDataBackground(kp, data.panel = 2, col="#FFAACB")

# sort gene by p_val_adj and select top 50 genes. 
ordered <- hg.genes[order(hg.genes$p_val_adj, na.last = TRUE),]

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

sign.genes <- filtered.hg.genes[filtered.hg.genes$p_val_adj < 0.05,]
head(sign.genes)
kp <- plotKaryotype(plot.type=2, L1somes = c("chr5"))
kpDataBackground(kp, data.panel = 1, col="#AACBFF")
kpPoints(kp, data=sign.genes, y=sign.genes$log.pval, ymax=max(sign.genes$log.pval))
kp <- kpPlotMarkers(kp, ordered[1:50], labels = ordered_symbols, text.orientation = "vertical")

kp <- plotKaryotype(plot.type=2, L1somes = c("chr5"))
kpDataBackground(kp, data.panel = 1, col="#AACBFF")
fc.ymax <- ceiling(max(abs(range(sign.genes$avg_log2FC))))
fc.ymin <- -fc.ymax
kpPoints(kp, data=sign.genes, y=sign.genes$avg_log2FC, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin)
kpAddLabels(kp, labels = "avg_log2FC", srt=90, pos=1, label.margin = 0.06, ymax=fc.ymax, ymin=fc.ymin)
kp <- kpPlotMarkers(kp, ordered[1:50], labels = ordered_symbols, text.orientation = "vertical")

# plot significant DE genes across genome

fc.ymax <- ceiling(max(abs(range(sign.genes$avg_log2FC))))
fc.ymin <- -fc.ymax

pdf(file.path(figures_dir, "Karyotype_L1.clean.DEG.chr5.pdf"),
    width=30,height=30,paper='special')
kp <- plotKaryotype(plot.type=2)
kpDataBackground(kp, data.panel = 1, col="#AACBFF")
kpPoints(kp, data=sign.genes, y=sign.genes$avg_log2FC, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin)
kpAddLabels(kp, labels = "avg_log2FC", srt=90, pos=1, label.margin = 0.02, ymax=fc.ymax, ymin=fc.ymin)
dev.off()

#  add p value to the plot. The size of the dot represent the log p-value. Bigger the size = more significant
cex.val <- sqrt(sign.genes$log.pval)/3
kp <- plotKaryotype(genome="hg19")
kpPoints(kp, data=sign.genes, y=sign.genes$avg_log2FC, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin)
kpAddLabels(kp, labels = "avg_log2FC", srt=90, pos=1, label.margin = 0.02, ymax=fc.ymax, ymin=fc.ymin)

# add top 50 gene symbol
pdf(file.path(figures_dir, "Karyotype_L1.clean.DEG.chr5_top50.pdf"),
    width=30,height=30,paper='special')

kp <- plotKaryotype(genome="hg19")
kpPoints(kp, data=sign.genes, y=sign.genes$avg_log2FC, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin)
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

#Data panel 1
kp <- plotKaryotype(genome="hg19", plot.type=2)
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

pdf(file.path(figures_dir, "Karyotype_L1.clean.DEG.Chr5_top50_AllChrom2.pdf"),
    width=30,height=30,paper='special')

kp <- plotKaryotype(genome="hg19", plot.type=2, plot.params = pp)
kpAddMainTitle(kp, main = "Gene expression - Chr5 dupli True vs False")
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

# === Only plot chr5  ===
pdf(file.path(figures_dir, "Karyotype_L1.clean.DEG.chr5_only_top50.pdf"),
    width=10,height=8,paper='special')

kp <- plotKaryotype(plot.type=2, plot.params = pp, L1somes = c("chr5"))
kpAddMainTitle(kp, main = "Gene expression - chr5 dupli True vs False")
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
kpAddLabels(kp, labels = "avg_log2FC", srt=90, pos=1, label.margin = 0.06, ymax=fc.ymax, ymin=fc.ymin)
kpPlotMarkers(kp, top.genes, labels = ordered_symbols, text.orientation = "vertical", 
              r0=points.top, label.dist = 0.008, 
              label.color="#444444", 
              line.color = "#777777",
              cex=0.6)

## Data panel 2
kp <- kpPlotDensity(kp, data=hg.genes, window.size = 10e4, data.panel = 2)
kp <- kpPoints(kp, data=GRanges(seqnames="chr5", ranges=IRanges(start=79901053, end=79901054)), 
               r0=0.5, r1=0.5, cex=1.5, col="red", y=0.05)
kp <- kpText(kp, data=GRanges(seqnames="chr5", ranges=IRanges(start=79901053, end=79901054)), 
             labels="gRNA", cex=0.8, vjust=-1, y=0.1)
dev.off()



