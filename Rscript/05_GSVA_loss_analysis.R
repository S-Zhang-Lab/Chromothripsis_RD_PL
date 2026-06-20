############################################################## 
#          GSVA analysis for R02 CITE-seq                    # 
############################################################## 
# loading packages and gene sets 
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(GSEABase) 
library(GSVA) 
library(GSVAdata) 
library(Biobase) 
library(limma) 
library(RColorBrewer) 
library(parallel) 

setwd("/path/to/Chromothripsis_2")  # <-- set to your local path to the project root
Mets <- readRDS("day5_cells_EA.rds")
# later repeat the same process with day10_cells_EA.rds
Mets <- readRDS("day10_cells_EA.rds")

# get gene sets
gs=getGmt("h.all.v2024.1.Hs.symbols.gmt") 

# v7 GSEA MSigDB gene sets, gene symbols, All genesets 
summary(gs) 
gs.full <- gs 

# run ssGSEA on all genesets, if not previously run. Otherwise can load directly. See the next step.  
ssGSEA_result_gs.full = gsva(as.matrix(Mets@assays$SCT$data), gs.full, min.sz=5, max.sz=500, verbose=TRUE, method="ssgsea") 
write.csv(ssGSEA_result_gs.full, file = "ssGSEA_result_gsFull_day10.csv") 


# Min-Max scaling 
# Load necessary library
library(scales)

# Assuming ssGSEA_result_gs.full contains GSVA scores
# Apply Min-Max scaling to transform scores to range [0, 1]
min_max_scaled_scores <- apply(ssGSEA_result_gs.full, 2, function(x) {
  rescaled_x <- scales::rescale(x, to = c(0, 1), from = range(x, na.rm = TRUE))
  return(rescaled_x)
})


# Create a new Seurat object with the Min-Max scaled GSVA scores
gsvasc <- CreateSeuratObject(counts = min_max_scaled_scores, project = "ssGSEA_result")
gsvasc 

any(gsvasc@assays$RNA$counts < 0) # check whether data contain negative values
any(gsvasc@assays$RNA$counts == 0)

gsva.meta <- gsvasc@meta.data

# Merge the meta
gsvasc@meta.data <- Mets@meta.data

head(gsvasc@meta.data) 
colnames(gsvasc@meta.data)

# Normalization
gsvasc <- NormalizeData(gsvasc)

# find variable features
gsvasc <- FindVariableFeatures(gsvasc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(gsvasc)
gsvasc <- ScaleData(gsvasc, features = all.genes)

plot1 <- VariableFeaturePlot(gsvasc)
plot1

# Set seed and create copy of R object because we will reuse the same code for day 10
set.seed(42)
gsvasc.2 <- gsvasc
# gsvasc.2 can be day5 cells, gsvasc day10

# Examine top marker pathways for each met cluster
Idents(gsvasc) <- "seurat_clusters"
gs.cluster.markers <- FindAllMarkers(gsvasc, only.pos = TRUE) 
gs.cluster.markers  %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) -> gs.cluster.markers.top50 
write.csv(as.matrix(gs.cluster.markers.top50), file = "gs.cluster_markers_top50_day10.csv") 

# check individual pathway expression
Idents(gsvasc) <- "recoded_loss_chr3_2_bin"
VlnPlot(gsvasc, features = c("HALLMARK-P53-PATHWAY", "HALLMARK-TNFA-SIGNALING-VIA-NFKB", "HALLMARK-KRAS-SIGNALING-UP", "HALLMARK-INFLAMMATORY-RESPONSE"))

# since cluster 1, enriched in 3p loss, does not have any DEG pathways, I will check the same
# pathways for day10 as I did for day 5

Idents(gsvasc) <- "recoded_loss_chr3_2_bin"
day10_DEGs <- FindMarkers(gsvasc, ident.1 = "02_chr3_loss", ident.2 = "01_No_loss", verbose = TRUE, logfc.threshold = 0.05)
# View top results
head(day10_DEGs)

# For day5 DEG pathways between loss and no loss
day10_DEGs$log2FC <- day10_DEGs$avg_log2FC
day10_DEGs$neg_log10_pval <- -log10(day10_DEGs$p_val_adj)

library(ggplot2)

# save progress
saveRDS(gsvasc.2, file = "Day5_gsvasc.rds")
saveRDS(gsvasc, file = "Day10_gsvasc.rds")
 
write.csv(day10_DEGs, file = "DEG_day10_gsvasc.csv") 


###Plot GSVA data in barplot/cleveland plot style
#1. combine DEGS lists acquired - up in chr3 loss and down in chr3 loss
#2. then create barplot or cleveland plot to show GS enrichment for each condition between time points

library(ggplot2)
library(dplyr)
library(readxl)

#combine DEGS lists acquired for control & 3p loss
#import excel files generated when measuring significant DEGS between time points per conditions 

DEG_down <- read_excel("DEG_day10_gsvasc_down.xlsx")
DEG_up <- read_excel("DEG_day10_gsvasc_up.xlsx")


#change column name from "...1" to "gene_set"
#the "...1" column name only happens if you import a DEG table that you had exported after running FindMarkers
colnames(DEG_down)[colnames(DEG_down)=="...1"]<-"gene_set"
colnames(DEG_up)[colnames(DEG_up)=="...1"]<-"gene_set"

#add condition information to each table
DEG_down$Condition <- "01_No_Loss"
DEG_up$Condition <- "02_3p_Loss"

#combine tables
comb_DEGS <- rbind(DEG_down, DEG_up)

#create Cleveland/Lollipop plot with combined DEGS information

#reduce list for plot by filtering out DEGS with > abs value 0.1 avg_log2FC
# for expansion comparison
comb_DEGS <- comb_DEGS %>% filter(abs(avg_log2FC) >= 0.1)

# Define colors for each condition
condition_colors <- c("01_No_Loss" = "gray", "02_3p_Loss" = "red")

# Create the plot with specific colors for each condition
#you can edit this part of the code if you'd rather plot p val on the x-axis rather than Log FC
#alternatively, you could have Log FC on the x-axis but have the bars colored by p value (I tried this but thought it looked weird/not very infomative because all the p-values were very low for my filtered gene set list)
ggplot(comb_DEGS, aes(x = avg_log2FC, y = reorder(gene_set, -avg_log2FC))) +
  geom_segment(aes(xend = 0, yend = reorder(gene_set, avg_log2FC), color = Condition), linewidth = 0.5, size = 2.5) +
  geom_point(aes(color = Condition)) +
  scale_color_manual(values = condition_colors) +  # Assign specific colors for each condition
  facet_grid(Condition ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Average Log2FC", y = "Gene Set", color = "Condition") +
  theme_minimal()

### Hallmark only ###
#do whole object
##### SF METHOD for GSEA #
# Load required packages specific for pathway enrichment
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)

# Sort by log fold change to get up-regulated and down-regulated genes
upregulated_genes_p <- deg_3p_loss %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  head(180)

upregulated_genes_q <- deg_y_loss %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  head(180)

# or repeat this with 3p vs y loss comparison
upregulated_genes_p <- deg_3p_loss_vs_y %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  head(200)

upregulated_genes_q <- deg_3p_loss_vs_y %>%
  filter(avg_log2FC < 0, p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  head(200)

p_up <- upregulated_genes_p
q_up <- upregulated_genes_q
p_up$gene <- rownames(p_up)  # Correctly assign gene symbols
q_up$gene <- rownames(q_up)  # Correctly assign gene symbols

# Validate symbols and convert to Entrez IDs
valid_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
filtered_p <- intersect(p_up$gene, valid_symbols)
entrez_ids_p_up <- bitr(filtered_p, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_p_up <- na.omit(entrez_ids_p_up$ENTREZID)

filtered_q <- intersect(q_up$gene, valid_symbols)
entrez_ids_q_up <- bitr(filtered_q, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_q_up <- na.omit(entrez_ids_q_up$ENTREZID)

# Get Hallmark gene sets for mouse (or human, adjust species if needed)
# 1. Load Hallmark gene sets for human
library(msigdbr)
library(dplyr)

# Load Hallmark gene sets for human
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)


# 2. Convert your gene list from ENTREZ to SYMBOL (as you're doing)
entrez_p_up_genes <- bitr(entrez_p_up, fromType = "ENTREZID", 
                           toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL

entrez_q_up_genes <- bitr(entrez_q_up, fromType = "ENTREZID", 
                          toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL


# 3. Run enrichment analysis with SYMBOLs
p_up_hallmark <- enricher(
  gene = entrez_p_up_genes,
  TERM2GENE = hallmark_df,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
  
)

# 4. Plot the top pathways
dotplot(p_up_hallmark, showCategory = 10, title = "Top Hallmark Pathways (3p loss)")

q_up_hallmark <- enricher(
  gene = entrez_q_up_genes,
  TERM2GENE = hallmark_df,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
  
)

# 4. Plot the top pathways
dotplot(q_up_hallmark, showCategory = 10, title = "Top Hallmark Pathways (Yq loss)")






