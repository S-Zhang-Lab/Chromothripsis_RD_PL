### Additional analysis (polished figures)
### 
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

library(dplyr)

obj@meta.data <- obj@meta.data %>%
  mutate(
    sg_group = case_when(
      # Lane 1 (Day 5)
      HTO_classification %in% c("L1-H-HTO-5", "L1-H-HTO-6") ~ "sgNTC",
      HTO_classification %in% c("L1-H-HTO-2", "L1-H-HTO-3") ~ "sgYq",
      
      # Lane 2 (Day 5)
      HTO_classification %in% c("L2-H-HTO-2", "L2-H-HTO-3") ~ "sg3p",
      
      # Lane 3 (Day 10)
      HTO_classification %in% c("L3-H-HTO-2", "L3-H-HTO-3") ~ "sg3p",
      HTO_classification %in% c("L3-H-HTO-4", "L3-H-HTO-5") ~ "sgYq",
      HTO_classification %in% c("L3-H-HTO-1", "L3-H-HTO-6") ~ "sgNTC",
      
      TRUE ~ NA_character_
    ),
    
    timepoint = case_when(
      grepl("^L1|^L2", HTO_classification) ~ "day5",
      grepl("^L3", HTO_classification) ~ "day10",
      TRUE ~ NA_character_
    )
  )
obj$sg_group <- obj@meta.data$sg_group
obj$timepoint <- factor(obj@meta.data$timepoint, levels = c("day5", "day10"))

obj@meta.data$is_3p_loss <- ifelse(
  obj@meta.data$recoded_loss_chr3_2_bin == "02_chr3_loss",
  1,
  0
)

summary_df <- obj@meta.data %>%
  filter(!is.na(sg_group)) %>%
  group_by(sg_group, timepoint) %>%
  summarise(
    total_cells = n(),
    n_3p_loss = sum(is_3p_loss),
    percent_3p_loss = (n_3p_loss / total_cells) * 100,
    .groups = "drop"
  )

print(summary_df)

library(ggplot2)

summary_df$sg_group <- factor(summary_df$sg_group, 
                              levels = c("sgNTC", "sg3p", "sgYq"))

summary_df$sg_group <- factor(
  summary_df$sg_group,
  levels = c("sgNTC", "sg3p", "sgYq")
)

ggplot(summary_df,
       aes(x = sg_group,
           y = percent_3p_loss,
           fill = timepoint)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.7),
           width = 0.6) +
  scale_fill_manual(values = c("day5" = "red",
                               "day10" = "gray")) +
  labs(title = "Percentage of Cells with 3p Loss",
       x = "sgRNA",
       y = "% Cells with 3p Loss") +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank()
  )


summary_overall <- obj@meta.data %>%
  filter(!is.na(sg_group)) %>%
  group_by(sg_group, timepoint) %>%
  summarise(
    total_cells = n(),
    n_3p_loss = sum(recoded_loss_chr3_2_bin == "02_chr3_loss"),
    percent_3p_loss = (n_3p_loss / total_cells) * 100,
    .groups = "drop"
  )

summary_overall

write.csv(summary_overall,
          "Overall_percent_3p_loss_by_sg_and_timepoint.csv",
          row.names = FALSE)

summary_replicates <- obj@meta.data %>%
  filter(!is.na(sg_group)) %>%
  group_by(sg_group, timepoint, HTO_classification) %>%
  summarise(
    total_cells = n(),
    n_3p_loss = sum(recoded_loss_chr3_2_bin == "02_chr3_loss"),
    percent_3p_loss = (n_3p_loss / total_cells) * 100,
    .groups = "drop"
  )

summary_replicates

write.csv(summary_replicates,
          "Replicate_percent_3p_loss_by_HTO.csv",
          row.names = FALSE)

summary_replicates <- summary_replicates %>%
  mutate(
    lane = sub("-H-.*", "", HTO_classification),
    hto  = sub(".*HTO-", "HTO-", HTO_classification)
  )

summary_stats <- summary_replicates %>%
  group_by(sg_group, timepoint) %>%
  summarise(
    mean_percent = mean(percent_3p_loss),
    sd_percent = sd(percent_3p_loss),
    n_reps = n(),
    .groups = "drop"
  )

write.csv(summary_stats,
          "Replicate_summary_stats_3p_loss.csv",
          row.names = FALSE)

# Plot barchart with biological replicates

summary_stats$sg_group <- factor(
  summary_stats$sg_group,
  levels = c("sgNTC", "sg3p", "sgYq")
)

summary_stats$timepoint <- factor(
  summary_stats$timepoint,
  levels = c("day5", "day10")
)

summary_replicates$sg_group <- factor(
  summary_replicates$sg_group,
  levels = c("sgNTC", "sg3p", "sgYq")
)

summary_replicates$timepoint <- factor(
  summary_replicates$timepoint,
  levels = c("day5", "day10")
)

library(ggplot2)

dodge_width <- 0.7

ggplot(summary_stats,
       aes(x = sg_group,
           y = mean_percent,
           fill = timepoint)) +
  
  # Bars (mean of replicates)
  geom_bar(stat = "identity",
           position = position_dodge(width = dodge_width),
           width = 0.6,
           color = "black") +
  
  # Error bars (SD)
  geom_errorbar(aes(ymin = mean_percent - sd_percent,
                    ymax = mean_percent + sd_percent),
                width = 0.2,
                position = position_dodge(width = dodge_width),
                size = 0.8) +
  
  # Replicate dots
  geom_point(data = summary_replicates,
             aes(x = sg_group,
                 y = percent_3p_loss,
                 group = timepoint),
             position = position_jitterdodge(
               jitter.width = 0.05,
               dodge.width = dodge_width),
             color = "black",
             size = 2.5) +
  
  scale_fill_manual(values = c("day5" = "#d73027",   # bright red
                               "day10" = "gray70")) +
  
  labs(x = "sgRNA",
       y = "% Cells with 3p Loss") +
  
  theme_classic(base_size = 14) +   # removes gridlines
  
  theme(
    legend.title = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    legend.position = "right"
  )


## Venn diagrams (updated)
day5_obj <- subset(obj, subset = timepoint == "day5")
Idents(day5_obj) <- "sg_group"
day5_obj <- PrepSCTFindMarkers(day5_obj)

#DEG analysis
deg_sgNTC <- FindMarkers(
  day5_obj,
  ident.1 = "sgNTC",
  ident.2 = c("sg3p", "sgYq"),
  logfc.threshold = 0.25,
  min.pct = 0.1
)

deg_sg3p <- FindMarkers(
  day5_obj,
  ident.1 = "sg3p",
  ident.2 = c("sgNTC", "sgYq"),
  logfc.threshold = 0.25,
  min.pct = 0.1
)

deg_sgYq <- FindMarkers(
  day5_obj,
  ident.1 = "sgYq",
  ident.2 = c("sgNTC", "sg3p"),
  logfc.threshold = 0.25,
  min.pct = 0.1
)

pval_cutoff <- 0.05
logfc_cutoff <- 0.25

genes_sgNTC <- rownames(subset(deg_sgNTC,
                               p_val_adj < pval_cutoff &
                                 avg_log2FC > logfc_cutoff))

genes_sg3p <- rownames(subset(deg_sg3p,
                              p_val_adj < pval_cutoff &
                                avg_log2FC > logfc_cutoff))

genes_sgYq <- rownames(subset(deg_sgYq,
                              p_val_adj < pval_cutoff &
                                avg_log2FC > logfc_cutoff))


library(VennDiagram)
library(grid)

venn_list <- list(
  sgNTC = genes_sgNTC,
  sg3p  = genes_sg3p,
  sgYq  = genes_sgYq
)

venn.plot <- venn.diagram(
  x = venn_list,
  category.names = c("sgNTC", "sg3p", "sgYq"),
  fill = c("gray", "red", "darkred"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.5,
  filename = NULL,
  lwd = 2
)

grid.draw(venn.plot)

genes_sgNTC
genes_sg3p
genes_sgYq

# Pairwise overlaps
overlap_NTC_3p  <- intersect(genes_sgNTC, genes_sg3p)
overlap_NTC_Yq  <- intersect(genes_sgNTC, genes_sgYq)
overlap_3p_Yq   <- intersect(genes_sg3p, genes_sgYq)

# Triple overlap
overlap_all <- Reduce(intersect, list(genes_sgNTC, genes_sg3p, genes_sgYq))

unique_sgNTC <- setdiff(
  genes_sgNTC,
  union(overlap_NTC_3p, union(overlap_NTC_Yq, overlap_all))
)

unique_sg3p <- setdiff(
  genes_sg3p,
  union(overlap_NTC_3p, union(overlap_3p_Yq, overlap_all))
)

unique_sgYq <- setdiff(
  genes_sgYq,
  union(overlap_NTC_Yq, union(overlap_3p_Yq, overlap_all))
)

write.table(unique_sgNTC, "Day5_unique_sgNTC_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(unique_sg3p, "Day5_unique_sg3p_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(unique_sgYq, "Day5_unique_sgYq_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(overlap_NTC_3p, "Day5_overlap_sgNTC_sg3p.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(overlap_NTC_Yq, "Day5_overlap_sgNTC_sgYq.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(overlap_3p_Yq, "Day5_overlap_sg3p_sgYq.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(overlap_all, "Day5_overlap_all_three.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)


### Venn for 3p loss vs sg3p vs sgYq
day5_obj <- subset(obj, subset = timepoint == "day5")
table(day5_obj$sg_group)
table(day5_obj$recoded_loss_chr3_2_bin)
day5_obj <- PrepSCTFindMarkers(day5_obj)

Idents(day5_obj) <- "recoded_loss_chr3_2_bin"

deg_3p_loss_day5 <- FindMarkers(
  day5_obj,
  ident.1 = "02_chr3_loss",
  ident.2 = "01_No_loss",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
Idents(day5_obj) <- "sg_group"

deg_sg3p_day5 <- FindMarkers(
  day5_obj,
  ident.1 = "sg3p",
  ident.2 = "sgNTC",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
deg_sgYq_day5 <- FindMarkers(
  day5_obj,
  ident.1 = "sgYq",
  ident.2 = "sgNTC",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
pval_cutoff <- 0.05
logfc_cutoff <- 0.25

genes_3p_loss <- rownames(subset(
  deg_3p_loss_day5,
  p_val_adj < pval_cutoff & avg_log2FC > logfc_cutoff
))

genes_sg3p <- rownames(subset(
  deg_sg3p_day5,
  p_val_adj < pval_cutoff & avg_log2FC > logfc_cutoff
))

genes_sgYq <- rownames(subset(
  deg_sgYq_day5,
  p_val_adj < pval_cutoff & avg_log2FC > logfc_cutoff
))
length(genes_3p_loss)
length(genes_sg3p)
length(genes_sgYq)


library(VennDiagram)
library(grid)

venn_list <- list(
  "3p Loss" = genes_3p_loss,
  "sg3p"    = genes_sg3p,
  "sgYq"    = genes_sgYq
)

venn.plot <- venn.diagram(
  x = venn_list,
  fill = c("pink", "lightgreen", "gray"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.5,
  filename = NULL,
  lwd = 2
)

grid.draw(venn.plot)

# Pairwise overlaps
overlap_3p_sg3p <- intersect(genes_3p_loss, genes_sg3p)
overlap_3p_sgYq <- intersect(genes_3p_loss, genes_sgYq)
overlap_sg3p_sgYq <- intersect(genes_sg3p, genes_sgYq)

# Triple overlap
overlap_all <- Reduce(intersect, list(
  genes_3p_loss,
  genes_sg3p,
  genes_sgYq
))

# Unique sets
unique_3p_loss <- setdiff(
  genes_3p_loss,
  union(overlap_3p_sg3p,
        union(overlap_3p_sgYq, overlap_all))
)

unique_sg3p <- setdiff(
  genes_sg3p,
  union(overlap_3p_sg3p,
        union(overlap_sg3p_sgYq, overlap_all))
)

unique_sgYq <- setdiff(
  genes_sgYq,
  union(overlap_3p_sgYq,
        union(overlap_sg3p_sgYq, overlap_all))
)

#export
sig_3p_loss <- deg_3p_loss_day5 %>%
  tibble::rownames_to_column("gene") %>%
  filter(p_val_adj < pval_cutoff & avg_log2FC > logfc_cutoff) %>%
  arrange(p_val_adj)

sig_sg3p <- deg_sg3p_day5 %>%
  tibble::rownames_to_column("gene") %>%
  filter(p_val_adj < pval_cutoff & avg_log2FC > logfc_cutoff) %>%
  arrange(p_val_adj)

sig_sgYq <- deg_sgYq_day5 %>%
  tibble::rownames_to_column("gene") %>%
  filter(p_val_adj < pval_cutoff & avg_log2FC > logfc_cutoff) %>%
  arrange(p_val_adj)

write.csv(sig_3p_loss,
          "Day5_sig_ranked_3p_loss_venn2.csv",
          row.names = FALSE)

write.csv(sig_sg3p,
          "Day5_sig_ranked_sg3p_venn2.csv",
          row.names = FALSE)

write.csv(sig_sgYq,
          "Day5_sig_ranked_sgYq_venn2.csv",
          row.names = FALSE)

library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)

hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

create_rank_vector <- function(deg_table) {
  ranks <- deg_table$avg_log2FC
  names(ranks) <- rownames(deg_table)
  ranks <- sort(ranks, decreasing = TRUE)
  return(ranks)
}

ranks_3p_loss <- create_rank_vector(deg_3p_loss_day5)
ranks_sg3p    <- create_rank_vector(deg_sg3p_day5)
ranks_sgYq    <- create_rank_vector(deg_sgYq_day5)

fgsea_3p_loss <- fgsea(
  pathways = hallmark_sets,
  stats    = ranks_3p_loss,
  nperm    = 10000
)

fgsea_sg3p <- fgsea(
  pathways = hallmark_sets,
  stats    = ranks_sg3p,
  nperm    = 10000
)

fgsea_sgYq <- fgsea(
  pathways = hallmark_sets,
  stats    = ranks_sgYq,
  nperm    = 10000
)

plot_fgsea <- function(fgsea_res, title) {
  fgsea_res %>%
    arrange(padj) %>%
    slice(1:15) %>%
    mutate(pathway = reorder(pathway, NES)) %>%
    ggplot(aes(x = pathway, y = NES, fill = padj < 0.05)) +
    geom_col() +
    coord_flip() +
    theme_classic(base_size = 14) +
    labs(title = title,
         x = "",
         y = "Normalized Enrichment Score (NES)") +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray70")) +
    theme(legend.position = "none")
}

plot_fgsea(fgsea_3p_loss, "3p Loss vs No Loss (Hallmark)")
plot_fgsea(fgsea_sg3p,    "sg3p vs sgNTC (Hallmark)")
plot_fgsea(fgsea_sgYq,    "sgYq vs sgNTC (Hallmark)")
