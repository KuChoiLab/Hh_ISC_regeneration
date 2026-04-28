# 0. Setting ----
mem.maxVSize(60000)
set.seed(1)

BASE_PATH  <- "/Users/choi_lab_sk/Dropbox/2026_Hh_ISC_regeneration/scrnaseq"
RDATA_DIR  <- file.path(BASE_PATH, "RData")
OUTPUT_DIR <- file.path(BASE_PATH, "plot", "Extended_Data_Fig7")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

condition_colors <- c("ENR" = "#90979F", "SAG" = "#B04D9D", "PGE2" = "#DE9B27")

phase_order <- c("qG0", "G1", "Late G1", "S", "S/G2", "G2/M", "M/Early G1")
phase_colors <- c("qG0"        = "#F8766E",
                  "G1"         = "#D19301",
                  "Late G1"    = "#8AAC00",
                  "S"          = "#01B3F6",
                  "S/G2"       = "#9B8DFF",
                  "G2/M"       = "#D178FF",
                  "M/Early G1" = "#FB699F")

potency_order  <- c("Multipotent", "Oligopotent", "Unipotent", "Differentiated")
potency_colors <- c("Multipotent"    = "#F6836D",
                    "Oligopotent"    = "#BEBBDA",
                    "Unipotent"      = "#FFFFB4",
                    "Differentiated" = "#8CD4C8")

# 1. Libraries ----
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# 2. Load data ----
ISC_annotated_seurat <- readRDS(file.path(RDATA_DIR, "ISC_annotated_seurat.rds"))
meta_data <- ISC_annotated_seurat@meta.data

## Add cilia module score
cilia_genes <- c("Ift88", "Ift140", "Ift20", "Kif3a", "Kif3b",
                 "Ttbk2", "Dync2h1", "Arl13b", "Tulp3", "Dynlt1",
                 "Ift57", "Ift172", "Tuba4a")
ISC_annotated_seurat <- AddModuleScore(ISC_annotated_seurat,
                                       features = list(cilia_genes),
                                       name     = "cilia_score")

# 3. Figures ----

## Ext. Data Fig. 7a ----

pdf(file.path(OUTPUT_DIR, "ExtData_Fig7a_UMAP_ccAFv2_by_condition.pdf"),
    width = 18, height = 6)
DimPlot(ISC_annotated_seurat, group.by = "ccAFv2", split.by = "Sample_Name",
        cols = phase_colors, pt.size = 0.3) +
  theme(aspect.ratio = 1)
dev.off()

## Ext. Data Fig. 7b ----
phase_features <- c("qG0", "G1", "Late.G1", "S", "S.G2", "G2.M", "M.Early.G1")
phase_labels   <- c("qG0", "G1", "Late G1", "S", "S/G2", "G2/M", "M/Early G1")

for (i in seq_along(phase_features)) {
  pdf(file.path(OUTPUT_DIR,
                paste0("ExtData_Fig7b_ccAFv2_", phase_features[i],
                       "_VlnPlot_by_cluster.pdf")),
      width = 28, height = 5)
  print(VlnPlot(ISC_annotated_seurat, features = phase_features[i],
                split.by = "Sample_Name", pt.size = 0) +
          scale_fill_manual(values = condition_colors) +
          geom_boxplot(width = 0.5, alpha = 0.5,
                       position = position_dodge(0.9)) +
          labs(title = paste("ccAFv2 score:", phase_labels[i])) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
}

## Ext. Data Fig. 7c ----
library(CytoTRACE2)
annotation_df <- as.data.frame(Idents(ISC_annotated_seurat))
cytotrace2_plots_all <- plotData(cytotrace2_result = ISC_annotated_seurat,
                                 annotation = annotation_df,
                                 is_seurat = TRUE,
                                 pc_dims = 41,
                                 seed = 1)

ct2_plot_data <- cytotrace2_plots_all$CytoTRACE2_UMAP$data
ct2_plot_data$Sample_Name <- ISC_annotated_seurat$Sample_Name[rownames(ct2_plot_data)]

plots_split <- lapply(c("ENR", "SAG", "PGE2"), function(cond) {
  df <- ct2_plot_data[ct2_plot_data$Sample_Name == cond, ]
  p  <- cytotrace2_plots_all$CytoTRACE2_UMAP
  p$data <- df
  p + ggtitle(cond) + theme(aspect.ratio = 1)
})

pdf(file.path(OUTPUT_DIR, "ExtData_Fig7c_CytoTRACE2_Relative_UMAP_split.pdf"),
    width = 18, height = 6)
wrap_plots(plots_split, ncol = 3)
dev.off()

## Ext. Data Fig. 7d ----
cytotrace2_plots_list <- lapply(c("ENR", "SAG", "PGE2"), function(cond) {
  sub <- subset(ISC_annotated_seurat, Sample_Name == cond)
  ann <- as.data.frame(Idents(sub))
  plotData(cytotrace2_result = sub,
           annotation       = ann,
           is_seurat        = TRUE,
           pc_dims          = 41,
           seed             = 1)
})

names(cytotrace2_plots_list) <- c("ENR", "SAG", "PGE2")

pdf(file.path(OUTPUT_DIR, "ExtData_Fig7d_CytoTRACE2_Boxplot_byPheno_split.pdf"),
    width = 30, height = 10)
wrap_plots(
  lapply(c("ENR", "SAG", "PGE2"), function(cond) {
    cytotrace2_plots_list[[cond]]$CytoTRACE2_Boxplot_byPheno +
      ggtitle(cond)
  }),
  ncol = 3
)
dev.off()

## Ext. Data Fig. 7e-g ----
ISC_annotated_seurat_ENR  <- subset(ISC_annotated_seurat, Sample_Name == "ENR")
ISC_annotated_seurat_SAG  <- subset(ISC_annotated_seurat, Sample_Name == "SAG")
ISC_annotated_seurat_PGE2 <- subset(ISC_annotated_seurat, Sample_Name == "PGE2")

gene_panels <- c("e" = "Smo", "f" = "Ptch1", "g" = "Ihh")

for (panel in names(gene_panels)) {
  gene <- gene_panels[[panel]]
  
  p1 <- VlnPlot(ISC_annotated_seurat, features = gene,
                group.by = "Sample_Name", pt.size = 0, alpha = 0.3) +
    geom_boxplot(width = 0.3, alpha = 0.9) +
    scale_fill_manual(values = condition_colors) +
    NoLegend() + theme(aspect.ratio = 1) + ggtitle(gene)
  
  p2 <- VlnPlot(ISC_annotated_seurat_ENR,  features = gene, pt.size = 0) +
    NoLegend() + ggtitle("ENR")
  p3 <- VlnPlot(ISC_annotated_seurat_SAG,  features = gene, pt.size = 0) +
    NoLegend() + ggtitle("SAG")
  p4 <- VlnPlot(ISC_annotated_seurat_PGE2, features = gene, pt.size = 0) +
    NoLegend() + ggtitle("PGE2")
  
  p5 <- FeaturePlot(ISC_annotated_seurat_ENR,  features = gene,
                    order = TRUE, cols = c("lightgray", "darkred"),
                    pt.size = 0.3) +
    NoLegend() + theme(aspect.ratio = 1) + ggtitle("ENR")
  p6 <- FeaturePlot(ISC_annotated_seurat_SAG,  features = gene,
                    order = TRUE, cols = c("lightgray", "darkred"),
                    pt.size = 0.3) +
    NoLegend() + theme(aspect.ratio = 1) + ggtitle("SAG")
  p7 <- FeaturePlot(ISC_annotated_seurat_PGE2, features = gene,
                    order = TRUE, cols = c("lightgray", "darkred"),
                    pt.size = 0.3) +
    NoLegend() + theme(aspect.ratio = 1) + ggtitle("PGE2")
  
  pdf(file.path(OUTPUT_DIR,
                paste0("ExtData_Fig7", panel, "_", gene, ".pdf")),
      width = 45, height = 5)
  print(wrap_plots(p1, p2, p3, p4, p5, p6, p7,
                   nrow   = 1,
                   widths = c(3, 10, 10, 10, 4, 4, 4)))
  dev.off()
}

## Ext. Data Fig. 7h ----
pdf(file.path(OUTPUT_DIR, "ExtData_Fig7h_GSE148528_UMAP.pdf"), width = 18, height = 5)
FeaturePlot(seurat_GSE148528, features = c("Smo","Ptch1","Ihh"), cols = c("lightgray", "darkred")) & theme(aspect.ratio = 1)
dev.off()

## Ext. Data Fig. 7i ----
pdf(file.path(OUTPUT_DIR, "ExtData_Fig7i_cilia_score_VlnPlot_by_cluster.pdf"),
    width = 28, height = 5)
VlnPlot(ISC_annotated_seurat, features = "cilia_score1",
        split.by = "Sample_Name", pt.size = 0) +
  scale_fill_manual(values = condition_colors) +
  geom_boxplot(width = 0.5, alpha = 0.5,
               position = position_dodge(0.9)) +
  labs(title = "Cilia module score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf(file.path(OUTPUT_DIR, "ExtData_Fig7i_cilia_score_VlnPlot_by_condition.pdf"),
    width = 5, height = 5)
VlnPlot(ISC_annotated_seurat, features = "cilia_score1",
        group.by = "Sample_Name", pt.size = 0, alpha = 0.3) +
  geom_boxplot(width = 0.3, alpha = 0.9) +
  scale_fill_manual(values = condition_colors) +
  NoLegend() + theme(aspect.ratio = 1) + ggtitle("Cilia module score")
dev.off()