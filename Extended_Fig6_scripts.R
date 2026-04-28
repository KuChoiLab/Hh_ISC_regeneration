# 0. Setting ----
mem.maxVSize(60000)
set.seed(1)

BASE_PATH  <- "/Users/choi_lab_sk/Dropbox/2026_Hh_ISC_regeneration/scrnaseq"
RDATA_DIR  <- file.path(BASE_PATH, "RData")
OUTPUT_DIR <- file.path(BASE_PATH, "plot", "Extended_Data_Fig6")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

condition_colors <- c("ENR" = "#90979F", "SAG" = "#B04D9D", "PGE2" = "#DE9B27")

cluster_order <- c("Canonical stem cell",
                   "Proliferating stem cell (S)",
                   "Proliferating stem cell (G2M)",
                   "Terminally proliferated stem cell",
                   "Mito high stem cell",
                   "Absorptive progenitor",
                   "Ly6a+ enterocyte progenitor",
                   "Ly6a+ enterocyte (early)",
                   "Ly6a+ enterocyte (late)",
                   "Early enterocyte",
                   "Enterocyte",
                   "Enteroendocrine progenitor",
                   "Enteroendocrine",
                   "Secretory progenitor",
                   "Early secretory cell",
                   "Secretory cell",
                   "Tuft cell")

# 1. Libraries ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# 2. Load data ----
ISC_annotated_seurat <- readRDS(file.path(RDATA_DIR, "ISC_annotated_seurat.rds"))

# 3. Figures ----

## Ext. Data Fig. 6a ----
all_markers <- FindAllMarkers(ISC_annotated_seurat, logfc.threshold = 0.5, min.pct = 0.3)

all_markers <- all_markers %>% mutate(statistics = qnorm(p_val/2, lower.tail=F) * sign(avg_log2FC),
                                                  max_stat_second = max(statistics[statistics != Inf]),
                                                  min_stat_second = min(statistics[statistics != -Inf]),
                                                  stat = ifelse(statistics == Inf, max_stat_second + 1,
                                                                ifelse(statistics == -Inf,  min_stat_second - 1, statistics)))
top5_genes <- all_markers %>% filter(p_val_adj < 0.05) %>% 
  group_by(cluster) %>% 
  arrange(desc(stat), desc(avg_log2FC)) %>% 
  group_by(cluster) %>% slice_head(n=5)

pdf(file.path(OUTPUT_DIR, "ExtData_Fig6a_top5_markers_heatmap.pdf"),
    width = 16, height = 10)
DoHeatmap(ISC_annotated_seurat, features = top5_genes$gene, size = 3) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(7, "RdYlBu"))) + theme(aspect.ratio = 1)
dev.off()

## Ext. Data Fig. 6b ----
## UMAP split by condition

pdf(file.path(OUTPUT_DIR, "ExtData_Fig6b_UMAP_split_by_condition.pdf"),
    width = 6, height = 18)
DimPlot(ISC_annotated_seurat, split.by = "Sample_Name",
        label = TRUE, label.size = 2.5, pt.size = 0.2, repel = TRUE, ncol = 1) + theme(aspect.ratio = 1) + NoLegend()
dev.off()

## Ext. Data Fig. 6c ----
make_lineage_vlnplot <- function(obj, cluster_name, genes, out_file,
                                 width = 10, height = 5) {
  cells <- rownames(obj@meta.data)[obj@meta.data$idents_clusters == cluster_name]
  sub   <- obj[, cells]
  plots <- lapply(genes, function(g) {
    VlnPlot(sub, features = g, group.by = "Sample_Name",
            pt.size = 0.01, alpha = 0.9) +
      geom_boxplot(alpha = 0.9, width = 0.5,
                   position = position_dodge(width = 0.75)) +
      scale_fill_manual(values = condition_colors) +
      ggtitle(g) + NoLegend() + theme(aspect.ratio = 1)
  })
  pdf(out_file, width = width, height = height)
  print(wrap_plots(plots, ncol = length(genes)))
  dev.off()
}

## Enterocyte: Alpi, Hes1
make_lineage_vlnplot(ISC_annotated_seurat, "Enterocyte",
                     genes    = c("Alpi", "Hes1"),
                     out_file = file.path(OUTPUT_DIR,
                                          "ExtData_Fig6c_VlnPlot_Enterocyte.pdf"))

## Secretory cell — Goblet markers: Muc2, Spdef
make_lineage_vlnplot(ISC_annotated_seurat, "Secretory cell",
                     genes    = c("Muc2", "Spdef"),
                     out_file = file.path(OUTPUT_DIR,
                                          "ExtData_Fig6c_VlnPlot_Goblet.pdf"))

## Secretory cell — Paneth markers: Lyz1, Sox9
make_lineage_vlnplot(ISC_annotated_seurat, "Secretory cell",
                     genes    = c("Lyz1", "Sox9"),
                     out_file = file.path(OUTPUT_DIR,
                                          "ExtData_Fig6c_VlnPlot_Paneth.pdf"))

## Tuft: Dclk1, Pou2f3
make_lineage_vlnplot(ISC_annotated_seurat, "Tuft",
                     genes    = c("Dclk1", "Pou2f3"),
                     out_file = file.path(OUTPUT_DIR,
                                          "ExtData_Fig6c_VlnPlot_Tuft.pdf"))

## Enteroendocrine: Chga, Neurog3
make_lineage_vlnplot(ISC_annotated_seurat, "Enteroendocrine",
                     genes    = c("Chga", "Neurog3"),
                     out_file = file.path(OUTPUT_DIR,
                                          "ExtData_Fig6c_VlnPlot_Enteroendocrine.pdf"))