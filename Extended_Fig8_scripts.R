# 0. Setting ----
mem.maxVSize(60000)
set.seed(1)

BASE_PATH  <- "/Users/choi_lab_sk/Dropbox/2026_Hh_ISC_regeneration/scrnaseq"
RDATA_DIR  <- file.path(BASE_PATH, "RData")
OUTPUT_DIR <- file.path(BASE_PATH, "plot", "Extended_Data_Fig8")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

condition_colors <- c("ENR" = "#90979F", "SAG" = "#B04D9D", "PGE2" = "#DE9B27")
my_comparisons   <- list(c("ENR", "SAG"), c("SAG", "PGE2"), c("ENR", "PGE2"))

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

stem_clusters <- c("Canonical stem cell",
                   "Proliferating stem cell (S)",
                   "Proliferating stem cell (G2M)",
                   "Terminally proliferated stem cell",
                   "Absorptive progenitor",
                   "Enteroendocrine progenitor",
                   "Secretory progenitor")

# 1. Libraries ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(ggpubr)

# 2. Load data ----
ISC_annotated_seurat <- readRDS(file.path(RDATA_DIR, "ISC_annotated_seurat.rds"))
ISC_annotated_seurat$idents_clusters <- Idents(ISC_annotated_seurat)
ISC_annotated_seurat$idents         <- Idents(ISC_annotated_seurat)

## Gene sets ----
rISC_genes <- c("Tert", "Bmi1", "Hopx", "Lrig1", "Msi1", "Sox9")

isthmus_genes <- c("Fgfbp1", "Atad2", "Birc5", "Stmn1", "Lgr4", "Mki67")

SSC2a_genes <- c("Hmgb2", "Tuba1b", "Ncl", "Nap1l1", "Hmgb1", "Npm1", "Set",
                 "Dut", "Anp32b", "Hspd1", "Dek", "Ranbp1", "Ptma", "Ifitm3",
                 "Hspe1", "Hsp90aa1", "Lyar", "Ube2c", "Top2a", "Tubb5")
SSC2b_genes <- c("S100a11", "Tmsb4x", "Sprr2a3", "Tmsb10", "Rps26", "Rps20",
                 "Rpl36a", "S100a6", "Eif3e", "Rpl12", "Gabarap", "Rpl10")
SSC2c_genes <- c("Malat1", "Ly6d", "Clu", "Cldn4", "Sprr1a", "Areg", "Lamc2",
                 "Il1rn", "Il33", "Ly6a", "Plaur", "Itgb6", "Suox", "Anxa3",
                 "Krt7", "Birc5", "Stmn1", "Atad2")

aVEC_genes <- c("Clu", "Msln", "Sprr1a", "Adam8", "Plaur", "Ly6d", "Plat",
                "Acsbg1", "Mal", "Anxa1", "Isg20", "Anxa5", "Inhba", "Vill",
                "Cd44", "Emp2", "F3", "Hmox1", "Il11", "Krt7", "Ap1s2",
                "Bmp8b", "Pdlim7", "Aldh1a3", "Dok2", "Egfr", "Slc16a3",
                "Slc39a6", "Tgm2", "Serpine1", "Ppfibp1", "Dusp4", "Slc15a3",
                "Glis2", "Cxcr4", "Ephb2", "Cldn4", "Ier3", "Gcnt2", "Thbs1",
                "Trim15", "Bdh2", "Enc1", "Abl2", "Arf2", "Spire1", "Uchl1",
                "Mthfd1l", "Ecm1", "Parvb", "Trim29", "Myadm", "Itga2", "Crip2",
                "Fgd3", "Tnfrsf22", "Rap2b", "Lox", "Akr1b8", "Prss12",
                "Rbp1", "Tnfrsf12a", "Phlda1", "Arid5a", "Marcksl1", "Osmr",
                "Actn1", "Dusp1", "Pls3", "Phldb2", "Phlda3")

damage_resistance_genes <- c("Atm", "Atr", "Mgmt", "Sod1", "Sod2",
                             "Abcb1a", "Abcb1b", "Foxo3")
damage_response_genes   <- c("Trp53", "Bbc3", "Bax", "Noxa1",
                             "Gadd45a", "Gadd45b", "Gadd45g",
                             "Cdkn1a", "Chek1", "Chek2", "Ddb2", "Xpc")
HR_genes   <- c("Brca1", "Brca2", "Rad51")

## Add module scores ----
ISC_annotated_seurat <- AddModuleScore(ISC_annotated_seurat,
                                       features = list(rISC_genes),
                                       name     = "rISC_genes")
ISC_annotated_seurat <- AddModuleScore(ISC_annotated_seurat,
                                       features = list(isthmus_genes),
                                       name     = "isthmus_genes")
ISC_annotated_seurat <- AddModuleScore(ISC_annotated_seurat,
                                       features = list(SSC2a_genes),
                                       name     = "SSC2a_genes")
ISC_annotated_seurat <- AddModuleScore(ISC_annotated_seurat,
                                       features = list(SSC2b_genes),
                                       name     = "SSC2b_genes")
ISC_annotated_seurat <- AddModuleScore(ISC_annotated_seurat,
                                       features = list(SSC2c_genes),
                                       name     = "SSC2c_genes")
ISC_annotated_seurat <- AddModuleScore(ISC_annotated_seurat,
                                       features = list(aVEC_genes),
                                       name     = "aVEC_genes")
ISC_annotated_seurat <- AddModuleScore(ISC_annotated_seurat,
                                       features = list(damage_resistance_genes),
                                       name     = "damage_resistance")
ISC_annotated_seurat <- AddModuleScore(ISC_annotated_seurat,
                                       features = list(damage_response_genes),
                                       name     = "damage_response")
ISC_annotated_seurat <- AddModuleScore(ISC_annotated_seurat,
                                       features = list(HR_genes),
                                       name     = "HR_genes")

ISC_annotated_seurat_stem <- ISC_annotated_seurat[,
                                                  ISC_annotated_seurat$idents %in% stem_clusters]
meta_data      <- ISC_annotated_seurat@meta.data
meta_data_stem <- ISC_annotated_seurat_stem@meta.data

# 3. Figures ----

## Ext. Data Fig. 8a ----
for (gene in isthmus_genes) {
  pdf(file.path(OUTPUT_DIR,
                paste0("ExtData_Fig8a_Isthmus_", gene, "_FeaturePlot.pdf")),
      width = 18, height = 6)
  print(FeaturePlot(ISC_annotated_seurat, features = gene,
                    split.by = "Sample_Name", order = TRUE,
                    cols = c("lightgray", "darkred"), pt.size = 0.3) &
          theme(aspect.ratio = 1))
  dev.off()
}

for (gene in rISC_genes) {
  pdf(file.path(OUTPUT_DIR,
                paste0("ExtData_Fig8a_rISC_", gene, "_FeaturePlot.pdf")),
      width = 18, height = 6)
  print(FeaturePlot(ISC_annotated_seurat, features = gene,
                    split.by = "Sample_Name", order = TRUE,
                    cols = c("lightgray", "darkred"), pt.size = 0.3) &
          theme(aspect.ratio = 1))
  dev.off()
}

for (gene in aVEC_genes) {
  pdf(file.path(OUTPUT_DIR,
                paste0("ExtData_Fig8a_aVEC_", gene, "_FeaturePlot.pdf")),
      width = 18, height = 6)
  print(FeaturePlot(ISC_annotated_seurat, features = gene,
                    split.by = "Sample_Name", order = TRUE,
                    cols = c("lightgray", "darkred"), pt.size = 0.3) &
          theme(aspect.ratio = 1))
  dev.off()
}

for (gene in damage_resistance_genes) {
  pdf(file.path(OUTPUT_DIR,
                paste0("ExtData_Fig8a_damage_resistance_", gene, "_FeaturePlot.pdf")),
      width = 18, height = 6)
  print(FeaturePlot(ISC_annotated_seurat, features = gene,
                    split.by = "Sample_Name", order = TRUE,
                    cols = c("lightgray", "darkred"), pt.size = 0.3) &
          theme(aspect.ratio = 1))
  dev.off()
}

for (gene in HR_genes) {
  pdf(file.path(OUTPUT_DIR,
                paste0("ExtData_Fig8a_HR_", gene, "_FeaturePlot.pdf")),
      width = 18, height = 6)
  print(FeaturePlot(ISC_annotated_seurat, features = gene,
                    split.by = "Sample_Name", order = TRUE,
                    cols = c("lightgray", "darkred"), pt.size = 0.3) &
          theme(aspect.ratio = 1))
  dev.off()
}

## Ext. Data Fig. 8b ----
## Module score VlnPlots per cluster, split by condition
score_panels <- list(
  rISC       = "rISC_genes1",
  isthmus    = "isthmus_genes1",
  SSC2a      = "SSC2a_genes1",
  SSC2b      = "SSC2b_genes1",
  SSC2c      = "SSC2c_genes1",
  aVEC       = "aVEC_genes1",
  dmg_resist = "damage_resistance1",
  dmg_resp   = "damage_response1",
  HR         = "HR_genes1"
)

for (nm in names(score_panels)) {
  feat <- score_panels[[nm]]
  pdf(file.path(OUTPUT_DIR,
                paste0("ExtData_Fig8b_", nm, "_VlnPlot_by_cluster.pdf")),
      width = 6, height = 6)
  print(VlnPlot(ISC_annotated_seurat, features = feat, pt.size = 0) +
          geom_boxplot(width = 0.5, alpha = 0.5,
                       position = position_dodge(0.9)) +
          labs(title = nm) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 1) + NoLegend())
  dev.off()
}

## Ext. Data Fig. 8c ----
## Module score FeaturePlots (UMAP)
for (nm in names(score_panels)) {
  feat <- score_panels[[nm]]
  pdf(file.path(OUTPUT_DIR,
                paste0("ExtData_Fig8c_", nm, "_FeaturePlot.pdf")),
      width = 5, height = 5)
  print(FeaturePlot(ISC_annotated_seurat, features = feat,
                    order = TRUE, pt.size = 0.2) +
          scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
          ggtitle(nm) + theme(aspect.ratio = 1))
  dev.off()
}

## Ext. Data Fig. 8d ----
## Module score VlnPlots by condition (all cells and stem cells)
make_score_vln_by_condition <- function(obj, feat, title, suffix = "") {
  VlnPlot(obj, features = feat, group.by = "Sample_Name",
          pt.size = 0, alpha = 0.3) +
    geom_boxplot(width = 0.3, alpha = 0.9) +
    scale_fill_manual(values = condition_colors) +
    stat_compare_means(comparisons  = my_comparisons,
                       method       = "wilcox.test",
                       p.adjust.method = "BH",
                       label        = "p.signif",
                       tip.length   = 0.02, size = 4) +
    NoLegend() + theme(aspect.ratio = 1.2) +
    ggtitle(paste0(title, suffix))
}

for (nm in names(score_panels)) {
  feat <- score_panels[[nm]]
  
  ## All cells
  pdf(file.path(OUTPUT_DIR,
                paste0("ExtData_Fig8d_", nm, "_VlnPlot_by_condition_all.pdf")),
      width = 4, height = 5)
  print(make_score_vln_by_condition(ISC_annotated_seurat, feat, nm,
                                    " — all cells"))
  dev.off()
  
  ## Stem cells only
  pdf(file.path(OUTPUT_DIR,
                paste0("ExtData_Fig8d_", nm, "_VlnPlot_by_condition_stem.pdf")),
      width = 4, height = 5)
  print(make_score_vln_by_condition(ISC_annotated_seurat_stem, feat, nm,
                                    " — stem clusters"))
  dev.off()
}