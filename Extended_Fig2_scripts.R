# 0. setting ----
BASE_PATH  <- "/Users/choi_lab_sk/Dropbox/2026_Hh_ISC_regeneration/bulkrnaseq"
PREP_DIR   <- file.path(BASE_PATH, "data", "preprocessing")
OUTPUT_DIR <- file.path(BASE_PATH, "plot", "Extended_Data_Fig2")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

GOREA_SCRIPT <- "~/Library/CloudStorage/OneDrive-Personal/Choi_Lab_SK/sngsoo_R/20250401_GOrea_HJ/mouse_7.5.1/20250401_gorea_function_mouse_v7.5.1_hj_pdf.R"
GOREA_DIR    <- "~/Library/CloudStorage/OneDrive-Personal/Choi_Lab_SK/sngsoo_R/20250401_GOrea_HJ/"

palette_colors_list <- list(
  sampletype = c("Homeostasis"          = "#008d42",
                 "Lgr5_DTA"             = "#e80072",
                 "Lgr5_DTA_Cyclopamine" = "#ba7c37")
)

# 1. Libraries ----
library(DESeq2)
library(apeglm)
library(fgsea)
library(msigdbr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(patchwork)
library(EnhancedVolcano)
library(GO.db)
library(GOSemSim)                                                                                                                                                                  
library(org.Mm.eg.db)
library(simplifyEnrichment)
library(ComplexHeatmap)
library(WriteXLS)
library(colorRamp2)

# 2. Load preprocessed objects ----
dds_wald    <- readRDS(file.path(PREP_DIR, "dds_wald.rds"))
dds_relevel <- readRDS(file.path(PREP_DIR, "dds_relevel.rds"))  # DTA as reference

gene_lists     <- readRDS(file.path(PREP_DIR, "gsea_gene_lists.rds"))
gene_list_d_h  <- gene_lists$d_h   # Wald stat: DTA vs Homeostasis
gene_list_dc_h <- gene_lists$dc_h  # Wald stat: DC vs Homeostasis
gene_list_dc_d <- gene_lists$dc_d  # Wald stat: Lgr5_DTA_Cyclopamine vs DTA
gene_list_h_d <- gene_lists$h_d  # Wald stat: Homeostasis vs DTA
gene_list_d_dc <- gene_lists$d_dc  # Wald stat: DTA vs Lgr5_DTA_Cyclopamine

# 3. MSigDB pathway collections ----
pathways_H <- split(
  msigdbr(species = "Mus musculus", collection = "H")$gene_symbol,
  msigdbr(species = "Mus musculus", collection = "H")$gs_name)

pathways_GOBP <- split(
  msigdbr(species = "Mus musculus", collection = "C5",
          subcollection = "GO:BP")$gene_symbol,
  msigdbr(species = "Mus musculus", collection = "C5",
          subcollection = "GO:BP")$gs_name)

# 4. Figures ----
## Ext. Data Fig. 2a ----
# Convert DESeq2 results to data.frame with symbol column
res_d_h_tb  <- as.data.frame(res_d_h_shrunk)  %>%
  rownames_to_column("symbol")
res_d_dc_tb <- as.data.frame(res_d_dc_shrunk) %>%
  rownames_to_column("symbol")

p.cut <- 0.05
f.cut <- 1

Hh_activation_gene_set <- readRDS('RData/Hh_activation_gene_set.rds')

## Left: Lgr5-DTA vs Lgr5-GFP (Homeostasis)
pdf(file.path(OUTPUT_DIR, "ExtData_Fig2a_volcano_DTA_vs_GFP.pdf"),
    width = 10, height = 10)
EnhancedVolcano(res_d_h_tb,
                lab            = res_d_h_tb$symbol,
                selectLab      = Hh_activation_gene_set,
                boxedLabels    = TRUE,
                drawConnectors = TRUE,
                max.overlaps   = Inf,
                x              = "log2FoldChange",
                y              = "padj",
                pCutoff        = p.cut,
                pCutoffCol     = "padj",
                FCcutoff       = f.cut,
                labSize        = 3,
                legendPosition = "top",
                subtitle = bquote(paste("p.adjust < ", .(p.cut),
                                        " & |", "Log"[2], " FC| ≥ ", .(f.cut))),
                title = "Lgr5-DTA (TMX) vs Lgr5-GFP (TMX)") +
  theme(aspect.ratio = 1)
dev.off()

## Right: Lgr5-DTA vs Lgr5-DTA+Cyclopamine
pdf(file.path(OUTPUT_DIR, "ExtData_Fig2a_volcano_DTA_vs_DTA_CYC.pdf"),
    width = 10, height = 10)
EnhancedVolcano(res_d_dc_tb,
                lab            = res_d_dc_tb$symbol,
                selectLab      = Hh_activation_gene_set,
                boxedLabels    = TRUE,
                drawConnectors = TRUE,
                max.overlaps   = Inf,
                x              = "log2FoldChange",
                y              = "padj",
                pCutoff        = p.cut,
                pCutoffCol     = "padj",
                FCcutoff       = f.cut,
                labSize        = 3,
                legendPosition = "top",
                subtitle = bquote(paste("p.adjust < ", .(p.cut),
                                        " & |", "Log"[2], " FC| ≥ ", .(f.cut))),
                title = "Lgr5-DTA (TMX) vs Lgr5-DTA (TMX + CYC)") +
  theme(aspect.ratio = 1)
dev.off()

## Ext. Data Fig. 2b ----
pid <- "REACTOME_HEDGEHOG_ON_STATE"

pathways_C2_CP_REACTOME_final <- split(
  msigdbr::msigdbr(species = "Mus musculus", collection = "C2",
                   subcollection = "CP:REACTOME")$gene_symbol,
  msigdbr::msigdbr(species = "Mus musculus", collection = "C2",
                   subcollection = "CP:REACTOME")$gs_name)

fgsea_RE_d_h  <- fgseaMultilevel(pathways = pathways_C2_CP_REACTOME_final,
                                 stats    = gene_list_d_h)
fgsea_RE_d_dc <- fgseaMultilevel(pathways = pathways_C2_CP_REACTOME_final,
                                 stats    = gene_list_d_dc)

## Left: Lgr5-DTA vs Lgr5-GFP
pdf(file.path(OUTPUT_DIR, "ExtData_Fig2b_Hh_enrichment_DTA_vs_GFP.pdf"))
plotEnrichment(pathways_C2_CP_REACTOME_final[[pid]], gene_list_d_h) +
  ggtitle(sprintf("%s\nNES = %.2f,  padj = %.3g",
                  pid,
                  subset(fgsea_RE_d_h,  pathway == pid)$NES,
                  subset(fgsea_RE_d_h,  pathway == pid)$padj)) +
  theme(aspect.ratio = 1)
dev.off()

## Right: Lgr5-DTA vs Lgr5-DTA+CYC
pdf(file.path(OUTPUT_DIR, "ExtData_Fig2b_Hh_enrichment_DTA_vs_DTA_CYC.pdf"))
plotEnrichment(pathways_C2_CP_REACTOME_final[[pid]], gene_list_d_dc) +
  ggtitle(sprintf("%s\nNES = %.2f,  padj = %.3g",
                  pid,
                  subset(fgsea_RE_d_dc, pathway == pid)$NES,
                  subset(fgsea_RE_d_dc, pathway == pid)$padj)) +
  theme(aspect.ratio = 1)
dev.off()

## Ext. Data Fig. 2c ----
source(GOREA_SCRIPT)
gorea_enviromnet(GOREA_DIR)

## fgsea on GOBP for all 4 directions
run_fgsea_gobp <- function(gene_list) {
  res <- fgseaMultilevel(pathways = pathways_GOBP,
                         stats    = gene_list,
                         maxSize  = 500)
  res <- as.data.frame(res)[!is.na(res$padj) & res$padj < 0.05, ]
  res[order(res$padj), ]
}

message("Running GOBP fgsea — 4 comparisons (this may take ~10 min each)...")
fgsea_h_d  <- run_fgsea_gobp(gene_list_h_d)   # Homeostasis up vs DTA → top-left
fgsea_d_h  <- run_fgsea_gobp(gene_list_d_h)   # DTA up vs H           → top-right
fgsea_dc_d <- run_fgsea_gobp(gene_list_dc_d)  # CYC up vs DTA         → bottom-left
fgsea_d_dc <- run_fgsea_gobp(gene_list_d_dc)  # DTA up vs CYC         → bottom-right

## GOREA runner
run_gorea_panel <- function(fgsea_sig, out_dir,
                            k_val = 20,
                            heatmap_width = 40, heatmap_height = 30) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  old_wd <- setwd(out_dir)
  on.exit(setwd(old_wd), add = TRUE)
  
  GOID     <- GOID_TERM[GOID_TERM$GOTERM %in% fgsea_sig$pathway, ]
  merged   <- base::merge(GOID, fgsea_sig,
                          by.x = "GOTERM", by.y = "pathway", all = TRUE)
  input_df <- merged %>%
    dplyr::filter(NES > 0) %>%
    dplyr::select(GOID, NES)
  
  if (nrow(input_df) == 0) {
    message("  No positive-NES terms — skipping GOREA for ", out_dir)
    return(invisible(NULL))
  }
  
  w <- gorea_sim_mat(input = input_df, godata_GO = godata_GO)
  gorea_outlier_plot(w = w)
  
  gorea(input                          = input_df,
        k_val                          = k_val,
        godata_GO                      = godata_GO,
        cutoff                         = 0.85,
        outlier_detect                 = TRUE,
        min_cluster                    = 3,
        representative_term_level_cutoff = 1,
        GO_explain                     = 3,
        score                          = "NES",
        filename1                      = "gorea_goterm.xlsx",
        filename2                      = "gorea_ancestor.xlsx",
        heatmap_filename               = "gorea_plot.pdf",
        plot                           = TRUE,
        heatmap_width                  = heatmap_width,
        heatmap_height                 = heatmap_height,
        ancestor_annotation            = TRUE,
        right_annotation_font_size     = 10,
        cluster_font_size              = 4,
        top_ancestor_annotation        = TRUE,
        top_ancestor                   = 3,
        color                          = c("gold", "white"))
}

## Top-left: Lgr5-GFP (Homeostasis) up vs Lgr5-DTA
message("GOREA 1/4: Lgr5-GFP vs Lgr5-DTA (up in Homeostasis)")
run_gorea_panel(fgsea_h_d,
                out_dir = file.path(OUTPUT_DIR, "GOREA_1_GFP_vs_DTA"),
                k_val   = 20)

## Top-right: Lgr5-DTA up vs Lgr5-GFP
message("GOREA 2/4: Lgr5-DTA vs Lgr5-GFP (up in DTA)")
run_gorea_panel(fgsea_d_h,
                out_dir = file.path(OUTPUT_DIR, "GOREA_2_DTA_vs_GFP"),
                k_val   = 20)

## Bottom-left: Lgr5-DTA+CYC up vs Lgr5-DTA
message("GOREA 3/4: Lgr5-DTA+CYC vs Lgr5-DTA (up in Cyclopamine)")
run_gorea_panel(fgsea_dc_d,
                out_dir = file.path(OUTPUT_DIR, "GOREA_3_DC_vs_DTA"),
                k_val   = 20)

## Bottom-right: Lgr5-DTA up vs Lgr5-DTA+CYC
message("GOREA 4/4: Lgr5-DTA vs Lgr5-DTA+CYC (up in DTA)")
run_gorea_panel(fgsea_d_dc,
                out_dir = file.path(OUTPUT_DIR, "GOREA_4_DTA_vs_DC"),
                k_val   = 20)