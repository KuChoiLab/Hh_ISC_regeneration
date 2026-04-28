# 0. setting ----
BASE_PATH  <- "/Users/choi_lab_sk/Dropbox/2026_Hh_ISC_regeneration/bulkrnaseq"
PREP_DIR    <- file.path(BASE_PATH, "data", "preprocessing")
OUTPUT_DIR  <- file.path(BASE_PATH, "plot", "Fig1")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

palette_colors_list <- list(
  sampletype = c("Homeostasis"          = "#008d42",
                 "Lgr5_DTA"             = "#e80072",
                 "Lgr5_DTA_Cyclopamine" = "#ba7c37"))

# 1. Libraries ----
library(DESeq2)
library(DEGreport)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(fgsea)
library(msigdbr)
library(ComplexHeatmap)
library(circlize)

# 2. Load preprocessed objects ----
dds_wald     <- readRDS(file.path(PREP_DIR, "dds_wald.rds"))
dds_lrt      <- readRDS(file.path(PREP_DIR, "dds_lrt.rds"))
rld          <- readRDS(file.path(PREP_DIR, "rld.rds"))
gene_lists   <- readRDS(file.path(PREP_DIR, "gsea_gene_lists.rds"))
gene_list_d_h  <- gene_lists$d_h
gene_list_dc_h <- gene_lists$dc_h
gene_list_dc_d <- gene_lists$dc_d
gene_list_d_dc <- gene_lists$d_dc

tx2gene <- read.delim("/Users/Choi_Lab_sk/Library/CloudStorage/OneDrive-Personal/Choi_Lab_SK/project1_intestinal_crypt_bulk_scRNA-seq/bulk_ver2/metadata/tx2gene.gencode.vM34.csv")
colnames(tx2gene) <- c("tx_id","ensgene","symbol")
vM34annot <- tx2gene %>% dplyr::select(ensgene, symbol) %>% dplyr::distinct()

norm_counts <- as.data.frame(counts(dds_wald, normalized = TRUE))
norm_counts$symbol <- rownames(norm_counts)
rld_mat <- assay(rld)

meta <- data.frame(sampletype = colData(dds_wald)$sampletype,
                   row.names  = colnames(dds_wald))

# 3. Figures ----
## Fig. 1e ----
fig1e_genes <- readRDS('RData/fig1e_genes.rds')

rld_df <- as.data.frame(rld_mat)                                                                                                                                                   
rld_df$symbol <- rownames(rld_df)                                                                                                                                                  
fig1e_mat <- rld_df %>%                                                                                                                                                            
  dplyr::filter(symbol %in% fig1e_genes) %>%              
  dplyr::select(-symbol) %>%                                                                                                                                                       
  as.matrix()

cond_order <- c("Homeostasis", "Lgr5_DTA", "Lgr5_DTA_Cyclopamine")                                                                                                                 
col_idx    <- order(match(meta$sampletype, cond_order))                                                                                                                            

fig1e_mat_sorted <- fig1e_mat[, col_idx, drop = FALSE]                                                                                                                             
meta_sorted      <- meta[col_idx, , drop = FALSE]         

n_per_cond <- table(factor(meta_sorted$sampletype, levels = cond_order))

mat_ordered <- fig1e_mat_sorted[fig1e_genes, ]           

mat_scaled <- scale(t(mat_ordered))

pheatmap_colors <- colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(100)
col_fun <- colorRamp2(                                    
  seq(-2, 2, length.out = 100),                                                                                                                                                    
  pheatmap_colors)

palette_row <- c("Homeostasis"          = "#008d42",      
                 "Lgr5_DTA"             = "#e80072",                                                                                                                               
                 "Lgr5_DTA_Cyclopamine" = "#ba7c37")      

ha_left <- rowAnnotation(                                                                                                                                                          
  Condition = meta_sorted$sampletype,
  col       = list(Condition = palette_row),                                                                                                                                       
  show_annotation_name = FALSE,                           
  annotation_legend_param = list(title = "sampletype"))                                                                                                                                                                                  

row_split <- factor(meta_sorted$sampletype,               
                    levels = c("Homeostasis", "Lgr5_DTA", "Lgr5_DTA_Cyclopamine"))                                                                                                

pdf(file.path(OUTPUT_DIR, "Fig1e_Hh_regen_heatmap.pdf"), width = 7, height = 4)                                                                                                    
ht <- Heatmap(                                            
  mat_scaled,           
  col = col_fun,
  cluster_rows     = FALSE,                                                                                                                                                        
  cluster_columns  = FALSE,                               
  row_gap          = unit(2, "mm"),                                                                                                                                                
  row_title        = NULL,
  left_annotation  = ha_left,                                                                                                                                                      
  show_row_names   = FALSE,                               
  column_names_rot = 45,
  name             = "z-score",                           
  heatmap_legend_param = list(                                                                                                                                                     
    title    = "z-score",                                 
    at       = c(-2, -1, 0, 1, 2),                                                                                                                                                 
    direction = "vertical"))                                                         
draw(ht, heatmap_legend_side = "right")
dev.off()

## Fig. 1f ----
go_bp <- msigdbr(species = "Mus musculus", collection = "C5",
                 subcollection = "GO:BP")
pathways_gobp <- split(go_bp$gene_symbol, go_bp$gs_name)

set.seed(1)
fgsea_gobp_d_h <- fgseaMultilevel(pathways = pathways_gobp,
                                  stats    = gene_list_d_h,
                                  maxSize  = 500,
                                  nPermSimple = 10000)

sig_gobp <- fgsea_gobp_d_h %>%
  as_tibble() %>%
  dplyr::filter(padj < 0.05) %>%
  arrange(NES) %>%
  mutate(direction = ifelse(NES > 0, "NES > 0", "NES < 0"),
         short_name = gsub("^GOBP_", "", pathway) %>%
           gsub("_", " ", .) %>%
           tolower() %>%
           tools::toTitleCase())

top_sig <- bind_rows(
  sig_gobp %>% dplyr::filter(NES > 0) %>% arrange(desc(NES)) %>% head(10),
  sig_gobp %>% dplyr::filter(NES < 0) %>% arrange(NES)       %>% head(10)
) %>%
  arrange(NES)

pdf(file.path(OUTPUT_DIR, "Fig1f_GSEA_GOBP_NES_bar.pdf"), width = 9, height = 7)
ggplot(top_sig, aes(x = reorder(short_name, NES), y = NES,
                    fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("NES > 0"   = "#e80072",
                               "NES < 0" = "#008d42")) +
  labs(x = NULL, y = "Normalized Enrichment Score",
       fill  = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", aspect.ratio = 1.2)
dev.off()

## Fig. 1g ----
res_lrt <- results(dds_lrt)
res_lrt_tb <- res_lrt %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()
res_lrt_tb

res_lrt_tb <- merge(res_lrt_tb, vM34annot, by.x="gene", by.y="symbol")
siglrt_genes <- res_lrt_tb %>% dplyr::filter(padj < 0.05)

nrow(siglrt_genes) # 6902

library(DEGreport)
cluster_rlog <- rld_mat[siglrt_genes$gene, ]
clusters_lrt <- degPatterns(cluster_rlog, metadata=meta, time = "sampletype", col=NULL)

pdf(file.path(OUTPUT_DIR, "Fig1g_lrt_degPatterns.pdf"), width = 10, height = 8)
clusters_lrt$plot + facet_wrap(~ cluster, ncol = 4) + theme(legend.position = 'none')
dev.off()

## Save cluster gene lists
lrt_df <- clusters_lrt$df
for (grp in unique(lrt_df$cluster)) {
  genes_grp <- lrt_df %>% dplyr::filter(cluster == grp) %>% pull(genes)
  write.csv(data.frame(gene = genes_grp),
            file.path(OUTPUT_DIR,
                      paste0("Fig1g_lrt_group", grp, "_genes.csv")),
            row.names = FALSE)
}

## Fig. 1h ----
injury_signature <- readRDS('RData/injury_signature.rds')

res_inj_d_h <- fgseaSimple(pathways = injury_signature,
                           stats    = gene_list_d_h,
                           nperm    = 150000)

pdf(file.path(OUTPUT_DIR,
              "Fig1h_left_injury_signature_Lgr5DTA_vs_Homeostasis.pdf"),
    width = 4.5, height = 4)
plotEnrichment(injury_signature[[1]], gene_list_d_h) +
  labs(title = "Injury-associated regenerative signature",
       subtitle = sprintf("Lgr5-DTA vs. Homeostasis\nNES = %.2f, padj = %.3g",
                          res_inj_d_h$NES, res_inj_d_h$padj)) +
  theme_minimal(base_size = 10) + theme(aspect.ratio = 1)
dev.off()

## Right panel: Lgr5-DTA (TMX) vs. Lgr5-DTA+CYC (TMX)
res_inj_d_dc <- fgseaSimple(pathways = injury_signature,
                            stats    = gene_list_d_dc,
                            nperm    = 150000)

pdf(file.path(OUTPUT_DIR,
              "Fig1h_right_injury_signature_Lgr5DTA_vs_Lgr5DTA_CYC.pdf"),
    width = 4.5, height = 4)
plotEnrichment(injury_signature[[1]], gene_list_d_dc) +
  labs(title = "Injury-associated regenerative signature",
       subtitle = sprintf("Lgr5-DTA vs. Lgr5-DTA+Cyclopamine\nNES = %.2f, padj = %.3g",
                          res_inj_d_dc$NES, res_inj_d_dc$padj)) +
  theme_minimal(base_size = 10) + theme(aspect.ratio = 1)
dev.off()