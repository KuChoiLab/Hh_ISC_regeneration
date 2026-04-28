# 0. Setting ----
mem.maxVSize(60000)
set.seed(1)

BASE_PATH  <- "/Users/choi_lab_sk/Dropbox/2026_Hh_ISC_regeneration/scrnaseq"
RDATA_DIR  <- file.path(BASE_PATH, "RData")
OUTPUT_DIR <- file.path(BASE_PATH, "plot", "Fig4")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

condition_colors <- c("ENR" = "#90979F", "SAG" = "#B04D9D", "PGE2" = "#DE9B27")

# 1. Libraries ----
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)


# 2. Load data ----
ISC_annotated_seurat <- readRDS(file.path(RDATA_DIR, "ISC_annotated_seurat.rds"))
ISC_annotated_seurat$idents_clusters <- Idents(ISC_annotated_seurat)
stem_clusters <- c("Canonical stem cell",
                   "Proliferating stem cell (S)",
                   "Proliferating stem cell (G2M)",
                   "Terminally proliferated stem cell",
                   "Absorptive progenitor",
                   "Enteroendocrine progenitor",
                   "Secretory progenitor")

table(Idents(ISC_annotated_seurat))

## Add module scores for Fig. 4l
aVEC_genes <- c("Clu", "Msln", "Sprr1a", "Adam8", "Plaur", "Ly6d", "Plat",
                "Acsbg1", "Mal", "Anxa1", "Isg20", "Anxa5", "Inhba", "Vill",
                "Cd44", "Emp2", "F3", "Hmox1", "Il11", "Krt7")
SSC2a_genes <- c("Hmgb2", "Tuba1b", "Ncl", "Nap1l1", "Hmgb1", "Npm1", "Set",
                 "Dut", "Anp32b", "Hspd1", "Dek", "Ranbp1", "Ptma", "Ifitm3",
                 "Hspe1", "Hsp90aa1", "Lyar", "Ube2c", "Top2a", "Tubb5")
SSC2b_genes <- c("S100a11", "Tmsb4x", "Sprr2a3", "Tmsb10", "Rps26", "Rps20",
                 "Rpl36a", "S100a6", "Eif3e", "Rpl12", "Gabarap", "Rpl10")
SSC2c_genes <- c("Malat1", "Ly6d", "Clu", "Cldn4", "Sprr1a", "Areg", "Lamc2",
                 "Il1rn", "Il33", "Ly6a", "Plaur", "Itgb6", "Suox", "Anxa3")
rISC_genes <- c("Tert", "Bmi1", "Hopx", "Lrig1", "Msi1", "Sox9")
isthmus_genes <- c("Fgfbp1", "Atad2", "Birc5", "Stmn1", "Lgr4", "Mki67")

ISC_annotated_seurat <- AddModuleScore(ISC_annotated_seurat, features = list(aVEC_genes),   name = "aVEC_genes")
ISC_annotated_seurat <- AddModuleScore(ISC_annotated_seurat, features = list(SSC2a_genes),  name = "SSC2a_genes")
ISC_annotated_seurat <- AddModuleScore(ISC_annotated_seurat, features = list(SSC2b_genes),  name = "SSC2b_genes")
ISC_annotated_seurat <- AddModuleScore(ISC_annotated_seurat, features = list(SSC2c_genes),  name = "SSC2c_genes")
ISC_annotated_seurat <- AddModuleScore(ISC_annotated_seurat, features = list(rISC_genes),   name = "rISC_genes")
ISC_annotated_seurat <- AddModuleScore(ISC_annotated_seurat, features = list(isthmus_genes),name = "isthmus_genes")

meta_data <- ISC_annotated_seurat@meta.data

# 3. Figures ----
## Fig. 4c ----
pdf(file.path(OUTPUT_DIR, "Fig4c_UMAP_by_condition.pdf"), width = 5, height = 5)
DimPlot(ISC_annotated_seurat, group.by = "Sample_Name", cols = condition_colors) +
  theme(aspect.ratio = 1) + ggtitle("Fig. 4c")
dev.off()

## Fig. 4d ----
pdf(file.path(OUTPUT_DIR, "Fig4d_UMAP_by_cluster.pdf"), width = 7, height = 6)
DimPlot(ISC_annotated_seurat) + theme(aspect.ratio = 1) + ggtitle("Fig. 4d")
dev.off()

## Fig. 4e ----
distribution_table_cluster <- ISC_annotated_seurat@meta.data %>%                                                                                                                   
  group_by(Sample_Name, idents_clusters) %>%                                                                                                                                       
  dplyr::summarise(count = n(), .groups = "drop") %>%                                                                                                                              
  tidyr::pivot_wider(names_from  = Sample_Name,                                                                                                                                    
                     values_from = count,
                     values_fill = 0)                                                                                                                                              

total_counts <- ISC_annotated_seurat@meta.data %>%                                                                                                                                 
  group_by(Sample_Name) %>%
  dplyr::summarise(total = n(), .groups = "drop") %>%                                                                                                                              
  tibble::deframe()          # named numeric vector       

pairwise_comparisons <- list(                                                                                                                                                      
  c("ENR", "SAG"),                                        
  c("SAG", "PGE2"),                                                                                                                                                                
  c("ENR", "PGE2")
)                                                                                                                                                                                  

pairwise_results_cluster <- lapply(seq_len(nrow(distribution_table_cluster)), function(i) {                                                                                        
  cluster  <- as.character(distribution_table_cluster$idents_clusters[i])
  counts   <- distribution_table_cluster[i, ]                                                                                                                                      
  p_values <- sapply(pairwise_comparisons, function(grps) {
    prop.test(                                                                                                                                                                     
      x = c(as.integer(counts[[grps[1]]]), as.integer(counts[[grps[2]]])),
      n = c(total_counts[grps[1]],         total_counts[grps[2]])                                                                                                                  
    )$p.value                                                                                                                                                                      
  })
  data.frame(Cluster    = cluster,                                                                                                                                                 
             Comparison = c("ENR vs SAG", "SAG vs PGE2", "ENR vs PGE2"),
             p_value    = p_values)                                                                                                                                                
})

pairwise_summary_cluster <- do.call(rbind, pairwise_results_cluster) %>%
  mutate(
    padj_BH      = p.adjust(p_value, method = "BH"),                                                                                                                               
    significance = case_when(
      padj_BH < 0.001 ~ "***",                                                                                                                                                     
      padj_BH < 0.01  ~ "**",                             
      padj_BH < 0.05  ~ "*",                                                                                                                                                       
      TRUE            ~ "ns"
    )                                                                                                                                                                              
  )                                                       

print(pairwise_summary_cluster[, c("Cluster", "Comparison", "padj_BH", "significance")])                                                                                           

cond_order <- c("ENR", "SAG", "PGE2")                     
dodge_w    <- 0.9                                                                                                                                                                  
step       <- dodge_w / length(cond_order) 
offsets    <- setNames(                                                                                                                                                            
  (seq_along(cond_order) - (length(cond_order) + 1) / 2) * step,
  cond_order                                                                                                                                                                       
)   # ENR = -0.3, SAG = 0, PGE2 = +0.3                    

x_pos <- setNames(seq_along(cluster_order), cluster_order)                                                                                                                         

max_pct_by_cluster <- cluster_pct %>%                                                                                                                                              
  group_by(idents_clusters) %>%
  summarise(max_pct = max(pct), .groups = "drop")                                                                                                                                  

sig_brackets <- pairwise_summary_cluster %>%                                                                                                                                       
  filter(significance != "ns") %>%
  mutate(                                                                                                                                                                          
    g1        = sub(" vs .*", "", Comparison),            
    g2        = sub(".* vs ", "", Comparison),
    x_cluster = x_pos[Cluster],
    xmin      = x_cluster + offsets[g1],                                                                                                                                           
    xmax      = x_cluster + offsets[g2]
  ) %>%                                                                                                                                                                            
  left_join(max_pct_by_cluster, by = c("Cluster" = "idents_clusters")) %>%
  arrange(Cluster, xmin) %>%                                                                                                                                                       
  group_by(Cluster) %>%
  mutate(y = max_pct * (1.08 + (row_number() - 1) * 0.12)) %>%                                                                                                                     
  ungroup()                                                                                                                                                                        

tick_h <- diff(range(cluster_pct$pct)) * 0.015

pdf(file.path(OUTPUT_DIR, "Fig4e_cluster_proportion_per_condition.pdf"),
    width = 14, height = 5)                                                                                                                                                        
ggplot(cluster_pct,                                       
       aes(x = factor(idents_clusters, levels = cluster_order),                                                                                                                    
           y = pct, fill = Sample_Name)) +                
  geom_bar(stat = "identity", position = position_dodge(width = dodge_w)) +                                                                                                        
  scale_fill_manual(values = condition_colors) +
  geom_segment(data = sig_brackets,                                                                                                                                                
               aes(x = xmin, xend = xmax, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.35) +                                                                                                                            
  geom_segment(data = sig_brackets,
               aes(x = xmin, xend = xmin, y = y - tick_h, yend = y),                                                                                                               
               inherit.aes = FALSE, linewidth = 0.35) +                                                                                                                            
  geom_segment(data = sig_brackets,                                                                                                                                                
               aes(x = xmax, xend = xmax, y = y - tick_h, yend = y),
               inherit.aes = FALSE, linewidth = 0.35) +                                                                                                                            
  geom_text(data = sig_brackets,                                                                                                                                                   
            aes(x = (xmin + xmax) / 2, y = y + tick_h * 0.5,                                                                                                                       
                label = significance),
            inherit.aes = FALSE, size = 2.2, vjust = 0) +                                                                                                                          
  labs(title = "Fig. 4e", x = NULL, y = "Cell proportion (%)") +                                                                                                                   
  theme_minimal(base_size = 9) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),                                                                                                                        
        legend.title = element_blank())                                                                                                                                            
dev.off()

## Fig. 4f ----
phase_order <- rev(c("qG0", "G1", "Late G1", "S", "S/G2", "G2/M", "M/Early G1"))
phase_colors <- c("qG0"        = "#F8766E",                                                                                                                                        
                  "G1"         = "#D19301",
                  "Late G1"    = "#8AAC00",                                                                                                                                        
                  "S"          = "#01B3F6",               
                  "S/G2"       = "#9B8DFF",                                                                                                                                        
                  "G2/M"       = "#D178FF",               
                  "M/Early G1" = "#FB699F")

phase_pct <- ISC_annotated_seurat@meta.data %>%                                                                                                                                    
  filter(ccAFv2 %in% phase_order) %>% 
  mutate(ccAFv2 = factor(ccAFv2, levels = phase_order)) %>%                                                                                                                        
  group_by(Sample_Name, ccAFv2) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%                                                                                                                              
  group_by(Sample_Name) %>%                               
  mutate(pct = count / sum(count) * 100) %>%                                                                                                                                       
  ungroup()                                                                                                                                                                        

pdf(file.path(OUTPUT_DIR, "Fig4f_ccAFv2_stacked_bar_all_cells.pdf"),                                                                                                               
    width = 4, height = 5)                                
ggplot(phase_pct, aes(x = Sample_Name, y = pct, fill = ccAFv2)) +                                                                                                                  
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = phase_colors,                                                                                                                                         
                    breaks = phase_order,                                                                                                                                          
                    limits = phase_order) +
  scale_y_continuous(breaks = seq(0, 100, 20),
                     limits = c(0, 100)) + 
  labs(title = "Fig. 4f", x = NULL, y = "Cell proportion (%)",                                                                                                                     
       fill = "ccAFv2 phase") +                           
  theme_minimal(base_size = 10) +                                                                                                                                                  
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme(aspect.ratio = 3)
dev.off()  


## Fig. 4g ----
low_prolif_phases     <- c("qG0", "G1", "Late G1")                                                                                                                                 
active_cycling_phases <- c("S", "S/G2", "G2/M", "M/Early G1")                                                                                                                      

sc_cycling <- ISC_annotated_seurat@meta.data %>%          
  mutate(                                                                                                                                                                          
    compartment = ifelse(idents_clusters %in% stem_clusters, "SC", "nonSC"),
    prolif = case_when(                                                                                                                                                            
      ccAFv2 %in% low_prolif_phases     ~ "Low",          
      ccAFv2 %in% active_cycling_phases ~ "Active",                                                                                                                                
      TRUE                              ~ NA_character_
    )                                                                                                                                                                              
  ) %>%                                                   
  filter(!is.na(prolif)) %>%
  group_by(Sample_Name, compartment, prolif) %>%                                                                                                                                   
  dplyr::summarise(count = n(), .groups = "drop") %>%
  group_by(Sample_Name) %>%                                                                                                                                                        
  mutate(pct = count / sum(count) * 100) %>%              
  ungroup() %>%                                                                                                                                                                    
  mutate(
    group = factor(paste(compartment, prolif, sep = "_"),                                                                                                                              
                   levels = c("nonSC_Active", "nonSC_Low",                                                                                                                             
                              "SC_Active",    "SC_Low")),
    Sample_Name = factor(Sample_Name, levels = c("ENR", "SAG", "PGE2"))
  )                                                                                                                                                                                

sc_boundary <- sc_cycling %>%                                                                                                                                                      
  filter(compartment == "SC") %>%                         
  group_by(Sample_Name) %>%
  summarise(y = sum(pct), .groups = "drop") %>%
  mutate(Sample_Name = factor(Sample_Name, levels = c("ENR", "SAG", "PGE2")))                                                                                                      

ref_y <- sc_boundary %>% filter(Sample_Name == "PGE2") %>% pull(y)

group_colors <- c("SC_Low"      = "#9E9E9E",
                  "SC_Active"   = "#43A047",                                                                                                                                       
                  "nonSC_Low"   = "#9E9E9E",              
                  "nonSC_Active"= "#43A047")                                                                                                                                       

pdf(file.path(OUTPUT_DIR, "Fig4g_SC_vs_nonSC_cycling.pdf"), width = 5, height = 6)                                                                                                 
ggplot(sc_cycling, aes(x = Sample_Name, y = pct, fill = group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +                                                                                                                   
  geom_line(data    = sc_boundary,                                                                                                                                                 
            mapping = aes(x = Sample_Name, y = y, group = 1),                                                                                                                      
            inherit.aes = FALSE,                                                                                                                                                   
            linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = group_colors,                                                                                                                                         
                    breaks = c("SC_Low", "SC_Active"),
                    labels = c("qG0 - Late G1", "S - M/Early G1")) +                                                                                                               
  scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 100.5)) +                                                                                                               
  annotate("segment", x = 3.45, xend = 3.45, y = 0,     yend = 100,                                                                                                                
           color = "black", linewidth = 0.4) +                                                                                                                                     
  annotate("segment", x = 3.40, xend = 3.45, y = ref_y, yend = ref_y,                                                                                                              
           color = "black", linewidth = 0.4) +                                                                                                                                     
  annotate("segment", x = 3.40, xend = 3.45, y = 0,     yend = 0,                                                                                                                  
           color = "black", linewidth = 0.4) +                                                                                                                                     
  annotate("segment", x = 3.40, xend = 3.45, y = 100,   yend = 100,                                                                                                                
           color = "black", linewidth = 0.4) +                                                                                                                                     
  annotate("text", x = 3.55, y = ref_y + (100 - ref_y) / 2,                                                                                                                        
           label = "non-SC\ncluster", size = 3, hjust = 0, lineheight = 0.9) +
  annotate("text", x = 3.55, y = ref_y / 2,                                                                                                                                        
           label = "SC\ncluster", size = 3, hjust = 0, lineheight = 0.9) +
  coord_cartesian(clip = "off") +                                                                                                                                                  
  labs(title = "Fig. 4g", x = NULL, y = "Cell proportion (%)", fill = NULL) +
  theme_minimal(base_size = 10) +                                                                                                                                                  
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.margin = margin(5, 70, 5, 5, "pt"))                                                                                                                                   
dev.off() 

## Fig. 4h ----
ISC_annotated_seurat_stem    <- subset(ISC_annotated_seurat, idents %in% stem_clusters)
ISC_annotated_seurat_nonstem <- subset(ISC_annotated_seurat, !idents %in% stem_clusters)

get_sig_labels <- function(seurat_obj, y_pos = 0.72) {                                  
  meta   <- seurat_obj@meta.data                                                                                                                                                   
  scores <- split(meta$CytoTRACE2_Score, meta$Sample_Name)                                                                                                                         
  
  p_vals <- c(                                                                                                                                                                     
    SAG  = wilcox.test(scores[["ENR"]], scores[["SAG"]])$p.value,                                                                                                                  
    PGE2 = wilcox.test(scores[["ENR"]], scores[["PGE2"]])$p.value
  )                                                                                                                                                                                
  padj <- p.adjust(p_vals, method = "BH")
  
  data.frame(                                             
    group = names(padj),
    padj  = padj,                                                                                                                                                                  
    label = case_when(
      padj < 0.001 ~ "***",                                                                                                                                                        
      padj < 0.01  ~ "**",                                
      padj < 0.05  ~ "*",
      TRUE         ~ ""                                                                                                                                                            
    ),
    y = y_pos,                                                                                                                                                                     
    stringsAsFactors = FALSE                              
  ) %>%
    dplyr::filter(label != "")
}

sig_all <- get_sig_labels(ISC_annotated_seurat)                                                                                                                                    
sig_sc  <- get_sig_labels(ISC_annotated_seurat_stem)
sig_nsc <- get_sig_labels(ISC_annotated_seurat_nonstem)                                                                                                                            

make_vln <- function(seurat_obj, title, sig_df) {                                                                                                                                  
  VlnPlot(seurat_obj, features = "CytoTRACE2_Score",                                                                                                                               
          group.by = "Sample_Name", pt.size = 0) +
    geom_boxplot(alpha = 0.5, width = 0.2) +                                                                                                                                       
    scale_fill_manual(values = condition_colors) +                                                                                                                                 
    scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75),
                       limits = c(0, 0.75)) +                                                                                                                                      
    geom_text(data = sig_df,                              
              aes(x = group, y = y, label = label),                                                                                                                                
              inherit.aes = FALSE, size = 4, vjust = 0) + 
    NoLegend() + ggtitle(title) +                                                                                                                                                  
    theme(aspect.ratio = 1.2)
}                                                                                                                                                                                  

p_all <- make_vln(ISC_annotated_seurat,         "All clusters",    sig_all)
p_sc  <- make_vln(ISC_annotated_seurat_stem,    "SC clusters",     sig_sc)                                                                                                         
p_nsc <- make_vln(ISC_annotated_seurat_nonstem, "Non-SC clusters", sig_nsc)                                                                                                        

pdf(file.path(OUTPUT_DIR, "Fig4h_CytoTRACE2_violin.pdf"), width = 10, height = 4)                                                                                                  
p_all | p_sc | p_nsc                                                                                                                                                               
dev.off() 

## Fig. 4i ----

potency_order  <- c("Multipotent", "Oligopotent", "Unipotent", "Differentiated")
potency_colors <- c("Multipotent"    = "#F6836D",                                                                                                                                  
                    "Oligopotent"    = "#BEBBDA",
                    "Unipotent"      = "#FFFFB4",                                                                                                                                  
                    "Differentiated" = "#8CD4C8")   

potency_pct <- ISC_annotated_seurat@meta.data %>%
  dplyr::filter(!is.na(CytoTRACE2_Potency)) %>%
  mutate(CytoTRACE2_Potency = factor(CytoTRACE2_Potency, levels = potency_order)) %>%
  group_by(Sample_Name, CytoTRACE2_Potency) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  group_by(Sample_Name) %>%
  mutate(pct = count / sum(count) * 100) %>%
  ungroup()

pdf(file.path(OUTPUT_DIR, "Fig4i_CytoTRACE2_classification_stacked_bar.pdf"),
    width = 4, height = 5)
ggplot(potency_pct, aes(x = Sample_Name, y = pct, fill = CytoTRACE2_Potency)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = potency_colors, breaks = potency_order) +
  labs(title = "Fig. 4i", x = NULL, y = "Cell proportion (%)",
       fill = "CytoTRACE2 class") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
dev.off()

## Fig. 4j ----
# processed in jupyter notebook 

## Fig. 4k ----
library(viridis)                                                                                                                                                                   

### 1) ENR ----                                                                                                                                                                    
df <- read.table("data/PAGA/dynamic_ENR_circos_input.tsv", header = TRUE, sep = "\t")                                                                                              
df2 <- df                                                                                                                                                                          
df2$source <- c(15, 1, 10, 13, 8, 1, 6, 12, 2, 14, 17, 6, 15, 7, 8, 14)
df2$target <- c(1, 10, 11, 12, 9, 5, 3, 3, 3, 3, 3, 4, 4, 4, 16, 16)                                                                                                               
all_clusters <- sort(unique(c(df2$source, df2$target)))                                                                                                                            
df2$source <- factor(df2$source, levels = all_clusters)                                                                                                                            
df2$target <- factor(df2$target, levels = all_clusters)                                                                                                                            
cluster_colors <- setNames(rainbow(length(all_clusters)), all_clusters)                                                                                                            
df_ENR <- df2                                                                                                                                                                      

### 2) SAG ----                                                                                                                                                                    
df <- read.table("data/PAGA/dynamic_SAG_circos_input.tsv", header = TRUE, sep = "\t")
df2 <- df                                                                                                                                                                          
df2$source <- c(10, 9, 13, 10, 10, 6, 7, 2, 4, 14, 17, 1, 13, 5, 15, 13)
df2$target <- c(11, 11, 12, 8, 5, 3, 3, 3, 3, 3, 3, 2, 16, 16, 14, 14)                                                                                                             
all_clusters <- sort(unique(c(df2$source, df2$target)))                                                                                                                            
df2$source <- factor(df2$source, levels = all_clusters)                                                                                                                            
df2$target <- factor(df2$target, levels = all_clusters)                                                                                                                            
df_SAG <- df2                                                                                                                                                                      

### 3) PGE2 ----                                                                                                                                                                   
df <- read.table("data/PAGA/dynamic_PGE2_circos_input.tsv", header = TRUE, sep = "\t")
df2 <- df                                                                                                                                                                          
df2$source <- c(15, 6, 11, 8, 9, 17, 16, 7, 14, 15, 13, 2, 4, 14, 3, 17)
df2$target <- c(1, 10, 10, 10, 10, 13, 12, 8, 7, 5, 5, 5, 5, 3, 4, 16)                                                                                                             
all_clusters <- sort(unique(c(df2$source, df2$target)))                                                                                                                            
df2$source <- factor(df2$source, levels = all_clusters)                                                                                                                            
df2$target <- factor(df2$target, levels = all_clusters)                                                                                                                            
df_PGE2 <- df2                                                                                                                                                                     

pdf(file = file.path(OUTPUT_DIR, "Fig4k_circos_PAGA.pdf"), width = 18, height = 6)
par(mfrow = c(1, 3))                                                                                                                                                               

set.seed(1)                                                                                                                                                                        
chordDiagram(df_ENR,
             annotationTrack = c('name', 'grid', 'axis'),                                                                                                                          
             annotationTrackHeight = c(0.05, 0.1, 0.01),  
             grid.col = cluster_colors,                                                                                                                                            
             directional = 1,
             direction.type = c('diffHeight', 'arrows'))                                                                                                                           
title("ENR", line = -2)                                   
circos.clear()

set.seed(1)
chordDiagram(df_SAG,                                                                                                                                                               
             annotationTrack = c('name', 'grid', 'axis'), 
             annotationTrackHeight = c(0.05, 0.1, 0.01),
             grid.col = cluster_colors,                                                                                                                                            
             directional = 1,
             direction.type = c('diffHeight', 'arrows'))                                                                                                                           
title("SAG", line = -2)                                   
circos.clear()

set.seed(1)
chordDiagram(df_PGE2,                                                                                                                                                              
             annotationTrack = c('name', 'grid', 'axis'), 
             annotationTrackHeight = c(0.05, 0.1, 0.01),
             grid.col = cluster_colors,                                                                                                                                            
             directional = 1,
             direction.type = c('diffHeight', 'arrows'))                                                                                                                           
title("PGE2", line = -2)                                  
circos.clear()                                                                                                                                                                     

dev.off() 

## Fig. 4l ----

fig4l_genes <- c(                         
  "Hopx", "Bmi1", "Lrig1", "Msi1", "Sox9", "Tert",                  # rISCs             
  "Fgfbp1", "Atad2", "Lgr4", "Birc5", "Stmn1",                      # Isthmus                                                                                                      
  "Ly6a", "Clu1",                                                     # revSC                                                                                                      
  "Ly6a", "Clu1", "Itgb6", "Plaur", "Lamc2", "Il33",                 # aVEC                                                                                                        
  "Sprr1a", "Ly6d", "Suox", "Anxa3", "Krt7",                         # Fetal-like                                                                                                  
  "Atm", "Bax", "Chek1", "Chek2", "Ddb2", "Noxa1", "Trp53", "Xpc",  # Damage resistance                                                                                            
  "Brca1", "Brca2", "Rad51"                                           # HR                                                                                                         
) 

fig4l_genes <- unique(fig4l_genes)

ISC_annotated_seurat_cond <- ISC_annotated_seurat
Idents(ISC_annotated_seurat_cond) <- factor(ISC_annotated_seurat_cond$Sample_Name,
                                   levels = c("ENR", "SAG", "PGE2"))

pdf(file.path(OUTPUT_DIR, "Fig4l_DotPlot_gene_sets_by_condition.pdf"),
    width = 14, height = 4)
DotPlot(ISC_annotated_seurat_cond, features = fig4l_genes) +
  scale_colour_gradientn(colors = rev(brewer.pal(11, "RdBu"))) +
  scale_size(range = c(1, 6)) +
  RotatedAxis() +
  labs(title = "Fig. 4l") +
  theme(axis.text.x = element_text(size = 8))
dev.off()


## Fig. 4m ----
library(ggpubr)

mature_celltypes <- c("Enterocyte", "Secretory cell", "Enteroendocrine", "Tuft cell")                                                                                              

make_correlation_plot <- function(obj, cell_type) {                                                                                                                                
  cells <- rownames(obj@meta.data)[obj@meta.data$idents_clusters == cell_type]
  if (length(cells) < 10) {                                                                                                                                                        
    message("Skipping ", cell_type, " — insufficient cells.")
    return(NULL)                                                                                                                                                                   
  }                                                       
  sub <- obj[, cells]                                                                                                                                                              
  Idents(sub) <- sub$Sample_Name
  
  ENR_cells  <- WhichCells(sub, idents = "ENR")           
  SAG_cells  <- WhichCells(sub, idents = "SAG")
  PGE2_cells <- WhichCells(sub, idents = "PGE2")                                                                                                                                   
  
  sd_data  <- GetAssayData(sub, layer = "scale.data")                                                                                                                              
  ENR_avg  <- rowMeans(sd_data[, ENR_cells,  drop = FALSE])
  SAG_avg  <- rowMeans(sd_data[, SAG_cells,  drop = FALSE])                                                                                                                        
  PGE2_avg <- rowMeans(sd_data[, PGE2_cells, drop = FALSE])
  
  r_sag <- round(cor(ENR_avg, SAG_avg,  method = "pearson"), 2)                                                                                                                    
  r_pge <- round(cor(ENR_avg, PGE2_avg, method = "pearson"), 2)                                                                                                                    
  
  df_combined <- rbind(                                   
    data.frame(gene = names(ENR_avg), ENR = ENR_avg, Other = SAG_avg,
               comparison = "ENR vs SAG"),                                                                                                                                         
    data.frame(gene = names(ENR_avg), ENR = ENR_avg, Other = PGE2_avg,
               comparison = "ENR vs PGE2")                                                                                                                                         
  ) %>%                                                   
    mutate(comparison = factor(comparison, levels = c("ENR vs SAG", "ENR vs PGE2")))                                                                                               
  
  ggplot(df_combined, aes(x = ENR, y = Other, color = comparison)) +                                                                                                               
    geom_point(alpha = 0.4, size = 0.5) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +                                                                                                                       
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +                                                                                                 
    scale_color_manual(
      values = c("ENR vs SAG"  = unname(condition_colors["SAG"]),                                                                                                                  
                 "ENR vs PGE2" = unname(condition_colors["PGE2"])),
      labels = c("ENR vs SAG"  = sprintf("ENR vs SAG (r = %.2f)",  r_sag),                                                                                                         
                 "ENR vs PGE2" = sprintf("ENR vs PGE2 (r = %.2f)", r_pge))                                                                                                         
    ) +                                                                                                                                                                            
    coord_fixed() +                                                                                                                                                                
    theme_minimal(base_size = 9) +                                                                                                                                                 
    theme(aspect.ratio    = 1,                            
          legend.position = "bottom",
          legend.title    = element_blank()) +                                                                                                                                     
    ggtitle(cell_type) +
    xlab("ENR avg expression (scaled)") +                                                                                                                                          
    ylab("Other condition avg expression (scaled)")       
}                                                                                                                                                                                  

plot_list <- Filter(Negate(is.null),                                                                                                                                               
                    setNames(lapply(mature_celltypes, make_correlation_plot,
                                    obj = ISC_annotated_seurat),
                             mature_celltypes))

pdf(file.path(OUTPUT_DIR, "Fig4m_pseudobulk_correlation.pdf"), width = 12, height = 4)
wrap_plots(plot_list, ncol = 4)                                                                                                                                                    
dev.off()  
