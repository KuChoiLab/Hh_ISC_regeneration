# 0. Setting ----
mem.maxVSize(60000)
set.seed(1)

BASE_PATH  <- "/Users/choi_lab_sk/Dropbox/2026_Hh_ISC_regeneration/scrnaseq"
RDATA_DIR  <- file.path(BASE_PATH, "RData")
OUTPUT_DIR <- file.path(BASE_PATH, "data", "preprocessing")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# 1. Libraries ----
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(harmony)
library(scDblFinder)
library(SingleCellExperiment)
library(ccAFv2)
library(RColorBrewer)

# 2. Load BD Rhapsody processed objects ----
# batch 1
ISC_seurat_batch1 <- readRDS(file.path(BASE_PATH, "RData/ISC_seurat_batch1.rds"))
batch1_cells_to_filter <- readRDS(file.path(BASE_PATH, "RData/batch1_cells_to_filter.rds"))
ISC_seurat_batch1_filtered <- subset(ISC_seurat_batch1, cells = batch1_cells_to_filter, invert=TRUE)

# batch 2
ISC_seurat_batch2 <- readRDS(file.path(BASE_PATH, "RData/ISC_seurat_batch2.rds"))

# merge() function
ISC_seurat_batch1_filtered$orig.ident <- "batch1"
ISC_seurat_batch2$orig.ident <- "batch2"
merged_seurat <- merge(x= ISC_seurat_batch1_filtered,
                       y= ISC_seurat_batch2,
                       add.cell.ids = c("batch1","batch2"),
                       project = "ISC")

# 3. QC metrics ----
merged_seurat$mitoRatio <- PercentageFeatureSet(merged_seurat,
                                                pattern = "^mt-") / 100

## Filter cells
filtered_seurat <- subset(merged_seurat,
                          subset = nFeature_RNA  > 500  &
                            mitoRatio     < 0.5)

## Remove Multiplet and Undetermined barcodes
filtered_seurat <- subset(filtered_seurat,
                          subset = Sample_Name %in% c("ENR", "SAG", "PGE"))
filtered_seurat$Sample_Name <- gsub("PGE$", "PGE2",
                                    filtered_seurat$Sample_Name)
filtered_seurat$Sample_Name <- factor(filtered_seurat$Sample_Name,
                                      levels = c("ENR", "SAG", "PGE2"))

# 4. Doublet detection (scDblFinder) ----
table(filtered_seurat$orig.ident)
# batch1 batch2 
# 2955   6424 

set.seed(123456)

split_seurat <- SplitObject(filtered_seurat, split.by = "Sample_Name")
seurat_ENR <- split_seurat[[1]]
seurat_SAG <- split_seurat[[2]]
seurat_PGE <- split_seurat[[3]]

sce_ENR_2 <- as.SingleCellExperiment(seurat_ENR)
sce_ENR_2 <- scDblFinder(sce_ENR_2, dbr = 0.01, clusters = FALSE)
seurat_ENR$scDblFinder.class <- sce_ENR_2$scDblFinder.class
seurat_ENR$scDblFinder.score <- sce_ENR_2$scDblFinder.score

sce_SAG_2 <- as.SingleCellExperiment(seurat_SAG)
sce_SAG_2 <- scDblFinder(sce_SAG_2, dbr = 0.0085, clusters = FALSE)
seurat_SAG$scDblFinder.class <- sce_SAG_2$scDblFinder.class
seurat_SAG$scDblFinder.score <- sce_SAG_2$scDblFinder.score

sce_PGE_2 <- as.SingleCellExperiment(seurat_PGE)
sce_PGE_2 <- scDblFinder(sce_PGE_2, dbr = 0.005, clusters = FALSE)
seurat_PGE$scDblFinder.class <- sce_PGE_2$scDblFinder.class
seurat_PGE$scDblFinder.score <- sce_PGE_2$scDblFinder.score

scdblfinder_ENR <- rownames(seurat_ENR@meta.data)[seurat_ENR$scDblFinder.class == 'doublet']
scdblfinder_SAG <- rownames(seurat_SAG@meta.data)[seurat_SAG$scDblFinder.class == 'doublet']
scdblfinder_PGE <- rownames(seurat_PGE@meta.data)[seurat_PGE$scDblFinder.class == 'doublet']

scdblfinder_all <- c(scdblfinder_ENR,scdblfinder_SAG,scdblfinder_PGE)
cells <- list('Doublet cell' = scdblfinder_all)

ENR_singlet <- subset(seurat_ENR, cells = scdblfinder_ENR, invert = T)
SAG_singlet <- subset(seurat_SAG, cells = scdblfinder_SAG, invert = T)
PGE_singlet <- subset(seurat_PGE, cells = scdblfinder_PGE, invert = T)

split_seurat_without_doublet <- list()
split_seurat_without_doublet$ENR <- ENR_singlet
split_seurat_without_doublet$SAG <- SAG_singlet
split_seurat_without_doublet$PGE <- PGE_singlet

merged_seurat_without_doublet <- merge(x= split_seurat_without_doublet[[1]],
                                       y = split_seurat_without_doublet[2:3],
                                       merge.data = TRUE)

# 5. Batch correction (Harmony) ----
merged_seurat_without_doublet$batch <- merged_seurat_without_doublet$orig.ident

merged_seurat_without_doublet <- NormalizeData(merged_seurat_without_doublet)
merged_seurat_without_doublet <- FindVariableFeatures(merged_seurat_without_doublet)
merged_seurat_without_doublet <- ScaleData(merged_seurat_without_doublet)
merged_seurat_without_doublet <- RunPCA(merged_seurat_without_doublet)
merged_seurat_without_doublet <- RunUMAP(merged_seurat_without_doublet, dims=1:41)
merged_seurat_without_doublet <- FindNeighbors(merged_seurat_without_doublet, dims = 1:41)
merged_seurat_without_doublet <- FindClusters(merged_seurat_without_doublet)

harmonized_seurat <- merged_seurat_without_doublet
VariableFeatures(harmonized_seurat) <- split(row.names(harmonized_seurat@meta.data), harmonized_seurat@meta.data$batch) %>% lapply(function(cells_use){
  harmonized_seurat[,cells_use] %>% FindVariableFeatures() %>% VariableFeatures()
}) %>% unlist() %>% unique()

harmonized_seurat <- RunHarmony(harmonized_seurat, group.by.vars = 'batch')
harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = 'harmony', dims = 1:41)
harmonized_seurat <- FindNeighbors(harmonized_seurat, reudction = 'harmony', dims = 1:41)
harmonized_seurat <- FindClusters(harmonized_seurat)
DimPlot(harmonized_seurat) + theme(aspect.ratio = 1)

# 6. UMAP & clustering ----
harmonized_seurat <- RunUMAP(harmonized_seurat,
                           reduction = "harmony",
                           dims      = 1:41)
harmonized_seurat <- FindNeighbors(harmonized_seurat,
                                 reduction = "harmony",
                                 dims      = 1:41)
harmonized_seurat <- FindClusters(harmonized_seurat)

DimPlot(harmonized_seurat, label = TRUE) + theme(aspect.ratio = 1)
DimPlot(harmonized_seurat, split.by = 'orig.ident') & theme(aspect.ratio = 1)

# 7. Annotation ----
ISC_processed_seurat <- readRDS("RData/ISC_processed_seurat.rds")
DimPlot(ISC_processed_seurat, label = TRUE) + theme(aspect.ratio = 1)

ISC_annotated_seurat <- RenameIdents(ISC_processed_seurat,
                                     "stem_0" = "Canonical stem cell",                                                     
                                     "stem_1" = "Proliferating stem cell (G2M)",
                                     "stem_2" = "Secretory progenitor",                                                                                                                                               
                                     "stem_3" = "Proliferating stem cell (S)",
                                     "stem_4" = "Absorptive progenitor",                                                                                                                                              
                                     "stem_5" = "Terminally proliferated stem cell",                                                                                                                                  
                                     "stem_6" = "Enteroendocrine progenitor",
                                     "2"      = "Early enterocyte",                                                                                                                                                   
                                     "3"      = "Secretory cell",                                                                                                                                                     
                                     "4"      = "Enteroendocrine",
                                     "5_0"    = "Ly6a+ enterocyte progenitor",                                                                                                                                        
                                     "5_1"    = "Ly6a+ enterocyte (late)",                                                                                                                                      
                                     "5_2"    = "Ly6a+ enterocyte (early)",
                                     "6_0"    = "Early secretory cell",                                                                                                                                               
                                     "6_1"    = "Mito high stem cell",                                                                                                                                                
                                     "7"      = "Tuft cell",
                                     "8"      = "Enterocyte")

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

Idents(ISC_annotated_seurat) <- factor(Idents(ISC_annotated_seurat),
                                   levels = cluster_order)
DimPlot(ISC_annotated_seurat) + theme(aspect.ratio = 1)

# 8. ccAFv2 cell cycle scoring -----
library(ccAFv2)
ISC_annotated_seurat <- PredictCellCycle(ISC_annotated_seurat, do_sctransform=TRUE, assay='SCT', species='mouse', gene_id='symbol')

# 9. CytoTRACE2 ----
ISC_annotated_seurat <- cytotrace2(ISC_annotated_seurat, species = "mouse", is_seurat = TRUE, slot_type = "counts", full_model = TRUE, batch_size = 9000, smooth_batch_size = 1000, parallelize_models = TRUE, parallelize_smoothing = TRUE, ncores = NULL, max_pcs = 200, seed = 14)

# 10. Prepare scVelo ----
library(Matrix)

## * ISC_annotated_seurat_ENR ----
ISC_annotated_seurat_ENR <- subset(ISC_annotated_seurat, Sample_Name %in% 'ENR')

# save metadata table:
ISC_annotated_seurat_ENR$barcode <- colnames(ISC_annotated_seurat_ENR)
ISC_annotated_seurat_ENR$UMAP_1 <- ISC_annotated_seurat_ENR@reductions$umap@cell.embeddings[,1]
ISC_annotated_seurat_ENR$UMAP_2 <- ISC_annotated_seurat_ENR@reductions$umap@cell.embeddings[,2]
write.csv(ISC_annotated_seurat_ENR@meta.data, file='ISC_annotated_seurat_ENR_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
counts_matrix <- GetAssayData(ISC_annotated_seurat_ENR, assay='RNA', layer='counts')
writeMM(counts_matrix, file='ISC_annotated_seurat_ENR_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(ISC_annotated_seurat_ENR@reductions$pca@cell.embeddings, file='ISC_annotated_seurat_ENR_pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='ISC_annotated_seurat_ENR_gene_names.csv',
  quote=F,row.names=F,col.names=F
)
#...............................................................................
## * ISC_annotated_seurat_SAG ----
ISC_annotated_seurat_SAG <- subset(ISC_annotated_seurat, Sample_Name %in% 'SAG')

# save metadata table:
ISC_annotated_seurat_SAG$barcode <- colnames(ISC_annotated_seurat_SAG)
ISC_annotated_seurat_SAG$UMAP_1 <- ISC_annotated_seurat_SAG@reductions$umap@cell.embeddings[,1]
ISC_annotated_seurat_SAG$UMAP_2 <- ISC_annotated_seurat_SAG@reductions$umap@cell.embeddings[,2]
write.csv(ISC_annotated_seurat_SAG@meta.data, file='ISC_annotated_seurat_SAG_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
counts_matrix <- GetAssayData(ISC_annotated_seurat_SAG, assay='RNA', layer='counts')
writeMM(counts_matrix, file='ISC_annotated_seurat_SAG_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(ISC_annotated_seurat_SAG@reductions$pca@cell.embeddings, file='ISC_annotated_seurat_SAG_pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='ISC_annotated_seurat_SAG_gene_names.csv',
  quote=F,row.names=F,col.names=F
)
#...............................................................................
## * ISC_annotated_seurat_PGE2 ----
ISC_annotated_seurat_PGE22 <- subset(ISC_annotated_seurat, Sample_Name %in% 'PGE2')

# save metadata table:
ISC_annotated_seurat_PGE2$barcode <- colnames(ISC_annotated_seurat_PGE2)
ISC_annotated_seurat_PGE2$UMAP_1 <- ISC_annotated_seurat_PGE2@reductions$umap@cell.embeddings[,1]
ISC_annotated_seurat_PGE2$UMAP_2 <- ISC_annotated_seurat_PGE2@reductions$umap@cell.embeddings[,2]
write.csv(ISC_annotated_seurat_PGE2@meta.data, file='ISC_annotated_seurat_PGE2_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
counts_matrix <- GetAssayData(ISC_annotated_seurat_PGE2, assay='RNA', layer='counts')
writeMM(counts_matrix, file='ISC_annotated_seurat_PGE2_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(ISC_annotated_seurat_PGE2@reductions$pca@cell.embeddings, file='ISC_annotated_seurat_PGE2_pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='ISC_annotated_seurat_PGE2_gene_names.csv',
  quote=F,row.names=F,col.names=F
)