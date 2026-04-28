# 0. setting ----
BASE_PATH  <- "/Users/choi_lab_sk/Dropbox/2026_Hh_ISC_regeneration/bulkrnaseq"
OUTPUT_DIR <- file.path(BASE_PATH, "data", "preprocessing")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

set.seed(1)

# 1. Load libraries ----
library(tximport)
library(DESeq2)
library(readr)
library(dplyr)
library(tidyverse)
library(apeglm)
library(DEGreport)

# 2. Load Salmon quantification files ----
salmon_dir <- file.path(BASE_PATH,
                        "GSE302173")
samples <- list.files(salmon_dir, full.names = TRUE, pattern = "_quant$")
files   <- setNames(file.path(samples, "quant.sf"),
                    basename(samples) %>%
                      str_replace("_quant$", ""))

# 3. tx2gene table (GENCODE vM34) ----
tx2gene_path <- file.path(BASE_PATH, "meta/tx2gene.gencode.vM34.csv")
tx2gene      <- read.delim(tx2gene_path)
colnames(tx2gene) <- c("tx_id","ensgene","symbol")

txi <- tximport(files, type = "salmon", 
                tx2gene = tx2gene[,c("tx_id","ensgene")], 
                countsFromAbundance = "lengthScaledTPM")

names(txi) <- c("abundance","counts","length","countsFromAbundance")

data <- txi$counts %>% 
  round() %>% data.frame()

sampletype <- factor(c(rep(c("Homeostasis","Lgr5_DTA","Lgr5_DTA_Cyclopamine"),2)))
meta <- data.frame(sampletype, row.names = colnames(txi$counts))
meta

all(colnames(txi$counts) %in% rownames(meta)) #TRUE
all(colnames(txi$counts) == rownames(meta))   #TRUE

meta$sampletype <- factor(meta$sampletype,
                             levels = c("Homeostasis",
                                        "Lgr5_DTA",
                                        "Lgr5_DTA_Cyclopamine"))

# 4. DESeq2 object & pre-filtering ----
dds  <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sampletype)
keep <- rowSums(counts(dds) >= 10) >= 2
dds  <- dds[keep, ]

# 5. Normalization & rlog transformation ----
dds               <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=T) %>% data.frame() %>% rownames_to_column(var="gene")
rld               <- rlog(dds, blind = TRUE)

vM34annot <- tx2gene %>% 
  dplyr::select(ensgene, symbol) %>% dplyr::distinct()

normalized_counts <- merge(normalized_counts, vM34annot, by.x="gene",by.y="ensgene") %>% as_tibble()
head(normalized_counts)

rownames(dds) %>% head
normalized_counts %>% head

rownames(dds) %>% tail
normalized_counts %>% tail

rownames(dds) <- normalized_counts$symbol

write.csv(as.data.frame(normalized_counts),
          file.path(OUTPUT_DIR, "normalized_counts.csv"))

# 6. Wald test ----
dds_wald <- DESeq(dds, test = "Wald")

## Unshrunken results
res_d_h_raw  <- results(dds_wald,
                        contrast = c("sampletype", "Lgr5_DTA", "Homeostasis"))
res_dc_h_raw <- results(dds_wald,
                        contrast = c("sampletype", "Lgr5_DTA_Cyclopamine", "Homeostasis"))

dds_relevel            <- dds
dds_relevel$sampletype <- relevel(dds_relevel$sampletype, ref = "Lgr5_DTA")
dds_relevel            <- DESeq(dds_relevel, test = "Wald")
res_dc_d_raw <- results(dds_relevel,
                        contrast = c("sampletype", "Lgr5_DTA_Cyclopamine", "Lgr5_DTA"))

res_h_d_raw <- results(dds_relevel,
                        contrast = c("sampletype", "Homeostasis", "Lgr5_DTA"))

dds_relevel2            <- dds
dds_relevel2$sampletype <- relevel(dds_relevel2$sampletype, ref = "Lgr5_DTA_Cyclopamine")
dds_relevel2            <- DESeq(dds_relevel2, test = "Wald")
res_d_dc_raw <- results(dds_relevel2,
                        contrast = c("sampletype", "Lgr5_DTA", "Lgr5_DTA_Cyclopamine"))

## Shrunken LFC (apeglm) for volcano plots and DEG tables
res_d_h_shrunk  <- lfcShrink(dds_wald,
                             coef = "sampletype_Lgr5_DTA_vs_Homeostasis",
                             type = "apeglm")
res_dc_h_shrunk <- lfcShrink(dds_wald,
                             coef = "sampletype_Lgr5_DTA_Cyclopamine_vs_Homeostasis",
                             type = "apeglm")
res_dc_d_shrunk <- lfcShrink(dds_relevel,
                             coef = "sampletype_Lgr5_DTA_Cyclopamine_vs_Lgr5_DTA",
                             type = "apeglm")
res_h_d_shrunk <- lfcShrink(dds_relevel,
                             coef = "sampletype_Homeostasis_vs_Lgr5_DTA",
                             type = "apeglm")
res_d_dc_shrunk <- lfcShrink(dds_relevel2,
                             coef = "sampletype_Lgr5_DTA_vs_Lgr5_DTA_Cyclopamine",
                             type = "apeglm")
## DEG tables
to_df <- function(res) { df <- as.data.frame(res); df$symbol <- rownames(df); df }
res_d_h_df  <- to_df(res_d_h_shrunk)
res_dc_h_df <- to_df(res_dc_h_shrunk)
res_dc_d_df <- to_df(res_dc_d_shrunk)
res_h_d_df <- to_df(res_h_d_shrunk)
res_d_dc_df <- to_df(res_d_dc_shrunk)

write.csv(subset(res_d_h_df,  padj < 0.05 & !is.na(padj)),
          file.path(OUTPUT_DIR, "DEG_Lgr5DTA_vs_Homeostasis.csv"))
write.csv(subset(res_dc_h_df, padj < 0.05 & !is.na(padj)),
          file.path(OUTPUT_DIR, "DEG_Lgr5DTA_Cyclopamine_vs_Homeostasis.csv"))
write.csv(subset(res_dc_d_df, padj < 0.05 & !is.na(padj)),
          file.path(OUTPUT_DIR, "DEG_Lgr5DTA_Cyclopamine_vs_Lgr5DTA.csv"))
write.csv(subset(res_h_d_df, padj < 0.05 & !is.na(padj)),
          file.path(OUTPUT_DIR, "DEG_Homeostasis_vs_Lgr5DTA.csv"))
write.csv(subset(res_d_dc_df, padj < 0.05 & !is.na(padj)),
          file.path(OUTPUT_DIR, "DEG_Lgr5DTA_vs_Lgr5DTA_Cyclopamine.csv"))

## GSEA-ready gene lists (Wald stat from unshrunken results)
make_gene_list <- function(res_raw) {
  df <- as.data.frame(res_raw)
  df$symbol <- rownames(df)
  df <- df[!is.na(df$stat) & df$symbol != "", ]
  stats <- setNames(df$stat, df$symbol)
  sort(stats[!duplicated(names(stats))], decreasing = TRUE)
}
gene_list_d_h  <- make_gene_list(res_d_h_raw)
gene_list_dc_h <- make_gene_list(res_dc_h_raw)
gene_list_dc_d <- make_gene_list(res_dc_d_raw)
gene_list_h_d <- make_gene_list(res_h_d_raw)
gene_list_d_dc <- make_gene_list(res_d_dc_raw)

# 7. Likelihood Ratio Test (LRT) ----
design(dds) = ~ sampletype
dds_LRT <- DESeq(dds, test="LRT", reduced = ~1)

# 8. Save RDS objects ----
saveRDS(dds_wald,   file.path(OUTPUT_DIR, "dds_wald.rds"))
saveRDS(dds_relevel,file.path(OUTPUT_DIR, "dds_relevel.rds"))
saveRDS(dds_lrt,    file.path(OUTPUT_DIR, "dds_lrt.rds"))
saveRDS(rld,        file.path(OUTPUT_DIR, "rld.rds"))
saveRDS(list(d_h  = gene_list_d_h,
             dc_h = gene_list_dc_h,
             dc_d = gene_list_dc_d,
             h_d = gene_list_h_d,
             d_dc = gene_list_d_dc),
        file.path(OUTPUT_DIR, "gsea_gene_lists.rds"))