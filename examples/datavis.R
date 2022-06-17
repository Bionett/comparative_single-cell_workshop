# Data vis & gene expression modeling ------------------------------------------

# Load depencies
stopifnot(require(data.table),
          require(dplyr),
          require(ggplot2),
          require(DESeq2))

# Set working directory
setwd("~/comparative_single-cell_workshop")

# difftissue pipeline visualization --------------------------------------------
# Cluster based differential composition ---------------------------------------

# inputs
sccoda_dir <- "examples/tissue_diffcomp/diffcomp.dir/sccoda.dir"
output_file <- "examples/tissue_diffcomp/diffcomp.dir/sccoda.dir/summary.png"

# run
sccoda_res <- list.files(sccoda_dir, pattern = ".tsv.gz", full.names = TRUE)
sccoda_res <- sccoda_res[grep("sccoda_out_ref_", sccoda_res)]

names(sccoda_res) <- gsub("sccoda_out_ref_|.tsv.gz", "", basename(sccoda_res))
sccoda_res

sccoda_tab <- do.call(rbind,
                      lapply(names(sccoda_res), function(refcell) {
                        sccoda <- fread(sccoda_res[[refcell]])
                        sccoda[["refcell"]] <- refcell
                        sccoda
                      })
)

sccoda_tab$Covariate <- gsub("\\s|C|\\(|condition|,|Treatment|ontrol|\\)|\\[|T|\\.|\\]|\\'",
                             "", sccoda_tab$Covariate)
head(sccoda_tab)

sccoda_tab %>%
  mutate('level_celltype' = paste0(Covariate, '_', `Cell Type`)) %>%
  group_by(level_celltype) %>%
  summarise('mean_log2FC' = mean(`log2-fold change`),
            'sd' = sd(`log2-fold change`),
            'credible' = (length(which(credible))/n())*100,
            'cell_type' = unique(`Cell Type`),
            'covar' = unique(Covariate)) -> sccoda_sum

unique(sccoda_sum$covar)

gg <- ggplot(sccoda_sum, aes(x = cell_type, y = mean_log2FC,
                             color = credible)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_log2FC-`sd`, ymax = mean_log2FC+`sd`),
                width = 0.2) +
  theme_classic() +
  scale_color_gradient(low = "grey50", high = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.1) +
  ylim(c(-2, 4)) +
  coord_flip() +
  facet_wrap(~covar, nrow = 1)

plot(gg)

pdf(gsub(".png", paste0("_scCODA_condition", ".pdf"), output_file),
    width = 4, height = 4)
  plot(gg)
dev.off()

gg <- ggplot(sccoda_sum, aes(x = cell_type, y = mean_log2FC, fill = credible)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean_log2FC-`sd`, ymax = mean_log2FC+`sd`),
                width = 0.2) +
  theme_classic() +
  scale_fill_gradient(low = "grey50", high = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.1) +
  ylim(c(-2, 4)) +
  coord_flip() +
  facet_wrap(~covar, nrow = 1)

plot(gg)

pdf(gsub(".png", paste0("_scCODA_condition_barplot.pdf"), output_file),
    width = 4, height = 4)
  plot(gg)
dev.off()

# Cluster independent differential composition ---------------------------------
# inputs
cna_dir <- "examples/tissue_diffcomp/diffcomp.dir/cna.dir"
umap_file <- "data/metadata/umap.tsv.gz"
meta_file <- "data/metadata/cell_metadata.tsv.gz"
output_file <- "examples/tissue_diffcomp/diffcomp.dir/cna.dir/summary.png"

cna_res <- list.files(cna_dir, pattern = "cell_corr", full.names = TRUE)

names(cna_res) <- gsub("cna_out_|_cell_corr.tsv.gz", "", basename(cna_res))
cna_res

cna_tab <- do.call(rbind,
                      lapply(names(cna_res), function(condition) {
                        cna <- fread(cna_res[[condition]], header = T)
                        colnames(cna)[2] <- "cell_corr"
                        cna[["condition"]] <- condition
                        cna %>%
                          select(barcode_id, condition, cell_corr)
                      })
)

umap <- fread(umap_file)
head(umap)

cna_tab <- merge(cna_tab, umap, by.x = "barcode_id", by.y = "barcode")
head(cna_tab)

gg <- ggplot(cna_tab, aes(x = UMAP_1, y = UMAP_2, color = cell_corr)) +
  geom_point(size = 0.05) +
  theme_classic() +
  scale_color_gradient2(low = "darkblue", high = "darkred", midpoint = 0) +
  facet_wrap(~condition)
  
plot(gg)

pdf(gsub(".png", paste0("_cna_cellcorr_umap.pdf"), output_file),
    width = 7, height = 4)
  plot(gg)
dev.off()

# genes [?]
head(fread("examples/tissue_diffcomp/diffcomp.dir/cna.dir/cna_out_IPD_gene_corr.tsv.gz"), 20)

# diffgene pipeline visualization ----------------------------------------------
# DESeq2 framework -------------------------------------------------------------

pbulk_file <- "examples/celltype_diffexp/genevar.dir/pbulk.dir/pbulk.tsv.gz"
meta_file <- "examples/celltype_diffexp/genevar.dir/pvca.dir/pvca_sample_metada.tsv"
cellvar <- "cell_ontology"
mincount <- 200
output_file <- "examples/celltype_diffexp/genevar.dir/deg.tsv"


pbulk <- fread(pbulk_file)
pbulk[1:3, 1:3]

pbulk_mat <- data.matrix(select(pbulk, -1))
rownames(pbulk_mat) <- pbulk[["gene_name"]]
pbulk_mat[1:3, 1:3]

meta <- read.delim(meta_file)
rownames(meta) <- meta[["sample_cell"]]
head(meta)

celltypes <- unique(meta[[cellvar]])
celltypes

lapply(celltypes, function(celltype) {
  
  print(celltype)
  
  meta %>%
    filter(!!rlang::sym(cellvar) == celltype) -> sub_meta
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = pbulk_mat[, rownames(sub_meta)],
                                colData = sub_meta,
                                design = ~condition)
  
  dds <- scran::computeSumFactors(dds)
  
  keep <- rowSums(counts(dds)) > mincount
  dds <- dds[keep, ]
  
  dds <- DESeq(dds)
  resultsNames(dds) # lists the coefficients
  res <- results(dds)
  
  write.table(res, gsub(".tsv", paste0("_", celltype, "_minmodel.tsv"), 
                        output_file),
              sep = "\t", row.names = F, quote = F)
  
})