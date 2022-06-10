# ------------------------------------------------------------------------------
# Prepare tables for scCODA statistical analysis
# Inputs: a cell metadata table, a sample variable name, and a cell cluster 
# variable
# ------------------------------------------------------------------------------

stopifnot(
  require(optparse),
  require(futile.logger),
  require(BiocParallel),
  require(R.utils),
  require(dplyr),
  require(data.table)
)

option_list <- list(
  make_option(c("--cellmeta"), default="cell.table.tsv.gz",
              help="File with the barcode metadata. Output of the fetch cells
              pipeline."),
  make_option(c("--cellvar"), default="barcode_id",
              help="Sample index/column/var-name in samplemeta file."),
  make_option(c("--curated_cellmeta"), default=NULL,
              help="File with curated cell metadata."),
  make_option(c("--cur_cellvar"), default="barcode_id",
              help="Sample index/column/var-name in samplemeta file."),
  make_option(c("--samplevariable"), default="library_id",
              help="Sample variable in cell metadata."),
  make_option(c("--samplemetavars"), default=NULL,
              help="Sample variable in cell metadata. If more than one,
              seprate them with ,"),
  make_option(c("--clustervariable"), default=NULL,
              help="Cell variable in cell metadata to gather cells."),
  make_option(c("--numcores"), default=4,
              help="Number of cores used to ..."),
  make_option(c("--log_filename"), 
              default="./diffcomp.dir/sccoda.dir/prepare_sccoda_tables.log"),
  make_option(c("--outfile"), 
              default="./diffcomp.dir/sccoda.dir/sccoda_input.tsv")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Logger ------

# Prepare the logger and the log file location 
# the logger is by default writing to ROOT logger and the INFO threshold level -
###################
flog.threshold(INFO)

# now set to append mode -------------------------------------------------------
###################
flog.appender(appender.file(opt$log_filename))
flog.info("Running with parameters:", opt, capture = TRUE)

multicoreParam <- MulticoreParam(workers = opt$numcores)

# Read in data -----------------------------------------------------------------
# Load cell metadata -----------------------------------------------------------

if(!file.exists(opt$cellmeta)) {
  stop("Non existing cell metadata file provided.")
}

meta <- fread(opt$cellmeta)
meta[meta == ""] <- NA
print(head(meta,2))
print(dim(meta))

# Sample centered sequencing & mapping metrics ---------------------------------
# If curated sample metadata provided ------------------------------------------

if(!file.exists(opt$curated_cellmeta)) {
  flog.info("Non existing curated cell metadata file provided.")
} else {
  cmeta <- fread(opt$curated_cellmeta)
  cmeta[cmeta == ""] <- NA
  print(head(cmeta,2))
  print(dim(cmeta))
  
  meta <- merge(meta, select(cmeta, c(!!rlang::sym(opt$clustervariable), barcode_id)), 
                by = "barcode_id")
  
  flog.info("Merged curated and cell metadata:",  head(meta,3), capture = TRUE)
}

# Cell count -------------------------------------------------------------------

meta %>%
  group_by(!!rlang::sym(opt$samplevariable), !!rlang::sym(opt$clustervariable)) %>%
  summarise('ncell' = n()) %>%
  tidyr::spread(!!rlang::sym(opt$clustervariable), ncell) %>%
  arrange(!!rlang::sym(opt$samplevariable)) %>%
  data.table -> ncell

ncell[is.na(ncell)] <- 0
head(ncell)

write.table(ncell, gsub("\\.tsv", "_cell_count_mat.tsv", opt$outfile), 
            sep = "\t", quote = FALSE, row.names = FALSE)
gzip(gsub("\\.tsv", "_cell_count_mat.tsv", opt$outfile))

# Cell metadata summary table --------------------------------------------------

metavars <- unlist(strsplit(opt$samplemetavars, ","))
metavars

sample_meta <- Reduce(function(x, y) merge(x, y, by = opt$samplevariable), 
                      lapply(metavars, function(v) {
                        meta %>%
                          group_by(!!rlang::sym(opt$samplevariable)) %>%
                          summarise(v = unique(!!rlang::sym(v)))
                        })
) 

colnames(sample_meta)[2:ncol(sample_meta)] <- metavars
rownames(sample_meta) <- sample_meta[[opt$samplevariable]]
head(sample_meta)

write.table(sample_meta[ncell[[opt$samplevariable]], ], 
            gsub("\\.tsv", "_sample_meta.tsv", opt$outfile), 
            sep = "\t", quote = FALSE, row.names = FALSE)

gzip(gsub("\\.tsv", "_sample_meta.tsv", opt$outfile))
