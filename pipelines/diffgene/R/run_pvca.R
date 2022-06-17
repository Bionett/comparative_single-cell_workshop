# ------------------------------------------------------------------------------
# Run PVCA given a pseudobulk gene expression matrix *.tsv.gz
# ------------------------------------------------------------------------------

stopifnot(
  require(optparse),
  require(futile.logger),
  require(BiocParallel),
  require(R.utils),
  require(dplyr),
  require(data.table),
  require(ggfortify),
  require(DESeq2),
  require(pvca),
  require(ggradar),
  require(ComplexHeatmap),
  require(circlize)
)

option_list <- list(
  make_option(c("--pbulk"), default="./diffexp.dir/pbulk.dir/pbulk.tsv.gz",
              help="File with the pbulk [gene, celltype-sample] matrix. *tsv.gz. Gene index in column 1"),
  make_option(c("--cellmeta"), default="cell.table.tsv.gz",
              help="File with the barcode metadata. Output of the fetch cells
              pipeline."),
  make_option(c("--cell_var"), default=NULL,
              help="Cell-cluster variable."),
  make_option(c("--sample_var"), default=NULL,
              help="Celltype-sample variables separated by comma to explain gene expression variation."),
  make_option(c("--cellvars"), default=NULL,
              help="Celltype-sample variables separated by comma to explain gene expression variation."),
  make_option(c("--sample_meta"), default=NULL,
              help="File with curated sample-celltype metadata."),
  make_option(c("--mingenecount"), default=200,
              help="Minimum number of UMI counts to consider gene."),
  make_option(c("--mincellsample"), default=10,
              help="Minimum number of cells per sample-celltype."),
  make_option(c("--numcores"), default=1,
              help="Number of cores used to ..."),
  make_option(c("--log_filename"), default="./diffexp.dir/pvca.dir/pvca.log"),
  make_option(c("--outfile"),
              default="./diffexp.dir/pvca.dir/metadata_var_weights.tsv")
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

### ---- Aux functions 
flog.info("Load aux functions ... ")
get_vst_mat <- function(pbulk, samplemeta, mingenecount, samplevar, drawpca) {
  
  # DESeq2
  dds <- DESeqDataSetFromMatrix(countData = pbulk,
                                colData = samplemeta,
                                design = ~1)
  
  dds <- scran::computeSumFactors(dds)
  
  # Subset low-expressed genes
  print(quantile(rowSums(counts(dds))))
  keep <- rowSums(counts(dds)) > mingenecount
  dds <- dds[keep, ]
  print(quantile(rowSums(counts(dds))))
  
  # Variance stabalizing normalization
  print(dim(dds))
  if(nrow(dds) < 1000) {
    vsd <- varianceStabilizingTransformation(dds) 
  } else {
    vsd <- vst(dds, blind=TRUE)
  }
  
  print(dim(assay(vsd)))
  vsd_mat <- assay(vsd)
  colnames(vsd_mat) <- colnames(pbulk)
  vs <- colnames(samplemeta)
  
  if(drawpca) {
    tr_vsd <- data.matrix(t(vsd_mat))
    print(tr_vsd[1:2, 1:2])
    
    # PCA plots colouring by meta variable
    pdf(gsub(".tsv", paste0("_", "pbulk_pca_vst_norm_", samplevar, ".pdf"), 
             opt$outfile), width = 13, height = 13)
    for(v in vs) {
      plot(
        autoplot(prcomp(tr_vsd), 
                 data=samplemeta,
                 colour = v,
                 label = FALSE, 
                 label.size = 3) +
          theme_classic() +
          theme(legend.position = "bottom")
      )
    }
    dev.off()
  }
  
  vsd_mat
  
}
get_pheno_weights <- function(norm_mat, samplemeta, samplevar, celltype=NULL) {
  
  pheno <- data.frame(samplemeta)
  
  print(head(pheno, 2))
  
  my_set = new("ExpressionSet", exprs=norm_mat, 
               phenoData=AnnotatedDataFrame(pheno))
  
  pvca_res <- tryCatch({
    
    pvcaBatchAssess(my_set, colnames(pheno), threshold = 0.75)
    
  }, error = function(e) {
    
    flog.info("PVCA calculation failed.")
    flog.info("error: ", e, capture = TRUE)

    flog.info("For var: ", samplevar, capture = TRUE)
    return(NULL)
    
  })
  
  if(!is.null(pvca_res)) {
    
    pvca_df <- data.frame('Variable' = pvca_res$label,
                          'pct_variance' = as.vector(pvca_res$dat)) %>%
      arrange(pct_variance)
    
    print(head(pvca_df))
    
    if(!is.null(celltype)) {
      samplevar <- paste0(samplevar[1], "_", samplevar[2], "_", celltype)
    } else {
      samplevar <- samplevar[2]
    }
    
    write.table(pvca_df, gsub(".tsv", 
                              paste0("_", "pbulk_pvca_deseq2_VSTnorm_varRank_",
                                     samplevar, 
                                     ".tsv"), opt$outfile),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    gg_pvca <- ggplot(pvca_df, aes(x = Variable, y = pct_variance*100)) +
      geom_bar(stat = "identity") +
      scale_x_discrete(limits = pvca_df$Variable) +
      theme_classic() +
      coord_flip()
    
    pdf(gsub(".tsv", paste0("_", "pseudobulk_pvca_barplot_", 
                            samplevar, 
                            ".pdf"),
             opt$outfile), width = 3, height = 3)
      plot(gg_pvca)
    dev.off()
    
    pvca_df[["sample_var"]] <- samplevar
    return(pvca_df)
    
  } else {
    NULL
  }
}

# Read in data -----------------------------------------------------------------
# Load cell metadata -----------------------------------------------------------

if(!file.exists(opt$pbulk)) {
  stop("Non existing pbulk metadata file provided.")
}

mtx <- fread(opt$pbulk)
mat <- data.matrix(mtx[, -1])
rownames(mat) <- mtx[[1]]
flog.info("Pbulk matrix:", mat[1:5, 1:2])

if(!file.exists(opt$cellmeta)) {
  stop("Non existing cell metadata file provided.")
}

meta <- fread(opt$cellmeta)

if(opt$sample_var %in% colnames(meta)) {
  if(opt$cell_var %in% colnames(meta)) {
    meta %>%
      mutate("sample_cell" = paste0(!!rlang::sym(opt$sample_var), "_", 
                                    !!rlang::sym(opt$cell_var))) -> meta
    flog.info('Sample #cells distribution:', 
              quantile(table(meta[["sample_cell"]])), 
              capture = TRUE)
    
    isamples <- names(which(table(meta[["sample_cell"]]) > opt$mincellsample))
    meta %>%
      filter(sample_cell %in% isamples) -> meta
    
    flog.info('Sample #cells distribution (filtered):', 
              quantile(table(meta[["sample_cell"]])), 
              capture = TRUE)
    
  } else {
    stop("cell_var provided not in cell metadata.")
  } 
} else {
  stop("sample_var provided not in cell metadata.")
}


if(!is.null(opt$cellvars)) {
  cellvars <- unlist(strsplit(opt$cellvars, ","))
  flog.info("Celltype-sample variables to explore: ", cellvars, capture = TRUE)
  cellvars <- cellvars[cellvars %in% colnames(meta)]
  if(length(cellvars) == 0) {
    stop("Provided celltype-sample variables no present in cell metadata.")
  } else {
    flog.info("Celltype-sample variables to explore: ", cellvars, capture = TRUE)
  }
} else {
  stop("No sample variable provided to explain gene expression variance.")
}

cellvars <- append(cellvars, opt$sample_var)
cellvars <- append(cellvars, opt$cell_var)

sample_ms <- lapply(cellvars, function(cv) {
  if(length(unique(meta[[cv]])) > 1) {
    meta %>%
      group_by(sample_cell) %>%
      summarise('var' = unique(!!rlang::sym(cv))) -> sample_meta
    colnames(sample_meta)[2] <- cv
    sample_meta
  } else {
    flog.info("Sample variable only has one level:", cv, capture = TRUE)
    NULL
  }
})

sample_ms <- Filter(Negate(is.null), sample_ms)

sample_meta <- Reduce(function(x, y) merge(x, y, by = "sample_cell"), sample_ms)

flog.info("Sample metadata:", head(sample_meta, 3), capture = TRUE)

write.table(sample_meta, gsub(".tsv", "_sample_metada.tsv", opt$outfile), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# PVCA
cellvars <- colnames(sample_meta[-1])
flog.info("Variables to rank:", cellvars, capture = TRUE)

pvcas <- lapply(cellvars, function(cv) {
  print(cv)
  flog.info("PVCA... ", cv, capture = TRUE)
  
  if(cv == opt$cell_var) {
    
    return(NULL) 
    
  } else {
    
    sample_meta %>%
      select(c(!!rlang::sym(opt$cell_var), 
               !!rlang::sym(cv))) %>%
      data.frame -> sub_sample_meta
    
    rownames(sub_sample_meta) <- sample_meta[["sample_cell"]]
    
    if(cv == opt$sample_var) {
      
      sub_sample_meta[["random"]] <- sample(LETTERS[1:2],
                                            nrow(sub_sample_meta), 
                                            replace=TRUE, 
                                            prob=c(0.5, 0.5))
      print(head(sub_sample_meta, 3))
      print(table(sub_sample_meta[[opt$sample_var]]))
      cellvar <- "random"
      
    } else {
      
      cellvar <- opt$cell_var
      
    }
    
    cellvar <- gsub(" ", ".", cellvar)
    
    sub_sample_meta %>%
      group_by(!!rlang::sym(cellvar),
               !!rlang::sym(cv)) %>%
      summarise('nsamples' = n()) %>%
      arrange(nsamples) -> nsample
    
    flog.info("#sample: ", head(nsample), capture = TRUE)
       
    if(1 %in% nsample$nsamples) {
      
      flog.info(cv, " variable create single-sample levels combined with the cell-type variable.",
                capture = TRUE)
      flog.info("#sample distribution:", quantile(nsample$nsamples), capture = TRUE)
      
      flog.info("Random 2 level variable as reference.")
      
      sub_sample_meta[["random"]] <- sample(LETTERS[1:2],
                                            nrow(sub_sample_meta), 
                                            replace=TRUE, 
                                            prob=c(0.5, 0.5))
      print(head(sub_sample_meta, 3))
      print(table(sub_sample_meta[[opt$sample_var]]))
      cellvar <- "random"
      
      sub_sample_meta %>%
        group_by(!!rlang::sym(cellvar),
                 !!rlang::sym(cv)) %>%
        summarise('nsamples' = n()) %>%
        arrange(nsamples) -> nsample
      
      flog.info("# sample after new var created", 
                head(nsample), capture = TRUE)

      if(1 %in% nsample$nsamples) { 
        
        return(NULL)
        
      } else {
        
        # DESeq2 Normalize pbulk matrix
        flog.info("DESeq2 ... VST normalization.")
        vst_mat <- get_vst_mat(pbulk = mat[, rownames(sub_sample_meta)], 
                               samplemeta = sub_sample_meta, 
                               mingenecount = opt$mingenecount, 
                               samplevar = cv, 
                               drawpca = TRUE)
      
        # PVCA
        flog.info("PVCA ...")
        get_pheno_weights(norm_mat = vst_mat[, rownames(sub_sample_meta)],
                          samplemeta = sub_sample_meta,
                          samplevar = c(cellvar, cv))
      }
      
    } else {
      
      # DESeq2 Normalize pbulk matrix
      flog.info("DESeq2 ... VST normalization.")
      vst_mat <- get_vst_mat(pbulk = mat[, rownames(sub_sample_meta)], 
                             samplemeta = sub_sample_meta, 
                             mingenecount = opt$mingenecount, 
                             samplevar = cv, 
                             drawpca = TRUE)
      
      # PVCA
      flog.info("PVCA ...")
      get_pheno_weights(norm_mat = vst_mat[, rownames(sub_sample_meta)],
                        samplemeta = sub_sample_meta,
                        samplevar = c(cellvar, cv))
    }
  }
})

pvcas <- Filter(Negate(is.null), pvcas)
pvca_res <- do.call(rbind, pvcas)
flog.info("Savce PVCA results: ", head(pvca_res, 3), capture = TRUE)

write.table(pvca_res, gsub(".tsv", 
                           "_pvca_var_loads.tsv", opt$outfile),
            sep = "\t", quote = FALSE, row.names = FALSE)

flog.info("Saved global pvca_loads.tsv")

################################################################################

flog.info("Run PVCA per cell-type ... ")

celltypes <- unique(sample_meta[[opt$cell_var]])

pvcas <- lapply(celltypes, function(ct) {
  
  print(ct)
  
  sample_meta %>%
    filter(!!rlang::sym(opt$cell_var) == ct) -> sub_sample_meta
  
  if(nrow(sub_sample_meta) > 3) {
    
    pvcas <- lapply(cellvars, function(cv) {
      
      print(cv)
      flog.info("PVCA... ", cv, capture = TRUE)
      
      if(cv %in% c(opt$cell_var, opt$sample_var)) {
        
        return(NULL) 
        
      } else if(length(unique(sub_sample_meta[[cv]])) <= 1){
        
        return(NULL)
        
      } else {
        
        sub_sample_meta %>%
          select(c(!!rlang::sym(cv))) %>%
          data.frame -> ss_sample_meta
        
        rownames(ss_sample_meta) <- sub_sample_meta[["sample_cell"]]
        
        ss_sample_meta[["random"]] <- sample(LETTERS[1:2],
                                              nrow(ss_sample_meta), 
                                              replace=TRUE, 
                                              prob=c(0.5, 0.5))
        print(head(ss_sample_meta, 3))
        print(table(ss_sample_meta[[opt$sample_var]]))
        cellvar <- "random"
        
        ss_sample_meta %>%
          group_by(!!rlang::sym(cellvar),
                   !!rlang::sym(cv)) %>%
          summarise('nsamples' = n()) %>%
          arrange(nsamples) -> nsample
          
        print(head(nsample))
        
        if(1 %in% nsample$nsamples) {
          
          flog.info(cv, " variable create single-sample levels combined with the cell-type variable.",
                    capture = TRUE)
          flog.info("#sample distribution:", quantile(nsample$nsamples), capture = TRUE)
          
          return(NULL)
          
        } else {
          
          # DESeq2 Normalize pbulk matrix
          flog.info("DESeq2 ... VST normalization.")
          vst_mat <- get_vst_mat(pbulk = mat[, rownames(ss_sample_meta)], 
                                 samplemeta = ss_sample_meta, 
                                 mingenecount = opt$mingenecount, 
                                 samplevar = paste0(cv, "_", gsub("/", "", ct)), 
                                 drawpca = TRUE)
          
          # PVCA
          flog.info("PVCA ...")
          print(head(ss_sample_meta))
          get_pheno_weights(norm_mat = vst_mat[, rownames(ss_sample_meta)],
                            samplemeta = ss_sample_meta,
                            samplevar = c(cellvar, cv),
                            celltype = gsub("/", "", ct))
        }
      }
    })
    
    pvcas <- Filter(Negate(is.null), pvcas)
    
    if(length(pvcas) > 0) {
      
      pvca_res <- do.call(rbind, pvcas)
      pvca_res[["cell_type"]] <- ct
      
      flog.info("Savce PVCA results: ", head(pvca_res, 3), capture = TRUE)
      
      write.table(pvca_res, gsub(".tsv", 
                                 paste0("_pvca_var_loads_", ct, ".tsv"), 
                                 opt$outfile),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      pvca_res
      
    } else {
      
      NULL
      
    }
    
  
  } else {
    
    NULL
    
  }
})

pvcas <- Filter(Negate(is.null), pvcas)

if(length(pvcas) > 0) {
  pvca_res <- do.call(rbind, pvcas)
  flog.info("Save PVCA cell-type results: ", head(pvca_res, 3), capture = TRUE)
  
  write.table(pvca_res, gsub(".tsv", 
                             paste0("_pvca_var_loads_per_celltype.tsv"), 
                             opt$outfile),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

flog.info("PVCA done.")
