# Differential tissue composition & gene expression pipelines

These two [cgat-core](https://cgat-core.readthedocs.io/) based pipelines, [difftissue](difftissue/) and [diffgene](diffgene/) implement a a handful of state-of-the-art methods to identify differences in tissue cellular composition and gene expression.

## [`difftissue`](difftissue/)
This pipeline implements two methodologies for tissue compositional differences:

* a cluster-independent: **[Co-varying neighborhood analysis](https://www.nature.com/articles/s41587-021-01066-4)** and 
* a cluster-dependent: **[scCODA](https://www.nature.com/articles/s41467-021-27150-6)**

## [`diffgene`](diffgene/)
This pipeline performs pseudo-bulk-based differential gene expression based on the [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) modeling framework.
