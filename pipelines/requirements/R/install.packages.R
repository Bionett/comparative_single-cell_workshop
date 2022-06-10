## Default repo
local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.r-project.org" 
       options(repos=r)
})

if(!rlang::is_installed("devtools")) {
    install.packages("devtools")
}
if(!rlang::is_installed("BiocManager")) {
    install.packages("BiocManager")
}

install.package <- function(pkg, src="CRAN", github_repo="") {
    if(!!rlang::is_installed(pkg)) {

        message("installing missing package: ", pkg)

        if(src %in% c("CRAN", "bioconductor")) {

                BiocManager::install(pkg)

            } else if (source == "github") {

                devtools::install_github(github_repo)
            } else {

		stop("Invalid package source")
            }
    }
}

cran_bio_pkgs <- c("optparse",
                   "futile.logger",
                   "BiocParallel",
                   "R.utils",
                   "data.table",
                   "dplyr",
                   "ggplot2",
                   "ggfortify",
                   "ComplexHeatmap",
                   "ggradar",
                   "circlize",
                   "DESeq2",
                   "pvca")


message("Installing cran/bioconductor packages")
for(p in cran_bio_pkgs)
{
    install.package(p, src="CRAN")
}
message("Done, machine ready to work ;)")
