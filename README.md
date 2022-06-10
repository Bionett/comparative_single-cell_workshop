# Introduction to comparative single-cell transcriptomics data-analysis

This repository hosts the hands-on component of the [introductory workshop on single-cell gene expression data analysis](https://www.medsci.ox.ac.uk/divisional-services/support-services-1/business-partnerships-office/industry-fellows-and-postdocs-network/attend-our-themed-workshops). It intends to guide participants on how to configure their computing machines and to present a couple of pipelined implementations to address two key comparative questions:

* What are the **tissue composition differences** between experimental conditions?
* Which are the **cell-type-specific gene expression signatures** associated with each experimenal condition?

## Getting started
To perform the workshop hands-on section, first, make sure you satisfy the systems requirements and then proceed with the R and Python specific dependencies.

### System requirements
Participants are expected to have access to a computing machine with the following software:

* [R version >=4.0](https://www.r-project.org/)
* [Rstudio](https://www.rstudio.com/products/rstudio/download/)
* [Python version >=3.8](https://www.python.org/about/gettingstarted/)
* [Jupyter notebooks](https://jupyter.org/install)

For installation guidance please click on each software.

Then, clone this repository locally and install its dependencies.

```bash
git clone https://github.com/Bionett/comparative_single-cell_workshop.git
```

#### Note: to run the [example](examples/) comparative analysis, it is necessary to install some dependencies:

##### For python:
```bash
pip install -r pipelines/requirements/python/requirements.txt
```

##### For R:

```bash
Rscript pipelines/requirements/R/install.packages.R
```
