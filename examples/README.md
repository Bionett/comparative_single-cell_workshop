# Comparative scRNAseq analysis
This hands-on session guides the participants along the data analysis workflow to identify cellular compositional changes and cell-type-specific differentially expressed genes starting from transcript count matrices of a public [dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157783) 
with multiple samples from two clinical conditions.  

## Dataset
This [study](https://academic.oup.com/brain/article/145/3/964/6469020) sequenced over 40.000 nuclei from 11 human midbrain samples, 5 Parkinson's disease (PD) and 6 control. The authors identified 12 major cell population making up the human midbrain. 

![Study - design](https://github.com/Bionett/comparative_single-cell_workshop/blob/main/docsrc/GSE157783_study_design.png?raw=true)

## Aim
Identify the cellular & gene expression changes in the PD midbrain.

## Tissue composition
#### New folder

```bash
mkdir tissue_diffcomp
cd tissue_diffcomp
```

#### Pipeline input parameters
```bash
python ../../pipelines/difftissue/pipeline_difftissue_comp.py config
```

#### Run pipeline
```bash
python ../../pipelines/difftissue/pipeline_difftissue_comp.py make full -v5 --no-cluster
```

#### Output availability
difftissue output can be found [here](https://drive.google.com/file/d/1yC6dBqd_Tmxh-MFhyJxUzNvSaN_iy6uS/view?usp=sharing) (.tar) or [here](https://drive.google.com/file/d/1dKIedhanGnHuSE6ctziLW7Q4lIgGCwfu/view?usp=sharing) (.zip).

#### Data visualization
Let's use Rstudio to explore the results. 

## Differential gene expression
#### New folder

```bash
mkdir celltype_diffexp
cd celltype_diffexp
```

#### Pipeline input parameters
```bash
python ../../pipelines/diffgene/pipeline_genevar.py config
```

#### Run pipeline
```bash
python ../../pipelines/diffgene/pipeline_genevar.py make full -v5 --no-cluster
```

#### Output availability
diffgene output can be found [here](https://drive.google.com/file/d/1HkFpd8eHNy4c7IfQwipGVYb0pjzLvONN/view?usp=sharing) (.tar) or [here](https://drive.google.com/file/d/1di2g6P0Dh6idnId2hfiJ_7yOAl65sueJ/view?usp=sharing) (.zip).

#### Model expression & data visualization
Let's use Rstudio to explore the results. 


