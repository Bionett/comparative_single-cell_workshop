# ------------------------------------------------------------------------------
# Run scCODA 
# ------------------------------------------------------------------------------

import os
import sys
import scanpy as sc
import pandas as pd
import anndata as ad
import logging
import argparse
import matplotlib.pyplot as plt
import warnings
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz
from sccoda.util import comp_ana as mod
import pickle as pkl

# ############################################################################ #
# ###################### Set up the logging ################################## #
# ############################################################################ #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("build_cellxgene_public.py")

# ############################################################################ #
# ######################## Parse the arguments ############################### #
# ############################################################################ #

parser = argparse.ArgumentParser()
parser.add_argument("--cell_count_file",
                    default="./diffcomp.dir/sccoda.dir/sccoda_input_cell_count_mat.tsv.gz",
                    type=str, help='Table with the cell numbers by sample.')
parser.add_argument("--sample_meta",
                    default="./diffcomp.dir/sccoda.dir/sccoda_input_sample_meta.tsv.gz",
                    type=str, help='Sample metadata table indexed by sample_id.')
parser.add_argument("--conditionvar", default=None, type=str,
                    help='The sample_meta variable to model.')
parser.add_argument("--conditionvar_level", default=None, type=str,
                    help='Level of the sample conditionvar to use as reference.')
parser.add_argument("--outfile",
                    default="./diffcomp.dir/sccoda.dir/sccoda_stats.tsv.gz",
                    type=str, help='Suffix for scCODA outputs.')

args = parser.parse_args()

L.info("args:")
print(args)

# ############################################################################ #
# ################################ scCODA #################################### #
# ############################################################################ #

if not os.path.exists(args.cell_count_file):
  sys.exit("Input cell count file non-existing.")
  
cell_counts = pd.read_csv(args.cell_count_file, sep='\t')
L.info(cell_counts)

adata = dat.from_pandas(cell_counts, covariate_columns=["sample_id"])
L.info(adata)

if not os.path.exists(args.sample_meta):
  sys.exit("Input sample metadata file non-existing.")
  
sample_meta = pd.read_csv(args.sample_meta, sep='\t')
L.info(sample_meta)

adata.obs = sample_meta
L.info(adata)

# Barplot per sample
viz.stacked_barplot(adata, feature_name="sample_id")
plt.savefig(args.outfile.replace(".tsv.gz", "_pct_cell_barplot_persample.pdf"))
# Barplot per condition
viz.stacked_barplot(adata, feature_name=args.conditionvar)
plt.savefig(args.outfile.replace(".tsv.gz", "_pct_cell_barplot_percondition.pdf"))
# Boxplots per condition
viz.boxplots(
    adata,
    feature_name=args.conditionvar,
    plot_facets=False,
    y_scale="relative",
    add_dots=False,
)
plt.savefig(args.outfile.replace(".tsv.gz", "_pct_cell_boxplot_percondition.pdf"))
# Cell cluster representation vs dispersion
viz.rel_abundance_dispersion_plot(
    data=adata,
    abundant_threshold=0.7
)
plt.savefig(args.outfile.replace(".tsv.gz", "_cell_scatter_prescenceVsdispersion.pdf"))

# Modelling
form = "C(" + args.conditionvar + ", Treatment('" + args.conditionvar_level + "'))"
L.info(form)

# Run with each cell type as the reference
cell_types = adata.var.index

for ct in cell_types:
    L.info(f"Reference: {ct}")

    # Run inference
    model = mod.CompositionalAnalysis(adata, formula=form,
                                      reference_cell_type=ct)
                                      
    results = model.sample_hmc(num_results=20000)
    
    res = results.effect_df.copy()
    cre = results.credible_effects().to_frame(name='credible').copy()
    res['credible'] = cre['credible'].copy()
    L.info(res)
    
    ct_name = ct.replace("/", "_")
    # save result stats as table
    res.to_csv(args.outfile.replace(".tsv.gz", "_ref_"+ct_name+".tsv.gz"), 
               sep='\t', compression='gzip')
    
    # saving all model results as pickle
    results.save(args.outfile.replace(".tsv.gz", "_ref_"+ct_name+".pickle"))

L.info("scCODA done.")
