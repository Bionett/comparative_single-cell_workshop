# ------------------------------------------------------------------------------
# Build a pbulk count matrix based on two cell metavariables (cell-type & sample) 
# save .tsv.gz matrix.
# ------------------------------------------------------------------------------

import os
import sys
import logging
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
import warnings
import scanpy as sc

# ############################################################################ #
# ###################### Set up the logging ################################## #
# ############################################################################ #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("build_cellxgene_public.py")

# ############################################################################ #
# ######################## Parse the arguments ############################### #
# ############################################################################ #

parser = argparse.ArgumentParser()
parser.add_argument("--adata_file", default=None, type=str,
                    help='anndata file with count layer.')
parser.add_argument("--cell_metadata", default="cell_metadata.tsv.gz", type=str,
                     help='A table indexed by `barcode_id`.')
parser.add_argument("--cell_var", default=None, type=str,
                     help='Cell metadata variable to denote the cell-type.')
parser.add_argument("--sample_var", default=None, type=str,
                     help='Cell metadata variable representing sample/library.')
parser.add_argument("--outfile", default="output_pbulk.tsv.gz", type=str,
                     help='*.tsv.gz output file.')

args = parser.parse_args()

L.info("args:")
print(args)

# ############################################################################ #
# ####################### Read h5ad + cell-meta ############################## #
# ############################################################################ #

if not os.path.exists(args.adata_file):
  sys.exit("Input anndata non-existing.")

adata = sc.read_h5ad(args.adata_file)
L.info(adata)

if not os.path.exists(args.cell_metadata):
  sys.exit("Input cell-metadata file non-existing.")

nobs = pd.read_csv(args.cell_metadata, sep = '\t', low_memory = False, 
                   index_col= "barcode_id")

adata.obs = nobs.loc[adata.obs.index]

adata.obs["sample_cell"] = adata.obs[args.sample_var].str.cat(others=adata.obs[args.cell_var], sep = "_")

L.info(adata.obs.columns)

sample_cells = adata.obs.sample_cell.unique()

pbulk = np.empty(shape=(adata.shape[1], sample_cells.shape[0]))
L.info(pbulk)

# ############################################################################ #
# ################################# pseudo-bulk ############################## #
# ############################################################################ #

i = 0
for sc in sample_cells:
    L.info(sc)
    sadata = adata[adata.obs["sample_cell"] == sc].copy()
    L.info(sadata.X.shape)
    pbulk[:,i] = sadata.layers['counts'].sum(0)
    i = i+1 
    
pbulk = pd.DataFrame(pbulk)
pbulk.columns = sample_cells
pbulk.index = adata.var_names

L.info(pbulk)

pbulk.to_csv(args.outfile)

L.info(args.outfile)
L.info("get_pbulk.py done.")
