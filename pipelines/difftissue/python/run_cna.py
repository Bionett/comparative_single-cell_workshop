# ------------------------------------------------------------------------------
# Run Covarying Neighborhood Analysis (CNA)
# ------------------------------------------------------------------------------

import os
import sys
import logging
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
from multianndata import MultiAnnData
import cna

# ############################################################################ #
# ###################### Set up the logging ################################## #
# ############################################################################ #

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
L = logging.getLogger("run_cna.py")

# ############################################################################ #
# ######################## Parse the arguments ############################### #
# ############################################################################ #

parser = argparse.ArgumentParser()
parser.add_argument("--adata_file", default=None, type=str, 
                    help='Scanpy like *.h5ad file with a cell-cell graph.')
parser.add_argument("--conditionvar", default=None, type=str,
                    help='The sample metadata variable to model.')
parser.add_argument("--covariate",
                    default=None, type=str, 
                    help='Sample metadata variable with 2 levels only.')
parser.add_argument("--batch",
                    default=None, type=str, 
                    help='Library construction batch (e.g. 10X well).')                    
parser.add_argument("--outfile",
                    default="./diffcomp.dir/cna.dir/cna_stats.tsv.gz",
                    type=str, help='Suffix for CNA outputs.')

args = parser.parse_args()

L.info("args:")
print(args)

# ############################################################################ #
# ########################## Read adata object ############################### #
# ############################################################################ #

if not os.path.exists(args.adata_file):
    sys.exit("Input h5ad file non-existing.")

adata = ad.read_h5ad(args.adata_file)
L.info("AnnDara")
L.info(adata)

adata = MultiAnnData(adata, sampleid='sample_id')
L.info("AnnDara multsample")
L.info(adata)

# ############################################################################ #
# ######################## Setting sample metadata ########################### #
# ############################################################################ #

L.info("Sample metadata")
L.info(adata.samplem)

samplevars = [args.conditionvar]

if not args.covariate == 'None':
  samplevars.append(args.covariate)

L.info(samplevars)

new_meta = []
for cv in samplevars:
  nmeta = pd.get_dummies(adata.obs[cv])
  if cv == args.conditionvar:
    varlevels = nmeta.columns
  new_meta.append(nmeta)

meta = pd.concat(new_meta)
L.info(meta)

if not args.batch == None:
  batch = pd.factorize(adata.obs[args.batch])[0].copy()
  batch = np.array(batch, dtype=np.int64)
  L.info(batch)
  meta['batch'] = batch
  L.info(meta)

meta['sample_id'] = adata.obs.loc[meta.index, 'sample_id']

L.info(meta)

adata.obs = meta.loc[adata.obs.index, :]
adata.obs_to_sample(meta.columns[:-1])

for col in adata.samplem.columns:
  adata.samplem[col] = np.array(adata.samplem[col], dtype=np.int64)

L.info(adata.samplem.head())

# ############################################################################ #
# ######################## Checking cell graph ############################### #
# ############################################################################ #

graph_keys = list(adata.obsp.keys())

for k in graph_keys:
  if 'connectivities' in k:
    adata.obsp['connectivities'] = adata.obsp[k]
  if 'distances' in k:
    adata.obsp['distances'] = adata.obsp[k]

L.info(adata.obsp)

# ############################################################################ #
# ######### Running CNA association per condition-variable level ############# #
# ############################################################################ #

for vlevel in varlevels:
  L.info(vlevel)
  if args.batch == None:
    res = cna.tl.association(adata,
                             adata.samplem[vlevel])
  else:
    res = cna.tl.association(adata,
                             adata.samplem[vlevel],
                             batches=adata.samplem.batch)
  L.info('global association p-value:', res.p)
  res.p.tofile(args.outfile.replace(".tsv.gz","_"+vlevel+"_global_pval.tsv.gz"))
  # Plot variance explained by each NAM PC
  plt.plot(adata.uns['NAM_svs']/adata.uns['NAM_svs'].sum(), marker='o', 
           linestyle='--')
  plt.xlabel('NAM PC')
  plt.ylabel('fraction of variance explained')
  plt.savefig(args.outfile.replace(".tsv.gz", "_"+vlevel+"_nam_pc_var.pdf"))
  # Save NAM principal components
  nam_pc = adata.uns['NAM_sampleXpc'].copy()
  nam_pc.to_csv(args.outfile.replace(".tsv.gz", "_"+vlevel+"_nam_pcs.tsv.gz"), 
                sep="\t")
  # Cell/neighbor association (corr) with condition-var level 
  cna.pl.umap_ncorr(adata,                           
              res,                               
              scatter0={'alpha':0.5, 's':20},    
              scatter1={'alpha':0.05, 's':20})
  plt.title('p = {:.2e}'.format(res.p))
  plt.savefig(args.outfile.replace(".tsv.gz", "_"+vlevel+"_corr_umap.png"))
  # 
  corr = pd.DataFrame(res.ncorrs)
  corr['barcode_id']=adata.uns['NAM_nbhdXpc'].index.copy()
  corr.to_csv(args.outfile.replace(".tsv.gz", "_"+vlevel+"_cell_corr.tsv.gz"), 
              sep="\t")
  # compute correlation between neighborhood coefficients and genes
  adata.var['corr_case'] = np.corrcoef(res.ncorrs.reshape(1,-1), 
                                       adata[adata.uns['NAM_nbhdXpc'].index, ].X, 
                                       rowvar=False)[0,1:]
  # Save genes strongly correlated genes
  gs = adata.var.sort_values('corr_case', ascending=False).copy()
  gs['gene'] = gs.index
  gs[['gene', 'corr_case']].to_csv(args.outfile.replace(".tsv.gz", 
                                                        "_"+vlevel+"_gene_corr.tsv.gz"), 
                                                        sep="\t")
                                                        
L.info("CNA done.")
