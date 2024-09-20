#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Predict cell abundancies')
parser.add_argument("visium", type=str,default=None,help='visium h5ad file')
parser.add_argument("ref", type=str,default=None,help='reference signatures in csv format')
parser.add_argument("output", type=str,default=None,help='folder to write output')
parser.add_argument("--batch_key", type=str,default=None,help='column in adata.obs to be used as batch')
parser.add_argument("--detection_alpha", type=float,default=None,help='value of alpha parameter')
parser.add_argument("--N_cells_per_location", type=int,default=None,help='value of alpha parameter')
parser.add_argument("--max_epochs", type=int,default=50000,help='number of epochs')
parser.add_argument("--batch_size", type=int,default=None,help='number of train batches (use with caution, non None values slows training and may produce weird results)')
parser.add_argument("--do_not_filter_empty", action='store_true',help='by defauls all non tissue spots (adata.obs.in_tissue==0) will be removed. This flag switchs removal off.')
parser.add_argument("--predict_gene_exp", action='store_true',help='Predict per celltype*spot gene expression. Be carefull, it results in extrmely large output h5ad.')
parser.add_argument("--seed", type=int,default=1,help='scvi seed value')

import sys
import os
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import torch
import cell2location
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel

import scvi

args = parser.parse_args()

################
os.mkdir(args.output)
sys.stdout = open(args.output+"/c2l.pred.log", "w")
print(args)
print("cuda avaliable: "+str(torch.cuda.is_available()))
scvi.settings.seed = args.seed

vis = sc.read_h5ad(args.visium)
inf_aver = pd.read_csv(args.ref,index_col=0)

if (inf_aver.index.value_counts()>1).sum() > 0:
  raise Exception("Reference gene names are not unique, please fix it in '"+args.ref+"'") 

if (vis.var_names.value_counts()>1).sum() > 0:
  raise Exception("Spatial gene names are not unique, please fix it in '"+args.ref+"'") 


intersect = np.intersect1d(vis.var_names, inf_aver.index)


if (not args.do_not_filter_empty) & ('in_tissue' in vis.obs):
    vis = vis[vis.obs.in_tissue.astype('str').astype('int')==1,]

vis = vis[:, intersect].copy()

inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=vis, 
                                                batch_key=args.batch_key)


mod = cell2location.models.Cell2location(
    vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=args.N_cells_per_location,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=args.detection_alpha
)
mod.view_anndata_setup()

mod.train(max_epochs=args.max_epochs, 
          batch_size=args.batch_size,
          train_size=1,
          scale_elbo=1,
          accelerator='gpu',
          #use_gpu=True, this option is outdated
          progress_bar_refresh_rate=0)

vis = mod.export_posterior(
    vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# plot ELBO loss history during training, removing first 100 epochs from the plot
fig, ax = plt.subplots()
mod.plot_history(1000)
fig.legend(labels=['full data training']);
fig.savefig(args.output+'/train.history.pdf') 

if args.predict_gene_exp:
    # Compute expected expression per cell type
    expected_dict = mod.module.model.compute_expected_per_cell_type(
        mod.samples["post_sample_q05"], mod.adata_manager
    )
    # Add to anndata layers
    for i, n in enumerate(mod.factor_names_):
        vis.layers[n] = expected_dict['mu'][i]


# Save model
mod.save(args.output+"/predmodel", overwrite=True)
# most likely I do not need it. 
# plus it can fail because of unallowed celltype names (slashes)
try:
    vis.write(args.output+"/predmodel/sp.h5ad")
except Exception as e:
    print(e)


for k in  vis.obsm_keys():
    if  'cell_abundance' in k:
        vis.obsm[k].to_csv(args.output+"/predmodel/" + k + '.csv')

mod.plot_QC()
plt.savefig(args.output+'/predict.QC.pdf')

cell2location.utils.list_imported_modules()
sys.stdout.close()



