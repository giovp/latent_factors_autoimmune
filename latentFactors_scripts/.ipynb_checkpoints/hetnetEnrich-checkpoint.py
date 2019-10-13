"""
Adapted from https://github.com/greenelab/BioBombe/blob/master/6.biobombe-projection/geneset_tracking.py

Apply BioBombe network projection approach to loading matrices to decipher gene signatures 
associated to latent factors
        
Output:
A single long dataframe storing the values and z_scores of all input gene set
analysis results for the given input dataset across one weight matrix for
all bottleneck dimensions and all compression algorithms.

Note:
Reason for feather as input/output is purely due to inability to use reticulate/rpy2 on the cluster/conda-forge channel, see:
    https://github.com/conda-forge/conda-forge.github.io/issues/234
    https://github.com/rstudio/reticulate/issues/428
"""
import re
import os
import sys
import pandas as pd
import feather
import biobombe_latent
import argparse
import glob
import numpy as np
import pickle

#save snakemake object
save_snakemake = re.sub(".log", ".pkl", str(snakemake.log["hetnet_log"]))
fileobj = open(save_snakemake, 'wb')
pickle.dump(snakemake, fileobj)
fileobj.close()

#import variables and params
gene_ids = snakemake.input["gene_ids"]
loadings_input = snakemake.input["loadings_best"]
path_hetnets = snakemake.params["path_hetnets"]
metaedge = snakemake.params["metaedge"]
shuffled = False
geneset_output = snakemake.output["hetnet_outputs"]
#read in gene table
rowData = feather.read_dataframe(gene_ids)
gene_list = rowData[["ENTREZID","symbol"]].dropna()
loadings_matrix = feather.read_dataframe(loadings_input)

loadings_matrix = loadings_matrix.set_index(gene_list["ENTREZID"].values.astype(str))
#takes long to load all hetnets
print("Loading hetnets...takes long")
hetnets = biobombe_latent.load_hetnets(
    hetnet_file = path_hetnets+'/interpret_hetnet_wikipath.json.bz2',
    permuted_directory = path_hetnets +'/permuted/',
    subset_genes = gene_list["ENTREZID"].values.astype(str),
    metaedge_abbrev = metaedge
   )

# Extract results from the weight weight matrices
weight_matrix_seed_results = []

# Build output file
if shuffled:
    out_base = 'geneset_scores_shuffled.feather'
else:
    out_base = 'geneset_scores.feather'
    
weight_matrix = loadings_matrix.T
    
mult_results = {}
all_list = []
weight_matrix_seed_results = []
for model in hetnets.keys():
    print('        and model {}'.format(model))
    hetnet = hetnets[model]
    mult_results[model] = weight_matrix @ hetnet 
    long_result = (
        mult_results[model]
        .reset_index()
        .melt(id_vars=['index'])
    )
    long_result = long_result.assign(model=model)
        
    if model != 'real':
        long_result = long_result.assign(model_type='permuted')
    else:
        long_result = long_result.assign(model_type='real')
            
    all_list.append(long_result)

all_df = pd.concat(all_list)
all_df.value = all_df.value.astype(float)

    # Obtain the z-scores for the real data against permutations
permuted_mean = (
    all_df
    .groupby(['model_type', 'index', 'variable'])
    .mean()
    .reset_index()
    .query("model_type == 'permuted'")
)

permuted_std = (
    all_df
    .groupby(['model_type', 'index', 'variable'])
    .std()
    .reset_index()
    .query("model_type == 'permuted'")
)

real_df = (
    all_df
    .groupby(['model_type', 'index', 'variable'])
    .mean()
    .reset_index()
    .query("model_type == 'real'")
)

z_score = (
    real_df.reset_index(drop=True).value - permuted_mean.value
) / permuted_std.value
    
real_df = real_df.reset_index(drop=True).assign(z_score=z_score)

# Process matrix for output
algorithm_feature_df = real_df['index'].str.split('_', expand=True)
algorithm_feature_df.columns = ['feature','nIter','kLatent']

real_df = pd.concat([real_df, algorithm_feature_df], axis='columns')
    #real_df = real_df.drop(['index'], axis='columns')
    #real_df = real_df.assign(
    #    z=z_dim,
    #    seed='seed'
    #)
weight_matrix_seed_results.append(real_df)
weight_seed_df = pd.concat(weight_matrix_seed_results)
weight_seed_df.value = weight_seed_df.value.astype(float)
weight_seed_df = weight_seed_df.sort_values(by='z_score')
feather.write_dataframe(weight_seed_df, geneset_output)