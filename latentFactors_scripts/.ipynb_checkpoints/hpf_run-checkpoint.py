import pickle
import re
import pandas as pd
import numpy as np
import scipy.io, scipy.sparse
import feather
import schpf as hpf
from sklearn.externals import joblib

save_snakemake = re.sub(".log", ".pkl", str(snakemake.log["hpf_log"]))
fileobj = open(save_snakemake, 'wb')
pickle.dump(snakemake, fileobj)
fileobj.close()

#input_variables
counts_input = snakemake.input["counts"]
#output_variables
cell_score_output = snakemake.output["cell_score"]
gene_score_output = snakemake.output["gene_score"]
maximum_overlap_output = snakemake.output["maximum_overlaps"]
ranked_genes_output = snakemake.output["ranked_genes"]
model_output = snakemake.output["model_output"]
metrics_output = snakemake.output["metrics_output"]
#params
seed = int(snakemake.wildcards["seed"])
k = int(snakemake.wildcards["k"])
n_trials = int(snakemake.params["n_trials"])
#set seed
np.random.seed(seed)

#read file
df = feather.read_dataframe(counts_input)
#start analysis
sparse_arr = scipy.sparse.coo_matrix(df.to_numpy())
model = hpf.run_trials(sparse_arr, nfactors=k, ntrials=n_trials, validation_data = None)

metrics = pd.DataFrame()
metrics["loss"] = model.loss[-10:]
metrics["k"] = np.repeat(k,10)
metrics["seed"] = np.repeat(seed, 10)

cell_score = model.cell_score()
gene_score = model.gene_score()
table = hpf.max_pairwise_table(gene_score, ntop_list=[50,100,150,200,250,300])


genes = df.iloc[:,[1]].to_numpy(dtype=str)
name_col = 0
ranks = np.argsort(gene_score, axis=0)[::-1]
ranked_genes = []
for i in range(gene_score.shape[1]):
    ranked_genes.append(genes[ranks[:,i], name_col])
ranked_genes = np.stack(ranked_genes).T

#transform and save
cell_score_df = pd.DataFrame(cell_score)
cell_score_df.columns = ["Factor"+str(i) for i in np.arange(1,cell_score.shape[1]+1)]
cell_score_df.to_feather(cell_score_output)

gene_score_df = pd.DataFrame(gene_score)
gene_score_df.columns = ["Factor"+str(i) for i in np.arange(1,gene_score.shape[1]+1)]
gene_score_df.to_feather(gene_score_output)

table.to_feather(maximum_overlap_output)
metrics.to_feather(metrics_output)
np.savetxt(ranked_genes_output, ranked_genes, fmt="%s", delimiter='\t')
joblib.dump(model, model_output)
