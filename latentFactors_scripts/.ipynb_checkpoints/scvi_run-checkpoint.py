import os
import re
import json
import anndata
import numpy as np
import pandas as pd
import feather
import torch
from scvi.dataset import GeneExpressionDataset, AnnDataset
from scvi.models import *
from scvi.inference import UnsupervisedTrainer
from scvi.inference.posterior import *

save_snakemake = re.sub(".log", ".pkl", str(snakemake.log["scvi_log"]))
fileobj = open(save_snakemake, 'wb')
pickle.dump(snakemake, fileobj)
fileobj.close()

#input    
input_rawcounts = snakemake.input["counts"]
input_meta = snakemake.input["metadata"]
input_var = snakemake.input["rowdata"]
#output
output_metrics = snakemake.output["metrics"]
output_model = snakemake.output["model"]
output_factors = snakemake.output["factors"]
output_loadings = snakemake.output["loadings"]
#params snakemake
seed = int(snakemake.wildcards["seed"])
k = int(snakemake.wildcards["k"])
batch_id = snakemake.params["batch_id"]
#seed
torch.manual_seed(seed)
np.random.seed(seed)

#params vae
size = 0.8
lr = 5e-4
latent = 16
layer = 1
hidden = 128
n_epochs_all = None
#params trainer
n_epochs=600 if n_epochs_all is None else n_epochs_all
dispersion="gene"
use_labels=False
use_cuda=False
reconstruction_loss="nb"


rawcounts = feather.read_dataframe(input_rawcounts)
meta = feather.read_dataframe(input_meta)
meta.index =  meta.loc[:,"cell_name"].values.astype(str)
var = feather.read_dataframe(input_var)
var.index = var.loc[:,"symbol"].values.astype(str)

annobj = anndata.AnnData(X=rawcounts)
annobj.obs = meta
annobj.var = var

X, local_mean, local_var, batch_indices, labels = GeneExpressionDataset.get_attributes_from_matrix(annobj.X)

geneExp = GeneExpressionDataset(X, local_mean, local_var, batch_indices, labels, gene_names=annobj.var.index)

if bool(batch_id) is not False:
    use_batches=True
    plates, plates_ids = pd.factorize(annobj.obs[batch_id])
    geneExp.batch_indices = plates.reshape(-1, 1)
    geneExp.n_batches = np.unique(plates.reshape(-1, 1)).size
else:
    use_batches = False

ldvae = LDVAE(geneExp.nb_genes, 
            n_batch=geneExp.n_batches * use_batches, 
            n_latent=latent,
            n_layers=layer,
            n_hidden=hidden,
            dispersion=dispersion,
            reconstruction_loss=reconstruction_loss)
            
trainer = UnsupervisedTrainer(ldvae,
                              geneExp,
                              train_size=size,
                              use_cuda=use_cuda,
                              frequency=5)

trainer.train(n_epochs=n_epochs, lr=lr)
                    
history = pd.DataFrame.from_dict(trainer.history).tail(50).reset_index()
history["train_size"] = size
history["learning_rate"] = lr
history["n_latent"] = latent
history["n_layers"] = layer
history["n_hidden"] = hidden
history["seed"] = seed

loadings = pd.DataFrame(list(ldvae.decoder.factor_regressor.parameters())[0].detach().numpy())
loadings.columns = ["Loading"+str(i) for i in np.arange(1,loadings.shape[1]+1)]

posterior = trainer.create_posterior(trainer.model, geneExp, indices=np.arange(len(geneExp)))
latent_posterior, batch_indices, labels = posterior.sequential().get_latent()
factors = pd.DataFrame(latent_posterior)
factors.columns = ["Factor"+str(i) for i in np.arange(1,loadings.shape[1]+1)]

#save
loadings.to_feather(output_loadings)
factors.to_feather(output_factors)
history.to_feather(output_metrics)
torch.save(trainer.model.state_dict(), output_model)





