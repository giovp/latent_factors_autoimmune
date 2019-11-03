# Snakemake implementation of latent factor models comparison

This pipeline was implemented to streamline the analysis of several latent factor models, their comparison and their gene set coverage.  
The models implemented are:
* **scCoGAPS** <a href="https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30146-2">paper</a> - [code](https://github.com/FertigLab/CoGAPS)  
* **LDA** [paper](https://www.biorxiv.org/content/10.1101/461228v2) - [code](https://github.com/kkdey/CountClust)  
* **scHPF** [paper](http://msb.embopress.org/content/15/2/e8557) - [code](https://github.com/simslab/scHPF)  
* **scVI** [paper](https://www.nature.com/articles/s41592-018-0229-2) - [linear-decoder paper](https://www.biorxiv.org/content/10.1101/737601v1) - [code](https://github.com/YosefLab/scVI)

To use this pipeline, you need to:
* Have Python > 3.6
* Have Conda: https://docs.conda.io/en/latest/miniconda.html
* Have Snakemake installed `conda install -c bioconda -c conda-forge snakemake`
* Have R>3.5 in your system. With packages:
  * `tidyverse`
  * `data.table`
  * `feather` (used to have an easy format for both R and python)
  * `reticulate`   
    *Plus the Method specific packages*
  * `CountClust` [github](https://github.com/kkdey/CountClust) and make sure to install `maptpx` also from github
  * `CoGAPS` [github](https://github.com/FertigLab/CoGAPS)
  * `scVI` [github](https://github.com/YosefLab/scVI). The Tutorial for the linear decoder can be found [here](https://nbviewer.jupyter.org/github/YosefLab/scVI/blob/master/tests/notebooks/linear_decoder.ipynb)
  * `HPF` [github](https://github.com/simslab/scHPF)

### Folders
To have a smooth start with the pipeline, create a new folder anywhere with the name of the analysis and within this folder, create `output` folder and change the paths in the `config.json` ad `cluster.json` files accordingly. Make sure that all the paths are consistent and corrects, this is also true for the `latentFactors_scripts` folder, that contains the environments specifics files.

### Files
You need several specific files to run this analysis. All files needed are listed in the `input_files` section of the `config.json` file. Here's the specifics:
* **counts**: dataframe cells (rows) v. genes (columns) of INTEGER counts. Because of method specificities, they cannot be log counts or CPM or any other format. The dataframe has to be saved in feather format
* **metadata**: the `SingleCellExperiment@colData` dataframe (in feather format). Here, you need couple of specific columns:
  * `cell_name` column, with the sample ID of the cells
  * `plate` column, in case you have batches (this is mostly needed by *scVI*
* **rowdata**: the `rowData(SingleCellExperiment)` dataframe (in feather format). Here, you need couple of specific columns:
  * `ENTREZID` column, with *entrezID* for the genes
  * `symbol` column, with *symbolsID* for the genes
* **SCE**: this is the `SingleCellExperiment` object. Here, you need to have bot counts and log counts slots, again because of specifics of the different methods.

Finally, in the `cluster.json` file there are the cluster specs for the job submission.

### Parameters
In the `config.json` file you will find the parameters of the analysis, change them accordingly to your needs.  

### Run
The pipeline requires [Snakemake](https://snakemake.readthedocs.io/en/stable/). Configuration files for parameters and cluster settings (in our case, SGE) are provided.

#### Additional files
The hetnets were pre-computed for the gene collections of interest. You might want to modify the gene set collections or generate new yourself. In that case, check [this amazing tutorial](https://github.com/greenelab/BioBombe/tree/master/3.build-hetnets) by the original authors!
