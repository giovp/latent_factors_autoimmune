OpenMPController::omp_set_num_threads(snakemake@threads)

library(dplyr)
library(CountClust)
library(SingleCellExperiment)
#snakemake <- readRDS("./RAsubset_pipeline/output/LOGS_AND_BENCHMARKS/lda_snakemake_55_32.rds")
logfile = gsub(".log",".rds",snakemake@log[["lda_log"]])
saveRDS(snakemake, logfile)
#below line is to load the snakemake for debugging purposes
#snakemake = readRDS("")

#input
scset_input = snakemake@input[["SCE"]]
k = as.integer(snakemake@wildcards[["k"]])
seeed = as.integer(snakemake@wildcards[["seed"]])
cores = as.integer(snakemake@threads)
n_trials = as.integer(snakemake@params[["n_trials"]])

#output
lda_output = snakemake@output[["file_output"]]

scset <- readRDS(scset_input)

set.seed(seeed)
lda_obj <- FitGoM(t(assays(scset)$counts),
                  options="BF",
                  K=k, tol=10, num_trials = n_trials)

compGoM_obj <- compGoM(t(counts(scset)), lda_obj$fit)

lda_obj$BIC <- compGoM_obj$BIC
lda_obj$loglik <- compGoM_obj$loglik
lda_obj$null_loglik<- compGoM_obj$null_loglik
lda_obj$seed <- seeed
lda_obj$k <- k


saveRDS(lda_obj, lda_output)