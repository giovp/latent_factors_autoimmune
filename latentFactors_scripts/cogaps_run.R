OpenMPController::omp_set_num_threads(snakemake@threads)

library(dplyr)
library(CoGAPS)
library(SingleCellExperiment)

logfile = gsub(".log",".rds",snakemake@log[["cogaps_log"]])
saveRDS(snakemake, logfile)
#below line is to load the snakemake for debugging purposes
#snakemake = readRDS("./latentFactors_pipeline/output/LOGS_AND_BENCHMARKS/cogaps_snakemake_1_3.rds")

#input
scset_input = snakemake@input[["SCE"]]
k = as.integer(snakemake@wildcards[["k"]])
seeed = as.integer(snakemake@wildcards[["seed"]])
cores = as.integer(snakemake@threads)
n_iter = as.integer(snakemake@params[["n_iter"]])

#output
cogas_output = snakemake@output[["file_output"]]

system(sprintf("taskset -p 0xffffffff %d", Sys.getpid()))
scset <- readRDS(scset_input)

params <- new("CogapsParams")
params <- setParam(params, "nPatterns", k)
params <- setParam(params, "nIterations", n_iter)
params <- setParam(params, "seed", seeed)
params <- setDistributedParams(params, nSets=cores)
#getParam(params, "nSets")

result <- scCoGAPS(assays(scset)$logcounts, 
                   params, 
                   distributed="single-cell", 
                   sparseOptimization=TRUE,
                   singleCell=TRUE,
                   BPPARAM=BiocParallel::MulticoreParam(workers=cores),
                   messages=T)

saveRDS(result, cogas_output)