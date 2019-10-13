import numpy as np
import pandas as pd
import feather
import cca_core as cca

seeds = [23,46,83]
input_files = ["/da/ATI/bioinformatics/projects/ATI9EDA/2018-11-05_AMP_scRNAseq_RA_LN/latentFactors_pipeline/output/evaluate_output/lda_4_factors.feather", "/da/ATI/bioinformatics/projects/ATI9EDA/2018-11-05_AMP_scRNAseq_RA_LN/latentFactors_pipeline/output/evaluate_output/lda_6_factors.feather", "/da/ATI/bioinformatics/projects/ATI9EDA/2018-11-05_AMP_scRNAseq_RA_LN/latentFactors_pipeline/output/evaluate_output/lda_8_factors.feather"]
threshold = 0.98
for file in input_files:
    df = feather.read_dataframe(file)
    df_list = []
    for seed in seeds:
        idx_seed = ["23" in i for i in df.columns.values]
        df_subset = df.loc[:,idx_seed].copy()
        df_list.append(df_subset)
        
    result = cca.robust_cca_similarity(df_list[1].T, df_list[2].T, verbose=False, threshold=threshold)
    output_list.append(np.mean(result["cca_coef1"]))
    