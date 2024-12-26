#!/bin/bash

# 执行每条语句前打印时间
printExecuteTime() {
    echo -e "\n-------------------------" start running $*: $(date "+%Y-%m-%d %H:%M:%S") "-------------------------"
}

#################################### SVD ############################################

printExecuteTime singular_value_decomposition.R
Rscript singular_value_decomposition.R 

printExecuteTime sample_loading_plot.R
Rscript sample_loading_plot.R

printExecuteTime s_value.R
Rscript s_value.R 

printExecuteTime gene_correlation.R
Rscript gene_correlation.R 

#################################### CCF&timing ##################################

printExecuteTime CCF_histone.R
Rscript CCF_histone.R 

printExecuteTime CCF_histone_perturbation.R
Rscript CCF_histone_perturbation.R 

printExecuteTime CCF_metabolome.R
Rscript CCF_metabolome.R 

printExecuteTime CCF_metabolome_perturbation.R
Rscript CCF_metabolome_perturbation.R 


#################################### contribution ##################################

printExecuteTime histone_contribution.R
Rscript histone_contribution.R 


#################################### enrichment ##################################

printExecuteTime gene_enrichment_analysis.R
Rscript gene_enrichment_analysis.R 

printExecuteTime format_enrichment_result.R
Rscript format_enrichment_result.R 


################################################################################

echo -e "\n-------------------------"end $*: $(date "+%Y-%m-%d %H:%M:%S")"-------------------------"

