#!/bin/bash

# nohup bash run.sh > run.log 2>&1 &

set -ue

# 执行每条语句前打印时间
printExecuteTime() {
    echo -e "\n-------------------------" start running $*: $(date "+%Y-%m-%d %H:%M:%S") "-------------------------"
}

#################################### main ############################################

printExecuteTime data_process.R
Rscript data_process.R 

printExecuteTime matrix_decomposition.R
Rscript matrix_decomposition.R 

printExecuteTime matrix_decomposition_metabolic.R
Rscript matrix_decomposition_metabolic.R 

printExecuteTime gene_enrichment_analysis.R
Rscript gene_enrichment_analysis.R 

printExecuteTime format_enrichment_result.R
Rscript format_enrichment_result.R 


################################################################################

echo -e "\n-------------------------"end $*: $(date "+%Y-%m-%d %H:%M:%S")"-------------------------"

