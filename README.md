# Introduction



# Directory description


## code folder

**singular_value_decomposition.R**  
This script reads the expression profiles of human and rat, performs singular value decomposition(SVD), and obtains the polarized gene and sample eigenvector. The polarized eigenvectors are saved in */results/loadings*

**gene_enrichment_analysis.R & format_enrichment_result.R**  
This script performs the gene set enrichment analysis by the Wilcoxon scoring method on the polarized gene eigenvector. The enrichment result is saved in */results/enrichment*

**rat_comparison.R**  
This script performs the comparison of GK/WST rat islet expressions of each week. The enrichment result is saved in */results/comparison*

**human_resampling.R**  
This script performs re-sampling of human islets to evaluate the stability of gene eigenvector. The result is saved in */results/human_resampling.csv*. It should be noted that the result run by CodeOcean is slightly different from the result file used to plot due to the different random seed.

**aggregation.R**  
This script performs Robust Rank Aggregation on the angiogenesis eigenvector of human and rat, and obtains the aggregated list. The result is saved in */results/aggregation.csv*

**concatenated_expression_profile_analysis.R**  
This script performs the analysis of the concatenated expression profile of human and rat. The result is saved in */results/concatenated_expression_profile_analysis*

**enrichment_utils.R**  
This script contains functions used to perform enrichment analysis.

**config.yaml**  
This is the configuration file that contains basic information about datasets and parameters.

## data folder

*/data/pathway* contains the pathway data of GO, KEGG, and Reactome. 

**sweeden_total_result_3.csv** and **reads_polished_reg_dp_fixed.csv** are the gene expression profile of human and rat. The rat RNA-seq sequencing reads of rats were downloaded from GEO under accession number GSE81811 and normalized by MUREN. The human Affymetrix microarray data were downloaded from GEO under accession number GSE38642 and normalized by sub-sub normalization.

**HuGene-1_0-st-v1.na36.hg19.probeset-single-name.csv** is the annotation file of the human microarray.

**rgd.gaf** and **goa_human.gaf** are the GO annotation files, which can be downloaded at <http://current.geneontology.org/products/pages/downloads.html>

**ensembl_convertID_102.tsv** is the annotation file of the conversion among Ensembl id, Entrez id, and gene symbol.

## Figures folder

This folder contains most of the figures in this article. Folders are arranged by the label of figures. We use R Notebook and Jupyter to display these figures. The figures of islet micrograph analysis, which include Figure 3C and Figure S8, are plotted by Python. Other figures in the folder are plotted by R.

