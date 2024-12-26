suppressMessages({
    library(openxlsx)
    library(rlist)
    library(tidyverse)
})


# Parameters setting ----------------------------------------------------------

dir_result <- "../results"

# read config file
config <- list.load("config.yaml")

minCommonGeneNum <- config$enrichment$minCommonGeneNum
maxCommonGeneNum <- config$enrichment$maxCommonGeneNum
signif_level <- config$enrichment$signif_level
number_dim <- config$number_dim
data_list <- config$data_transcriptome_list

species_list <- c("yeast")
family_list <- c("S.cerevisiae")
species_family_list <- list(
    list(species = "yeast", family = "S.cerevisiae")
)
pathway_list <- c(
    "go.bp",
    "go.cc",
    "go.mf"
)
gene_list <- str_c("gene", 1:(number_dim-1))


# load functions --------------------------------------------------------------

source("enrichment_utils.R")


# transform data -------------------------------------------------------------------

for (data_name in data_list) {

    # Directory setting
    enrichment.dir <- file.path(dir_result, "enrichment", data_name)
    if (!dir.exists(enrichment.dir)) dir.create(enrichment.dir, recursive = TRUE)

    # read data
    raw_result.dir <- file.path(dir_result, "temp", data_name, "enrichment")
    enrichment_result <- read_raw_result(raw_result.dir)

    outputCols <- c(
        "database",
        "id",
        str_c("pvalue.", gene_list),
        "description",
        "pathwayGeneNum",
        "num.of.targets",
        "num.of.non.targets",
        "pathwayGene_in_geneEigenvectors",
        "pathwayGene",
        NULL
    )
    bgFill <- list()
    bgFill$pos <- RColorBrewer::brewer.pal(9, "Reds")[c(2, 4, 6)] %>% rev()
    bgFill$neg <- RColorBrewer::brewer.pal(9, "Blues")[c(2, 4, 6)] %>% rev()
    for (s in species_list) {
        result <- formatResult3(enrichment_result, s)
        saveResult3(result, s)
    }

}

