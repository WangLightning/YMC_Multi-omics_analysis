suppressMessages({
    library(openxlsx)
    library(rlist)
    library(tidyverse)
})


# Parameters setting ----------------------------------------------------------

# read config file
config <- list.load("config.yaml")

dir_result <- "../results"
minCommonGeneNum <- config$enrichment$minCommonGeneNum
maxCommonGeneNum <- config$enrichment$maxCommonGeneNum
signif_level <- config$enrichment$signif_level
number_dim <- config$number_dim
data_list <- config$data_list


species_list <- c("yeast")
family_list <- c("S.cerevisiae")
species_family_list <- list(
    list(species = "yeast", family = "S.cerevisiae")
)
pathway_list <- c(
    "go.bp",
    "go.cc",
    "go.mf",
    "kegg",
    "reactome"
)
gene_list <- str_c("gene", 1:(number_dim-1))


# load functions --------------------------------------------------------------

source("enrichment_utils.R")


# transform data -------------------------------------------------------------------

for (data_name in data_list) {

    # Directory setting
    enrichment.dir <- file.path(dir_result, data_name, "enrichment")
    if (!dir.exists(enrichment.dir)) dir.create(enrichment.dir, recursive = TRUE)

    # read data
    raw_result.dir <- file.path(dir_result, "temp", data_name, "enrichment")
    enrichment_result <- read_raw_result(raw_result.dir)

    # result1
    outputCols <- c(
        "id",
        "pvalue.pos",
        "pvalue.neg",
        "description",
        "num.of.targets",
        "num.of.non.targets",
        "pathwayGeneNum",
        "pathwayGene_in_geneEigenvectors",
        "pathwayGene",
        NULL
    )
    bgFill <- list()
    bgFill$pos <- RColorBrewer::brewer.pal(9, "Reds")[c(2, 3, 4)] %>% rev()
    bgFill$neg <- RColorBrewer::brewer.pal(9, "Blues")[c(3, 4, 5)] %>% rev()
    for (s in species_list) {
        for (g in gene_list) {
            result <- formatResult1(enrichment_result, s, g)
            saveResult1(result, s, g)
        }
    }


    # result2
    outputCols <- c(
        "database",
        "id",
        "pvalue",
        "description",
        "num.of.targets",
        "num.of.non.targets",
        "pathwayGeneNum",
        "pathwayGene_in_geneEigenvectors",
        "pathwayGene",
        NULL
    )
    for (s in species_list) {
        for (g in gene_list) {
            result <- formatResult2(enrichment_result, s, g)
            saveResult2(result, s, g)
        }
    }

    # result 3
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

