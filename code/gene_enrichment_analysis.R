suppressMessages({
    library(rlist)
    library(tidyverse)
    library(snowfall)
})


# Parameters setting ----------------------------------------------------------

# read config file
config <- list.load("config.yaml")

dir_result <- "../results"
database.dir <- "../data/pathway_database_20230530"
CPUs <- config$CPUs
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


# read pathway data -------------------------------------------------------------------

pathway_data.list <- list()
for (f in family_list) {
    for (p in pathway_list) {
        pathway_data.list[[f]][[p]] <- read_tsv(
            file.path(
                database.dir,
                str_glue("{f}.{p}.tsv")
            ),
            col_types = cols(id = col_character())
        )
    }
}


# pathway data pre-processing --------------------------------------------------------------

for (f in family_list) {

    # GO
    minGeneNum <- config$enrichment$GO$minGeneNum
    maxGeneNum <- config$enrichment$GO$maxGeneNum
    for (p in c("go.bp", "go.cc", "go.mf")) {
        pathway_data.list[[f]][[p]] <- pathway_data.list[[f]][[p]] %>%
            filter(geneNum %>% between(minGeneNum, maxGeneNum))
    }

    # KEGG
    p <- "kegg"
    minGeneNum <- config$enrichment$KEGG$minGeneNum
    maxGeneNum <- config$enrichment$KEGG$maxGeneNum
    pathway_data.list[[f]][[p]] <- pathway_data.list[[f]][[p]] %>%
        filter(geneNum %>% between(minGeneNum, maxGeneNum))

    # Reactome
    p <- "reactome"
    minGeneNum <- config$enrichment$Reactome$minGeneNum
    maxGeneNum <- config$enrichment$Reactome$maxGeneNum
    pathway_data.list[[f]][[p]] <- pathway_data.list[[f]][[p]] %>%
        filter(geneNum %>% between(minGeneNum, maxGeneNum))

}


# enrichment ----------------------------------------------------------

outputCols <- c(
    "id",
    "pvalue.pos",
    "pvalue.neg",
    "description",
    "num.of.targets",
    "num.of.non.targets",
    "pathwayGeneNum",
    "pathwayGene_in_geneEigenvectors",
    "pathwayGene"
)

for (data_name in data_list) {

    cat("Running Enrichment for", data_name, "\n")

    par.list <- expand_grid(
        species_family_list %>% bind_rows(),
        gene = gene_list,
    ) %>%
        # transpose() %>%
        list.parse() %>%
        list.names(str_glue("{species} {gene}"))

    loading.dir <- file.path(dir_result, "temp", data_name, "loadings")

    sfInit(parallel = TRUE, cpus = CPUs)
    sfLibrary(openxlsx)
    sfLibrary(tidyverse)
    sfExportAll()
    result <- sfLapply(
        par.list,
        do_enrichment,
        loading.dir
    )
    sfStop()


    # transform result 
    enrichment_result <- list()
    for (s in species_list) {
        for (g in gene_list) {
            enrichment_result[[s]][[g]] <- result[[str_glue("{s} {g}")]]
        }
    }

    raw.result.dir <- file.path(dir_result, "temp", data_name, "enrichment")
    if (!dir.exists(raw.result.dir)) dir.create(raw.result.dir, recursive = TRUE)
    # write data 
    write_raw_result(enrichment_result, raw.result.dir)

}












