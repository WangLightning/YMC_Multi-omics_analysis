suppressMessages({
    library(rlist)
    library(openxlsx)
    library(tidyverse)
})


# parameters --------------------------------------------------------------

# read config file
config <- list.load("config.yaml")

dir_result <- "../results"
data_list <- c(
    'Metabolite_gc',
    'Metabolite_lc',
    'Metabolite_interpolated_gc',
    'Metabolite_interpolated_lc'
)
number_dim <- config$number_dim


# main ---------------------------------------------------------------


for (data_name in data_list) {

    # read data
    file_data  <- file.path("../data", str_glue("{data_name}.csv"))
    data <- read_csv(file_data, show_col_types = FALSE)

    name_gene <- data$Symbol
    name_sample <- names(data)[-1]
    data_mat <- data %>%
        select(-1) %>%
        as.matrix()

    # result directory
    dir_res <- file.path(dir_result, "temp", data_name)
    if (!dir.exists(dir_res)) dir.create(dir_res, recursive = TRUE)
    file_gene <- file.path(dir_res, "gene_name.txt")
    file_sample <- file.path(dir_res, "sample_name.txt")
    file_mat <- file.path(dir_res, "matrix.txt")

    # save result
    write.table(name_gene, file = file_gene, row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(name_sample, file = file_sample, row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(data_mat, file = file_mat, row.names = F, col.names = F, sep = "\t")

    # read data
    dir_temp <- file.path(dir_result, "temp", data_name)
    file_mat <- file.path(dir_temp, "matrix.txt")
    file_gene <- file.path(dir_temp, "gene_name.txt")
    file_sample <- file.path(dir_temp, "sample_name.txt")
    data_mat <- read_table(file_mat, col_names = F, show_col_types = FALSE) %>% as.matrix()
    data_gene <- read_table(file_gene, col_names = "Symbol", show_col_types = FALSE)
    data_sample <- read_table(file_sample, col_names = "Sample", show_col_types = FALSE)

    # singular value decomposition
    # genes*samples
    # u is for gene
    # v is for sample
    res <- svd(data_mat)
    s <- res$d
    u <- res$u[, 1:number_dim]
    v <- res$v[, 1:number_dim]


    if (data_name == "Metabolite_gc" ) {
        # u[,2] = -u[,2]
        # v[,2] = -v[,2]
        u[,3] = -u[,3]
        v[,3] = -v[,3]
    } else if (data_name == "Metabolite_lc") {
        u[,2] = -u[,2]
        v[,2] = -v[,2]
        u[,3] = -u[,3]
        v[,3] = -v[,3]
    } else if (data_name == "Metabolite_interpolated_gc") {
        u[,2] = -u[,2]
        v[,2] = -v[,2]
        u[,3] = -u[,3]
        v[,3] = -v[,3]
    } else if (data_name == "Metabolite_interpolated_lc") {
        # u[,2] = -u[,2]
        # v[,2] = -v[,2]
        u[,3] = -u[,3]
        v[,3] = -v[,3]
    }


    # save the matrix
    dir_svd <- file.path(dir_temp, "svd")
    if (!dir.exists(dir_svd)) dir.create(dir_svd, recursive = TRUE)
    file_s <- file.path(dir_svd, "S.txt")
    file_v <- file.path(dir_svd, "V.txt")
    file_u <- file.path(dir_svd, "U.txt")
    write.table(s, file = file_s, row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(v, file = file_v, row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(u, file = file_u, row.names = F, col.names = F, quote = F, sep = "\t")

    dir_loading_1 <- file.path(dir_result, data_name)
    if (!dir.exists(dir_loading_1)) dir.create(dir_loading_1, recursive = TRUE)
    dir_loading_2 <- file.path(dir_temp, "loadings")
    if (!dir.exists(dir_loading_2)) dir.create(dir_loading_2, recursive = TRUE)

    # save s value and ratio
    s_total <- sum(s[-1]^2)
    s_ratio <- s %>% 
        as_tibble() %>% 
        mutate(dim = row_number()-1,
               ratio = (value)^2/s_total) 
    write.xlsx(s_ratio, file.path(dir_loading_1, "s_values.xlsx"))
    write_csv(s_ratio, file.path(dir_loading_2, "s_values.csv"))

    # save loading data
    colnames(u) <- as.character(0:(number_dim-1))
    colnames(v) <- as.character(0:(number_dim-1))
    gene_loading <- data_gene %>%
        bind_cols(u)
    sample_loading <- data_sample %>% 
        bind_cols(v)
    file_sample_v <- file.path(dir_loading_1, "sample_loadings.xlsx")
    file_gene_u <- file.path(dir_loading_1, "gene_loadings.xlsx")
    write.xlsx(sample_loading, file_sample_v)
    write.xlsx(gene_loading, file_gene_u)

    # save loading data respectively
    for (i in 1:number_dim) {
        gene_loading <- data_gene %>%
            mutate(Loading = u[, i]) %>%
            arrange(desc(Loading)) 
        sample_loading <- data_sample %>%
            mutate(Loading = v[, i]) %>%
            arrange(desc(Loading))
        file_gene <- file.path(dir_loading_2, str_c("gene", i-1, "_loadings.txt"))
        file_sample <- file.path(dir_loading_2, str_c("sample", i-1, "_loadings.txt"))
        write_delim(gene_loading, file = file_gene, delim = "\t")
        write_delim(sample_loading, file = file_sample, delim = "\t")
    }
}



