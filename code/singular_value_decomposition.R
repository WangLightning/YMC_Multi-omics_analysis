suppressMessages({
    library(rlist)
    library(openxlsx)
    library(tidyverse)
})


# parameters --------------------------------------------------------------

dir_result <- "../results"

# read config file
config <- list.load("config.yaml")

number_dim <- config$number_dim
data_transcriptome_list <- config$data_transcriptome_list
data_histone_list <- config$data_histone_list
data_metabolome_list <- config$data_metabolome_list
# the gene annotation file
file_annotation <- config$file_annotation
data_anno <- read_tsv(file_annotation, col_names = FALSE, show_col_types = FALSE) %>% 
    select(X2, X4, X5) %>% 
    rename(GeneID = X2, GeneSymbol = X4, GeneTitle = X5)

dir_data <- config$dir_data

dir_res <- file.path(dir_result, "matrix_decomposition")
if (!dir.exists(dir_res)) dir.create(dir_res, recursive = TRUE)

# read data ---------------------------------------------------------------


for (data_name in c(data_transcriptome_list, data_histone_list, data_metabolome_list)) {

    # read data
    if (data_name %in% data_transcriptome_list) {
        file_data <- file.path(dir_data, "transcriptome", str_glue("{data_name}.csv"))
    } else if (data_name %in% data_histone_list) {
        file_data <- file.path(dir_data, "histone_modification", str_glue("{data_name}.csv"))
    } else if (data_name %in% data_metabolome_list) {
        file_data <- file.path(dir_data, "metabolome", str_glue("{data_name}.csv"))
    }

    data <- read_csv(file_data, show_col_types = FALSE)

    name_gene <- data$Symbol
    name_sample <- names(data)[-1]
    data_mat <- data %>%
        select(-1) %>%
        as.matrix()

    # result directory
    dir_temp <- file.path(dir_result, "temp", data_name)
    if (!dir.exists(dir_temp)) dir.create(dir_temp, recursive = TRUE)
    file_gene <- file.path(dir_temp, "gene_name.txt")
    file_sample <- file.path(dir_temp, "sample_name.txt")
    file_mat <- file.path(dir_temp, "matrix.txt")

    # save result
    write.table(name_gene, file = file_gene, row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(name_sample, file = file_sample, row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(data_mat, file = file_mat, row.names = F, col.names = F, sep = "\t")

}


# transcriptome & histone ---------------------------------------------------------------

data_list <- c(data_transcriptome_list, data_histone_list)

for (data_name in data_list) {

    if(data_name == "transcriptome_concatenate") next

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

    # correct sign
    data_ref <- data_list[1]
    if(data_name == data_ref){
        u[,1] = -u[,1]
        v[,1] = -v[,1]
        # u[,2] = -u[,2]
        # v[,2] = -v[,2]
        u[,3] = -u[,3]
        v[,3] = -v[,3]
        data_ref_u <- data_gene %>%
            bind_cols(u)
    }else{
        data_u <- data_gene %>%
            bind_cols(u)
        data_common_gene <- intersect(data_ref_u$Symbol, data_u$Symbol) %>%
            as_tibble() %>%
            set_names("Symbol")
        data_1 <- data_common_gene %>% 
            left_join(data_ref_u, by = "Symbol") %>%
            arrange(Symbol) %>% 
            distinct(Symbol, .keep_all = TRUE) %>%
            select(-Symbol) %>%
            as.matrix()
        data_2 <- data_common_gene %>% 
            left_join(data_u, by = "Symbol") %>%
            arrange(Symbol) %>% 
            distinct(Symbol, .keep_all = TRUE) %>%
            select(-Symbol) %>%
            as.matrix()

        for(i in 1:number_dim){
            if(data_1[,i] %*% data_2[,i] < 0){
                u[,i] = -u[,i]
                v[,i] = -v[,i]
            }
        }
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

    dir_loading <- file.path(dir_temp, "loadings")
    if (!dir.exists(dir_loading)) dir.create(dir_loading, recursive = TRUE)

    # save s value and ratio
    s_total <- sum(s[-1]^2)
    s_ratio <- s %>% 
        as_tibble() %>% 
        mutate(dim = row_number()-1,
               ratio = (value)^2/s_total,
               cum_ratio = ifelse(dim > 0, cumsum(ifelse(dim > 0, ratio, 0)), NA),
               distance = 1 - value / lag(value))
    write.xlsx(s_ratio, file.path(dir_res, str_glue("s_values_{data_name}.xlsx")))
    write_csv(s_ratio, file.path(dir_loading, "s_values.csv"))

    # save loading data
    colnames(u) <- as.character(0:(number_dim-1))
    colnames(v) <- as.character(0:(number_dim-1))
    gene_loading <- data_gene %>%
        rename(GeneID = Symbol) %>% 
        left_join(data_anno, by = "GeneID") %>% 
        bind_cols(u)
    sample_loading <- data_sample %>% 
        bind_cols(v)
    file_sample_v <- file.path(dir_res, str_glue("sample_loadings_{data_name}.xlsx"))
    file_gene_u <- file.path(dir_res, str_glue("gene_loadings_{data_name}.xlsx"))
    write.xlsx(sample_loading, file_sample_v)
    write.xlsx(gene_loading, file_gene_u)

    # save loading data respectively
    for (i in 1:number_dim) {
        # keep the maximum absolute value
        gene_loading <- data_gene %>%
            mutate(Loading = u[, i], a = abs(u[, i])) %>%
            arrange(desc(a)) %>%
            distinct(Symbol, .keep_all = T) %>%
            arrange(desc(Loading)) %>%
            select(Symbol, Loading)
        sample_loading <- data_sample %>%
            mutate(Loading = v[, i]) %>%
            arrange(desc(Loading))
        file_gene <- file.path(dir_loading, str_c("gene", i-1, "_loadings.txt"))
        file_sample <- file.path(dir_loading, str_c("sample", i-1, "_loadings.txt"))
        write_delim(gene_loading, file = file_gene, delim = "\t")
        write_delim(sample_loading, file = file_sample, delim = "\t")
    }
}

# transcriptome_concatenate ---------------------------------------------------------------

data_name <- "transcriptome_concatenate"
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

# correct sign
u[,1] = -u[,1]
v[,1] = -v[,1]
u[,2] = -u[,2]
v[,2] = -v[,2]

# save the matrix
dir_svd <- file.path(dir_temp, "svd")
if (!dir.exists(dir_svd)) dir.create(dir_svd, recursive = TRUE)
file_s <- file.path(dir_svd, "S.txt")
file_v <- file.path(dir_svd, "V.txt")
file_u <- file.path(dir_svd, "U.txt")
write.table(s, file = file_s, row.names = F, col.names = F, quote = F, sep = "\t")
write.table(v, file = file_v, row.names = F, col.names = F, quote = F, sep = "\t")
write.table(u, file = file_u, row.names = F, col.names = F, quote = F, sep = "\t")

dir_loading <- file.path(dir_temp, "loadings")
if (!dir.exists(dir_loading)) dir.create(dir_loading, recursive = TRUE)

# save s value and ratio
s_total <- sum(s^2)
s_ratio <- s %>% 
    as_tibble() %>% 
    mutate(dim = row_number(),
           ratio = (value)^2/s_total,
           cum_ratio = cumsum(ratio),
           distance = 1 - value / lag(value)) 
write.xlsx(s_ratio, file.path(dir_res, str_glue("s_values_{data_name}.xlsx")))
write_csv(s_ratio, file.path(dir_loading, "s_values.csv"))

# save loading data
colnames(u) <- as.character(1:number_dim)
colnames(v) <- as.character(1:number_dim)
gene_loading <- data_gene %>%
    rename(GeneID = Symbol) %>% 
    left_join(data_anno, by = "GeneID") %>% 
    bind_cols(u)
sample_loading <- data_sample %>% 
    bind_cols(v)
file_sample_v <- file.path(dir_res, str_glue("sample_loadings_{data_name}.xlsx"))
file_gene_u <- file.path(dir_res, str_glue("gene_loadings_{data_name}.xlsx"))
write.xlsx(sample_loading, file_sample_v)
write.xlsx(gene_loading, file_gene_u)

# save loading data respectively
for (i in 1:number_dim) {
    # keep the maximum absolute value
    gene_loading <- data_gene %>%
        mutate(Loading = u[, i], a = abs(u[, i])) %>%
        arrange(desc(a)) %>%
        distinct(Symbol, .keep_all = T) %>%
        arrange(desc(Loading)) %>%
        select(Symbol, Loading)
    sample_loading <- data_sample %>%
        mutate(Loading = v[, i]) %>%
        arrange(desc(Loading))
    file_gene <- file.path(dir_loading, str_c("gene", i, "_loadings.txt"))
    file_sample <- file.path(dir_loading, str_c("sample", i, "_loadings.txt"))
    write_delim(gene_loading, file = file_gene, delim = "\t")
    write_delim(sample_loading, file = file_sample, delim = "\t")
}

# metabolome ---------------------------------------------------------------

for (data_name in data_metabolome_list) {

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

    dir_loading <- file.path(dir_temp, "loadings")
    if (!dir.exists(dir_loading)) dir.create(dir_loading, recursive = TRUE)

    # save s value and ratio
    s_total <- sum(s[-1]^2)
    s_ratio <- s %>% 
        as_tibble() %>% 
        mutate(dim = row_number()-1,
               ratio = (value)^2/s_total,
               cum_ratio = ifelse(dim > 0, cumsum(ifelse(dim > 0, ratio, 0)), NA),
               distance = 1 - value / lag(value))
    write.xlsx(s_ratio, file.path(dir_res, str_glue("s_values_{data_name}.xlsx")))
    write_csv(s_ratio, file.path(dir_loading, "s_values.csv"))

    # save loading data
    colnames(u) <- as.character(0:(number_dim-1))
    colnames(v) <- as.character(0:(number_dim-1))
    gene_loading <- data_gene %>%
        bind_cols(u)
    sample_loading <- data_sample %>% 
        bind_cols(v)
    file_sample_v <- file.path(dir_res, str_glue("sample_loadings_{data_name}.xlsx"))
    file_gene_u <- file.path(dir_res, str_glue("gene_loadings_{data_name}.xlsx"))
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
        file_gene <- file.path(dir_loading, str_c("gene", i-1, "_loadings.txt"))
        file_sample <- file.path(dir_loading, str_c("sample", i-1, "_loadings.txt"))
        write_delim(gene_loading, file = file_gene, delim = "\t")
        write_delim(sample_loading, file = file_sample, delim = "\t")
    }
}
