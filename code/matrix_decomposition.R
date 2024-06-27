suppressMessages({
    library(rlist)
    library(openxlsx)
    library(tidyverse)
})


# parameters --------------------------------------------------------------

# read config file
config <- list.load("config.yaml")

dir_result <- "../results"
data_list <- config$data_list
number_dim <- config$number_dim


# main ---------------------------------------------------------------

# the annotation file
file_annotation <- "../data/annotation_SGD.tsv"
data_anno <- read_tsv(file_annotation, col_names = FALSE, show_col_types = FALSE) %>% 
    select(X2, X4, X5) %>% 
    rename(GeneID = X2, GeneSymbol = X4, GeneTitle = X5)

for (data_name in data_list) {

    if(data_name == "concatenate") next

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
        rename(GeneID = Symbol) %>% 
        left_join(data_anno, by = "GeneID") %>% 
        bind_cols(u)
    sample_loading <- data_sample %>% 
        bind_cols(v)
    file_sample_v <- file.path(dir_loading_1, "sample_loadings.xlsx")
    file_gene_u <- file.path(dir_loading_1, "gene_loadings.xlsx")
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
        file_gene <- file.path(dir_loading_2, str_c("gene", i-1, "_loadings.txt"))
        file_sample <- file.path(dir_loading_2, str_c("sample", i-1, "_loadings.txt"))
        write_delim(gene_loading, file = file_gene, delim = "\t")
        write_delim(sample_loading, file = file_sample, delim = "\t")
    }
}

# concatenate ---------------------------------------------------------------

data_name <- "concatenate"
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
for(i in 1:(number_dim-1)){
    if(data_1[,i+1] %*% data_2[,i] < 0){
        u[,i] = -u[,i]
        v[,i] = -v[,i]
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

dir_loading_1 <- file.path(dir_result, data_name)
if (!dir.exists(dir_loading_1)) dir.create(dir_loading_1, recursive = TRUE)
dir_loading_2 <- file.path(dir_temp, "loadings")
if (!dir.exists(dir_loading_2)) dir.create(dir_loading_2, recursive = TRUE)

# save s value and ratio
s_total <- sum(s^2)
s_ratio <- s %>% 
    as_tibble() %>% 
    mutate(dim = row_number(),
           ratio = (value)^2/s_total) 
write.xlsx(s_ratio, file.path(dir_loading_1, "s_values.xlsx"))
write_csv(s_ratio, file.path(dir_loading_2, "s_values.csv"))

# save loading data
colnames(u) <- as.character(1:number_dim)
colnames(v) <- as.character(1:number_dim)
gene_loading <- data_gene %>%
    rename(GeneID = Symbol) %>% 
    left_join(data_anno, by = "GeneID") %>% 
    bind_cols(u)
sample_loading <- data_sample %>% 
    bind_cols(v)
file_sample_v <- file.path(dir_loading_1, "sample_loadings.xlsx")
file_gene_u <- file.path(dir_loading_1, "gene_loadings.xlsx")
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
    file_gene <- file.path(dir_loading_2, str_c("gene", i, "_loadings.txt"))
    file_sample <- file.path(dir_loading_2, str_c("sample", i, "_loadings.txt"))
    write_delim(gene_loading, file = file_gene, delim = "\t")
    write_delim(sample_loading, file = file_sample, delim = "\t")
}

