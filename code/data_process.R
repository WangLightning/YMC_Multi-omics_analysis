suppressMessages({
        library(rlist)
        library(tidyverse)
})


# parameters --------------------------------------------------------------

# read config file
config <- list.load("config.yaml")

dir_result <- "../results"
data_list <- config$data_list

# functions ---------------------------------------------------------------

RemoveBaseline <- function(m){
    res <- svd(m)
    s1 <- res$d[1]
    u1 <- res$u[,1]
    v1 <- res$v[,1]
    m2 <- m - u1 %*% t(v1) * s1
    return(m2)
}

Fnorm_sacle <- function(m, coef = 1){
    # rescale by F-norm， the coefficient is p 
    m2 <- (m * coef)/norm(m, "F")
    return(m2)
}


# process concatenated data ---------------------------------------------------------------

file_data_m <- "../data/microarray.csv"
file_data_r1 <- "../data/rnaseq14.csv"

name_sample_m <- paste0("MT",1:36)
name_sample_r1 <- paste0("RT",1:16)

data_m <- read_csv(file_data_m, show_col_types = FALSE) %>% 
    filter(str_detect(Symbol, "^Y")) %>% 
    set_names(c('Symbol', name_sample_m))
data_r1 <- read_csv(file_data_r1, show_col_types = FALSE) %>% 
    filter(str_detect(Symbol, "^Y")) %>% 
    set_names(c('Symbol', name_sample_r1))
data_all <- data_m %>% 
    inner_join(data_r1, by = "Symbol")

name_gene <- data_all$Symbol
name_sample <- names(data_all)[-1]

data_mat_m <- data_all %>%
    select(all_of(name_sample_m)) %>%
    as.matrix() %>%
    RemoveBaseline() %>%
    Fnorm_sacle(6.4)
data_mat_r1 <- data_all %>% 
    select(all_of(name_sample_r1)) %>%
    as.matrix() %>% 
    RemoveBaseline() %>%
    Fnorm_sacle(6)

data_mat_all <- data_mat_m %>% 
    bind_cols(data_mat_r1) 

# result directory
dir_temp <- file.path(dir_result, "temp", "concatenate")
if (!dir.exists(dir_temp)) dir.create(dir_temp, recursive = TRUE)
file_gene <- file.path(dir_temp, "gene_name.txt")
file_sample <- file.path(dir_temp, "sample_name.txt")
file_mat <- file.path(dir_temp, "matrix.txt")

# save result
write.table(name_gene, file = file_gene, row.names = F, col.names = F, quote = F, sep = "\t")
write.table(name_sample, file = file_sample, row.names = F, col.names = F, quote = F, sep = "\t")
write.table(data_mat_all, file = file_mat, row.names = F, col.names = F, sep = "\t")


# process other data ---------------------------------------------------------------

for (data_name in data_list) {

    if(data_name == "concatenate") next

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

}