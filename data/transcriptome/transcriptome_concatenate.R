suppressMessages({
    library(tidyverse)
})

# parameters --------------------------------------------------------------

dir_data <- "./"

# function -----------------------------------------------------------------

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


# main -----------------------------------------------------------------

file_data_m <- file.path(dir_data, "microarray.csv")
file_data_r1 <- file.path(dir_data, "MUREN_log_YMC2014.csv")

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

data_res <- data_all %>% 
    select(Symbol) %>%
    bind_cols(data_mat_all)
    set_names(c("Symbol", name_sample))

write_csv(data_res, file.path(dir_data, "transcriptome_concatenate.csv"))