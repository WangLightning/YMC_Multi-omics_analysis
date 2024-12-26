suppressMessages({
    library(rlist)
    library(tidyverse)
    library(snowfall)
})


# parameters --------------------------------------------------------------

# read config file
config <- list.load("config.yaml")

dir_result <- "../results"
number_dim <- config$number_dim
CPUs <- config$CPUs

file_annotation <- config$file_annotation
data_anno <- read_tsv(file_annotation, col_names = FALSE, show_col_types = FALSE) %>% 
    select(X2, X4, X5) %>% 
    rename(ID = X2, Symbol = X4, GeneTitle = X5)


# functions ---------------------------------------------------------------

permutation_under_constraint <- function(S) {
    m <- nrow(S)  
    n <- ncol(S) 
    
    # step 1: generate vector v
    v <- c(rep(0, n), rnorm(m - n))
    # v <- c(rep(0, n), runif(m-n, min = -1, max = 1))
    # step 2: calculate b = -S^T * v
    b <- -t(S) %*% v
    # step 3: get the first n rows of S and transpose it
    A <- t(S[1:n, ])
    # step 4: solve the linear equation Ax = b
    x <- solve(A, b)
    # step 5: replace the first n elements of v with the elements in x
    v[1:n] <- x
    # step 6: normalize vector v and return
    v <- v / sqrt(sum(v^2))
    
    return(v)
}

permutation_2 <- function(j){
    
    vec_permutation <- permutation_under_constraint(S)
    vec_res <- data_mat %*% vec_permutation
    vec_unit_res <- vec_res / sqrt(sum(vec_res^2))
    
    return(vec_unit_res)
}

calc_ratio <- function(row, u_value) {
    if (u_value > 0) {
        mean(row >= u_value)
    } else {
        mean(row <= u_value)
    }
}


# concatenate ---------------------------------------------------------------

data_name = "concatenate"

# read data
dir_temp <- file.path(dir_result, "temp", data_name)
file_gene <- file.path(dir_temp, "gene_name.txt")
file_mat <- file.path(dir_temp, "matrix.txt")
file_v <- file.path(dir_temp, "svd", "V.txt")
file_u <- file.path(dir_temp, "svd", "U.txt")
data_mat <- read_table(file_mat, col_names = F, show_col_types = FALSE) %>% as.matrix()
data_v <-  read_table(file_v, col_names = F, show_col_types = FALSE) %>% as.matrix()
data_u <-  read_table(file_u, col_names = F, show_col_types = FALSE) %>% as.matrix()
data_gene <- read_table(file_gene, col_names = "ID") %>% 
    left_join(data_anno, by = "ID") 

# dim 1 
i = 1
S <- matrix(rep(1, nrow(data_v))) 
vec_u <- data_u[,i]

sfInit(parallel = TRUE, cpus = CPUs)
sfExportAll()
mat_permutation <- sfSapply(1:100000, permutation_2)
sfStop()

mat_list <- split(mat_permutation, row(mat_permutation))
pvalues <- mapply(FUN = calc_ratio, mat_list, vec_u)

data_res <- data_gene %>% 
    bind_cols(vec_u) %>% 
    bind_cols(pvalues) %>% 
    set_names(c("ID", "Symbol", "GeneTitle", "Loading", "p value")) %>% 
    arrange(desc(Loading))

file_res <- file.path(dir_result, data_name, paste0("gene ", i, " pvalue ver2.csv"))
write_csv(data_res, file_res)


# dim 2 
i = 2
S <- matrix(rep(1, nrow(data_v))) %>% 
    cbind(data_v[,1:i-1])
vec_u <- data_u[,i]

sfInit(parallel = TRUE, cpus = CPUs)
sfExportAll()
mat_permutation <- sfSapply(1:100000, permutation_2)
sfStop()

mat_list <- split(mat_permutation, row(mat_permutation))
pvalues <- mapply(FUN = calc_ratio, mat_list, vec_u)

data_res <- data_gene %>% 
    bind_cols(vec_u) %>% 
    bind_cols(pvalues) %>% 
    set_names(c("ID", "Symbol", "GeneTitle", "Loading", "p value")) %>% 
    arrange(desc(Loading))

file_res <- file.path(dir_result, data_name, paste0("gene ", i, " pvalue ver2.csv"))
write_csv(data_res, file_res)

