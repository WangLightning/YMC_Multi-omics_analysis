suppressMessages({
    library(snowfall)
    library(zoo)
    library(rlist)
    library(tidyverse)
})


# parameters --------------------------------------------------------------

dir_result <- "../results"

# read config file
config <- list.load("config.yaml")

number_dim <- config$number_dim
CPUs <- config$CPUs
data_list <- config$data_histone_list

# normalization does not affect the result
is_norm <- "raw"
# is_norm <- "Scale"
# is_norm <- "Norm"
# interpolation_method <- "Linear"
# interpolation_method <- "Spline"
interpolation_method <- "Spline_natural"

dir_res <- file.path(dir_result, "ccf_histone_Spline_natural")
if (!dir.exists(dir_res)) dir.create(dir_res, recursive = TRUE)


# functions ---------------------------------------------------------------

# get time for RNA-seq and ChIP-seq data
getRNASeqTime <- function(x) {
    switch(
        x,
        T1 = 0.356,
        T2 = 0.666,
        T3 = 0.782,
        T4 = 0.866,
        T5 = 1.017,
        T6 = 1.187,
        T7 = 1.357,
        T8 = 1.516,
        T9 = 1.686,
        T10 = 1.856,
        T11 = 2.031,
        T12 = 2.196,
        T13 = 2.366,
        T14 = 2.908,
        T15 = 3.427,
        T16 = 3.987
    )
}

getChipSeqTime <- function(x) {
    switch(
        x,
        T1 = 0.351,
        T2 = 0.618,
        T3 = 0.753,
        T4 = 0.836,
        T5 = 0.985,
        T6 = 1.160,
        T7 = 1.373,
        T8 = 1.600,
        T9 = 1.770,
        T10 = 1.996,
        T11 = 2.220,
        T12 = 2.414,
        T13 = 2.741,
        T14 = 3.183,
        T15 = 3.547,
        T16 = 4.033
    )
}

scalize <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)  
}

normalize <- function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}


# read data ---------------------------------------------------------------

# rnaseq data
dir_rnaseq <- file.path(dir_result, "temp", "MUREN_log_YMC2014")
file_rnaseq_sample <- file.path(dir_rnaseq, "sample_name.txt")
file_rnaseq_v <- file.path(dir_rnaseq, "svd/V.txt")
data_rnaseq_v <- read_table(file_rnaseq_v,
                            col_names = as.character(0:(number_dim - 1)),
                            show_col_types = F)
data_ref <- read_table(file_rnaseq_sample,
                          col_names = "Sample",
                          show_col_types = F) %>%
    mutate(time = sapply(Sample, getRNASeqTime)) %>%
    bind_cols(data_rnaseq_v)

# CCF ---------------------------------------------------------------

results <- tibble(Data = NA, Dim = NA, Max_ACF = numeric(), Max_Lag = numeric()) 
for (data_name in data_list) {

    # read data
    dir_data <- file.path(dir_result, "temp", data_name)
    file_sample <- file.path(dir_data, "sample_name.txt")
    file_v <- file.path(dir_data, "svd/V.txt")
    data_v <- read_table(file_v,
                         col_names = as.character(0:(number_dim - 1)),
                         show_col_types = F)
    data <- read_table(file_sample,
                       col_names = "Sample",
                       show_col_types = F) %>%
        mutate(time = sapply(Sample, getChipSeqTime)) %>%
        bind_cols(data_v)

    # interpolation time points
    time_min <- max(min(data_ref$time), min(data$time))
    time_max <- min(max(data_ref$time), max(data$time))
    timepoints <- seq(from = time_min, to = time_max, by = 0.001)

    # different levels
    for(i in 1:2){

        if(i == 1){
            data_ref_1 <- data_ref %>% 
                select(time, value = `1`)
            data_1 <- data %>% 
                select(time, value = `1`)
        } else if(i == 2){
            data_ref_1 <- data_ref %>% 
                select(time, value = `2`)
            data_1 <- data %>% 
                select(time, value = `2`)
        }

        # interpolation
        if (interpolation_method == "Linear") {
            # linear
            interp_ref <- tibble(
                time = timepoints,
                value = coredata(na.approx(zoo(data_ref_1$value, order.by = data_ref_1$time), xout = timepoints))
            )
            interp <- tibble(
                time = timepoints,
                value = coredata(na.approx(zoo(data_1$value, order.by = data_1$time), xout = timepoints))
            )
        } else if (interpolation_method == "Spline") {
            # spline
            interp_ref <- tibble(
                time = timepoints,
                value = spline(data_ref_1$time, data_ref_1$value, xout = timepoints)$y 
            )
            interp <- tibble(
                time = timepoints,
                value = spline(data_1$time, data_1$value, xout = timepoints)$y
            )
        } else if (interpolation_method == "Spline_natural") {
            # natural spline
            interp_ref <- tibble(
                time = timepoints,
                value = spline(data_ref_1$time, data_ref_1$value, xout = timepoints, method = "natural")$y 
            )
            interp <- tibble(
                time = timepoints,
                value = spline(data_1$time, data_1$value, xout = timepoints, method = "natural")$y
            )
        }

        # normalize
        if (is_norm == "Norm") {
            interp_ref$value <- normalize(interp_ref$value)
            interp$value <- normalize(interp$value)
        }  else if (is_norm == "Scale") {
            interp_ref$value <- normalize(interp_ref$value)
            interp$value <- normalize(interp$value)
        }

        # cross-correlation
        ccf_result <- ccf(interp$value, interp_ref$value, lag.max = 1000, plot = F)
        lag_max <- ccf_result$lag[which.max(abs(ccf_result$acf))]
        corr_max <- ccf_result$acf[which.max(abs(ccf_result$acf))]

        results <- rbind(results, 
            data.frame(Data = data_name, Dim = i, Max_ACF = corr_max, Max_Lag = lag_max*0.001*60))

        # plot the results of cross-correlation
        pdf(file.path(dir_res, str_glue("CCF_{data_name}_{i}.pdf")), height = 5, width = 7)
        plot(ccf_result, main = str_glue("Optimal time difference is {lag_max*0.001*60} minutes, max ACF is {corr_max}"), xlab = "Lag (0.06 minutes)"  )
        abline(v = lag_max, col = "grey", lty = 2)
        dev.off()
    }

    # save results
    write_csv(results, file.path(dir_res, "CCF_histone_summary.csv"))

}

