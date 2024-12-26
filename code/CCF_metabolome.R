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
data_list <- config$data_metabolome_list

# normalization does not affect the result
is_norm <- "raw"
# is_norm <- "Scale"
# is_norm <- "Norm"
# interpolation_method <- "Linear"
# interpolation_method <- "Spline"
interpolation_method <- "Spline_natural"

dir_res <- file.path(dir_result, "ccf_metabolome_Spline_natural")
if (!dir.exists(dir_res)) dir.create(dir_res, recursive = TRUE)


# functions ---------------------------------------------------------------

# get time for microarray and metabolome data
getMicroarrayTime <- function(x) {
    switch(
        x,
        T1 = 0.417,
        T2 = 0.833,
        T3 = 1.250,
        T4 = 1.667,
        T5 = 2.083,
        T6 = 2.500,
        T7 = 2.917,
        T8 = 3.333,
        T9 = 3.750,
        T10 = 4.167,
        T11 = 4.583,
        T12 = 5.000,
        T13 = 5.417,
        T14 = 5.833,
        T15 = 6.250,
        T16 = 6.667,
        T17 = 7.083,
        T18 = 7.500,
        T19 = 7.917,
        T20 = 8.333,
        T21 = 8.750,
        T22 = 9.167,
        T23 = 9.583,
        T24 = 10.000,
        T25 = 10.417,
        T26 = 10.833,
        T27 = 11.250,
        T28 = 11.667,
        T29 = 12.083,
        T30 = 12.500,
        T31 = 12.917,
        T32 = 13.333,
        T33 = 13.750,
        T34 = 14.167,
        T35 = 14.583,
        T36 = 15.000
    )
}

getlcTime <- function(x) {
    switch(
        x,
        T1 = 0.474,
        T2 = 0.830,
        T3 = 1.342,
        T4 = 1.727,
        T5 = 2.111,
        T6 = 2.539,
        T7 = 3.051,
        T8 = 3.521,
        T9 = 3.892,
        T10 = 4.319,
        T11 = 4.746,
        T12 = 5.202,
        T13 = 5.358,
        T14 = 5.757,
        T15 = 6.142,
        T16 = 6.569,
        T17 = 7.010,
        T18 = 7.395,
        T19 = 7.737,
        T20 = 8.036,
        T21 = 8.506,
        T22 = 8.962,
        T23 = 9.432,
        T24 = 9.902
    )
}

getgcTime <- function(x) {
    switch(
        x,
        T1 = 0.488,
        T2 = 0.801,
        T3 = 1.186,
        T4 = 1.656,
        T5 = 2.168,
        T6 = 2.681,
        T7 = 3.108,
        T8 = 3.322,
        T9 = 3.835,
        T10 = 4.347,
        T11 = 4.817,
        T12 = 5.159,
        T13 = 5.358,
        T14 = 5.800,
        T15 = 6.227,
        T16 = 6.654,
        T17 = 7.096,
        T18 = 7.594,
        T19 = 8.079,
        T20 = 8.349,
        T21 = 8.862,
        T22 = 9.346,
        T23 = 9.802,
        T24 = 10.101
    )
}

scalize <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)  
}

normalize <- function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}


# read data ---------------------------------------------------------------

# reference: microarray data
dir_microarray <- file.path(dir_result, "temp", "microarray")
file_microarray_sample <- file.path(dir_microarray, "sample_name.txt")
file_microarray_v <- file.path(dir_microarray, "svd/V.txt")
data_microarray_v <- read_table(file_microarray_v,
                                col_names = as.character(0:(number_dim - 1)),
                                show_col_types = F)
data_ref <- read_table(file_microarray_sample,
                              col_names = "Sample",
                              show_col_types = F) %>%
    mutate(time = sapply(Sample, getMicroarrayTime)) %>%
    bind_cols(data_microarray_v)


# CCF ---------------------------------------------------------------------

results <- tibble(Data = NA, Dim = NA, Max_ACF = numeric(), Max_Lag = numeric()) 
for (data_name in data_list) {
    # read data
    dir_data <- file.path(dir_result, "temp", data_name)
    file_sample <- file.path(dir_data, "sample_name.txt")
    file_v <- file.path(dir_data, "svd/V.txt")
    data_v <- read_table(file_v,
                         col_names = as.character(0:(number_dim - 1)),
                         show_col_types = F)

    if (data_name == "Metabolite_lc" | data_name == "Metabolite_interpolated_lc") {
        data <- read_table(file_sample,
                           col_names = "Sample",
                           show_col_types = F) %>%
            mutate(time = sapply(Sample, getlcTime)) %>%
            bind_cols(data_v)
    } else if (data_name == "Metabolite_gc" | data_name == "Metabolite_interpolated_gc") {
        data <- read_table(file_sample,
                           col_names = "Sample",
                           show_col_types = F) %>%
            mutate(time = sapply(Sample, getgcTime)) %>%
            bind_cols(data_v)
    }
                

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

        # plot the results of interpolation
        p <- ggplot() +
            geom_line(data = interp_ref, aes(x = time, y = value, color = "Transcriptome of Tu et al. 2005")) +
            geom_line(data = interp, aes(x = time, y = value, color = "Metabolome")) +
            labs(title = str_glue("{interpolation_method} interpolation"), x = "Time (h)", y = "Loading") +
            theme_minimal()
        ggsave(file.path(dir_res, str_glue("Interpolation_{data_name}_{i}.pdf")), p, height = 5, width = 10)

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
    write_csv(results, file.path(dir_res, "CCF_metabolome_summary.csv"))

}



