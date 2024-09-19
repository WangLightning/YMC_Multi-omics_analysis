suppressMessages({
    library(snowfall)
    library(zoo)
    library(svglite)
    library(rlist)
    library(tidyverse)
})


# parameters --------------------------------------------------------------

# read config file
config <- list.load("../../code/config.yaml")
number_dim <- config$number_dim
CPUs <- config$CPUs

# time_align <- "dtw"
time_align <- "ddtw"
# time_align <- "LinearTransform"
# is_norm <- "raw"
# is_norm <- "Scale"
is_norm <- "Norm"
interpolation_method <- "Linear"
# interpolation_method <- "Spline"
time_range <- "all"
# time_range <- "sub1"
# time_range <- "sub2"
# time_range <- "sub3"
align_method <- "ccf"
# align_method <- "distance1"
# align_method <- "distance2"

dir_res <- str_glue("./metabolite_lc_{time_align}_{is_norm}_{interpolation_method}_{align_method}_{time_range}")
if (!dir.exists(dir_res)) dir.create(dir_res, recursive = TRUE)


# functions ---------------------------------------------------------------

if (time_align == "ddtw") {
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

} 

scalize <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)  
}
normalize <- function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}


distance_ts <- function(ts1, ts2, lag = 1000) {
    
    calculate_distance <- function(i, ts1, ts2) {
        shifted_ts1 <- if (i < 0) {
            c(rep(NA, abs(i)), ts1[1:(length(ts1) + i)])
        } else {
            c(ts1[(i + 1):length(ts1)], rep(NA, i))
        }

        mean(abs(shifted_ts1 - ts2), na.rm = TRUE)
    }

    l <- min(lag, length(ts1))
    
    distances <- mapply(calculate_distance, -l:l, MoreArgs = list(ts1 = ts1, ts2 = ts2))
    
    distance_df <- data.frame(Time_Difference = -l:l, Distance = distances)
    
    return(distance_df)
}

distance2_ts <- function(ts1, ts2, lag = 1000) {
    
    calculate_distance <- function(i, ts1, ts2) {
        shifted_ts1 <- if (i < 0) {
            c(rep(NA, abs(i)), ts1[1:(length(ts1) + i)])
        } else {
            c(ts1[(i + 1):length(ts1)], rep(NA, i))
        }

        specific_times <- c(
            0.417,
            0.833,
            1.250,
            1.667,
            2.083,
            2.500,
            2.917,
            3.333,
            3.750,
            4.167,
            4.583,
            5.000,
            5.417,
            5.833,
            6.250,
            6.667,
            7.083,
            7.500,
            7.917,
            8.333,
            8.750,
            9.167,
            9.583,
            10.000,
            10.417,
            10.833,
            11.250,
            11.667,
            12.083,
            12.500,
            12.917,
            13.333,
            13.750,
            14.167,
            14.583,
            15.000
        )  
        
        selected_ts2 <- rep(NA, length(ts2))  
        indices <- match(specific_times, timepoints)
        valid_indices <- !is.na(indices)
        selected_ts2[indices[valid_indices]] <- ts2[indices[valid_indices]]
        
        mean(abs(shifted_ts1 - selected_ts2), na.rm = TRUE)
    }
    
    l <- min(lag, length(ts1))
    
    distances <- mapply(calculate_distance, -l:l, MoreArgs = list(ts1 = ts1, ts2 = ts2))
    
    distance_df <- data.frame(Time_Difference = -l:l, Distance = distances)
    
    return(distance_df)
}


# read data ---------------------------------------------------------------

# microarray data
file_microarray_sample <-
    "../../results/temp/microarray/sample_name.txt"
file_microarray_v <- "../../results/temp/microarray/svd/V.txt"
data_microarray_sample <-
    read_table(file_microarray_sample,
               col_names = "Sample",
               show_col_types = F) %>%
    mutate(Time = sapply(Sample, getMicroarrayTime))
data_microarray_v <- read_table(file_microarray_v,
                                col_names = as.character(0:(number_dim - 1)),
                                show_col_types = F)
data_microarray <- data_microarray_sample %>%
    bind_cols(data_microarray_v)

# # Metabolite_gc
# file_metabolite_sample <-
#     "../../results/temp/Metabolite_gc/sample_name.txt"
# file_metabolite_v <- "../../results/temp/Metabolite_gc/svd/V.txt"
# data_metabolite_sample <-
#     read_table(file_metabolite_sample,
#                col_names = "Sample",
#                show_col_types = F) %>%
#     mutate(Time = sapply(Sample, getgcTime))
# data_metabolite_v <-  read_table(file_metabolite_v,
#                                  col_names = as.character(0:(number_dim - 1)),
#                                  show_col_types = F)
# data_metabolite <- data_metabolite_sample %>%
#     bind_cols(data_metabolite_v)

# # Metabolite_lc
# file_metabolite_sample <-
#     "../../results/temp/Metabolite_lc/sample_name.txt"
# file_metabolite_v <- "../../results/temp/Metabolite_lc/svd/V.txt"
# data_metabolite_sample <-
#     read_table(file_metabolite_sample,
#                col_names = "Sample",
#                show_col_types = F) %>%
#     mutate(Time = sapply(Sample, getlcTime))
# data_metabolite_v <-  read_table(file_metabolite_v,
#                                  col_names = as.character(0:(number_dim - 1)),
#                                  show_col_types = F)
# data_metabolite <- data_metabolite_sample %>%
#     bind_cols(data_metabolite_v)

# # Metabolite_interpolated_gc
# file_metabolite_sample <-
#     "../../results/temp/Metabolite_interpolated_gc/sample_name.txt"
# file_metabolite_v <-
#     "../../results/temp/Metabolite_interpolated_gc/svd/V.txt"
# data_metabolite_sample <-
#     read_table(file_metabolite_sample,
#                col_names = "Sample",
#                show_col_types = F) %>%
#     mutate(Time = sapply(Sample, getgcTime))
# data_metabolite_v <-  read_table(file_metabolite_v,
#                                  col_names = as.character(0:(number_dim - 1)),
#                                  show_col_types = F)
# data_metabolite <- data_metabolite_sample %>%
#     bind_cols(data_metabolite_v)

# Metabolite_interpolated_lc
file_metabolite_sample <-
    "../../results/temp/Metabolite_interpolated_lc/sample_name.txt"
file_metabolite_v <-
    "../../results/temp/Metabolite_interpolated_lc/svd/V.txt"
data_metabolite_sample <-
    read_table(file_metabolite_sample,
               col_names = "Sample",
               show_col_types = F) %>%
    mutate(Time = sapply(Sample, getlcTime))
data_metabolite_v <-  read_table(file_metabolite_v,
                                 col_names = as.character(0:(number_dim - 1)),
                                 show_col_types = F)
data_metabolite <- data_metabolite_sample %>%
    bind_cols(data_metabolite_v)


# 插值时间点
if (time_range == "all") {
    time_min <- max(min(data_microarray$Time), min(data_metabolite$Time))
    time_max <- min(max(data_microarray$Time), max(data_metabolite$Time))
} else if (time_range == "sub1") {
    time_min <- 0.8
    time_max <- 9.9
}

timepoints <- seq(from = time_min, to = time_max, by = 0.001)


# dim 1 -------------------------------------------------------------------

data_microarray_1 <- data_microarray %>% 
    select(time = Time, value = `1`)
data_metabolite_1 <- data_metabolite %>% 
    select(time = Time, value = `1`)

# 插值计算
if (interpolation_method == "Linear") {
    # linear
    interp_microarray <- tibble(
        time = timepoints,
        value = coredata(na.approx(zoo(data_microarray_1$value, order.by = data_microarray_1$time), xout = timepoints))
    )
    interp_metabolite <- tibble(
        time = timepoints,
        value = coredata(na.approx(zoo(data_metabolite_1$value, order.by = data_metabolite_1$time), xout = timepoints))
    )
} else if (interpolation_method == "Spline") {
    # spline
    interp_microarray <- tibble(
        time = timepoints,
        value = spline(data_microarray_1$time, data_microarray_1$value, xout = timepoints)$y 
    )
    interp_metabolite <- tibble( 
        time = timepoints,
        value = spline(data_metabolite_1$time, data_metabolite_1$value, xout = timepoints)$y
    )
}


if (is_norm == "Norm") {
    interp_microarray$value <- normalize(interp_microarray$value)
    interp_metabolite$value <- normalize(interp_metabolite$value)
} else if (is_norm == "Scale") {
    interp_microarray$value <- scalize(interp_microarray$value)
    interp_metabolite$value <- scalize(interp_metabolite$value)
}

p <- ggplot() +
    geom_line(data = interp_microarray, aes(x = time, y = value, color = "Transcriptome of Tu et al. 2005")) +
    geom_line(data = interp_metabolite, aes(x = time, y = value, color = "Metabolome")) +
    labs(title = str_glue("{interpolation_method} interpolation"), x = "Time (h)", y = "Loading") +
    theme_minimal()
ggsave(file.path(dir_res, "Interpolation_dim1.pdf"), p, height = 5, width = 10)

if (align_method == "ccf") {
    
    
    ccf_microarray_metabolite <- ccf(interp_metabolite$value, interp_microarray$value, lag.max = 1500, plot = F)
    timing_microarray_metabolite <- ccf_microarray_metabolite$lag[which.max(ccf_microarray_metabolite$acf)]
    # plot
    # pdf(file.path(dir_res, "CCF_dim1.pdf"), height = 7, width = 10)
    svg(file.path(dir_res, "CCF_dim1.svg"), height = 5, width = 7)
    plot(ccf_microarray_metabolite, main = str_glue("The optimal time difference is {timing_microarray_metabolite*0.001*60} minutes"), xlab = "Lag (0.06 minutes)"  )
    abline(v = timing_microarray_metabolite, col = "grey", lty = 2)
    dev.off()
    
} else if (align_method == "distance1") {
    
    
    dis_microarray_metabolite <- distance_ts(interp_metabolite$value, interp_microarray$value)
    timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]
    
    p <- ggplot(dis_microarray_metabolite, aes(x = Time_Difference*0.001, y = Distance)) +
        geom_line() +
        geom_vline(xintercept = timing_microarray_metabolite*0.001, color = "grey", lty = 2) +
        labs(x = "Time difference (h)", y = "Distance", title = str_glue("The optimal time difference is {timing_microarray_metabolite*0.001*60} minutes")) +
        theme_classic()
    # ggsave(file.path(dir_res, "Distance_dim1.pdf"), p, height = 7, width = 10)
    ggsave(file.path(dir_res, "Distance_dim1.svg"), p, height = 5, width = 7, device = "svg")
    
} else if (align_method == "distance2") {
    
    
    dis_microarray_metabolite <- distance2_ts(interp_metabolite$value, interp_microarray$value)
    timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]
    
    p <- ggplot(dis_microarray_metabolite, aes(x = Time_Difference*0.001, y = Distance)) +
        geom_line() +
        geom_vline(xintercept = timing_microarray_metabolite*0.001, color = "grey", lty = 2) +
        labs(x = "Time difference (h)", y = "Distance", title = str_glue("The optimal time difference is {timing_microarray_metabolite*0.001*60} minutes")) +
        theme_classic()
    # ggsave(file.path(dir_res, "Distance_dim1.pdf"), p, height = 7, width = 10)
    ggsave(file.path(dir_res, "Distance_dim1.svg"), p, height = 5, width = 7, device = "svg")
    
}


ts1 <- interp_metabolite$value
ts2 <- interp_microarray$value
i <- timing_microarray_metabolite
shifted_ts1 <- if (i < 0) {
    c(rep(NA, abs(i)), ts1[1:(length(ts1) + i)])
} else {
    c(ts1[(i + 1):length(ts1)], rep(NA, i))
}
p <- ggplot() +
    geom_line(aes(x = 1:(length(ts1))*0.001, y = shifted_ts1, color = "Metabolome")) +
    geom_line(aes(x = 1:(length(ts1))*0.001, y = ts2, color = "Transcriptome of Tu et al. 2005")) +
    labs(x = "Time (h)", y = "Loading", title = "The Optimal alignment of curves") +
    theme_minimal()
ggsave(file.path(dir_res, "Optimal_alignment_dim1.pdf"), p, height = 5, width = 10)


# dim 2 -------------------------------------------------------------------

data_microarray_1 <- data_microarray %>% 
    select(time = Time, value = `2`)
data_metabolite_1 <- data_metabolite %>% 
    select(time = Time, value = `2`)


if (interpolation_method == "Linear") {
    # linear
    interp_microarray <- tibble(
        time = timepoints,
        value = coredata(na.approx(zoo(data_microarray_1$value, order.by = data_microarray_1$time), xout = timepoints))
    )
    interp_metabolite <- tibble(
        time = timepoints,
        value = coredata(na.approx(zoo(data_metabolite_1$value, order.by = data_metabolite_1$time), xout = timepoints))
    )
} else if (interpolation_method == "Spline") {
    # spline
    interp_microarray <- tibble(
        time = timepoints,
        value = spline(data_microarray_1$time, data_microarray_1$value, xout = timepoints)$y 
    )
    interp_metabolite <- tibble( 
        time = timepoints,
        value = spline(data_metabolite_1$time, data_metabolite_1$value, xout = timepoints)$y
    )
}


if (is_norm == "Norm") {
    interp_microarray$value <- normalize(interp_microarray$value)
    interp_metabolite$value <- normalize(interp_metabolite$value)
} else if (is_norm == "Scale") {
    interp_microarray$value <- scalize(interp_microarray$value)
    interp_metabolite$value <- scalize(interp_metabolite$value)
}


p <- ggplot() +
    geom_line(data = interp_microarray, aes(x = time, y = value, color = "Transcriptome of Tu et al. 2005")) +
    geom_line(data = interp_metabolite, aes(x = time, y = value, color = "Metabolome")) +
    labs(title = str_glue("{interpolation_method} interpolation"), x = "Time (h)", y = "Loading") +
    theme_minimal()
ggsave(file.path(dir_res, "Interpolation_dim2.pdf"), p, height = 5, width = 10)


if (align_method == "ccf") {
    
    
    ccf_microarray_metabolite <- ccf(interp_metabolite$value, interp_microarray$value, lag.max = 1000, plot = F)
    timing_microarray_metabolite <- ccf_microarray_metabolite$lag[which.max(ccf_microarray_metabolite$acf)]
    # plot
    # pdf(file.path(dir_res, "CCF_dim2.pdf"), height = 7, width = 10)
    svg(file.path(dir_res, "CCF_dim2.svg"), height = 5, width = 7)
    plot(ccf_microarray_metabolite, main = str_glue("Cross correlation between curves; the timing is {timing_microarray_metabolite*0.001*60} minutes"), xlab = "Lag (0.06 minutes)"  )
    abline(v = timing_microarray_metabolite, col = "grey", lty = 2)
    dev.off()
    
} else if (align_method == "distance1") {
    
    
    dis_microarray_metabolite <- distance_ts(interp_metabolite$value, interp_microarray$value)
    timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]
    
    p <- ggplot(dis_microarray_metabolite, aes(x = Time_Difference*0.001, y = Distance)) +
        geom_line() +
        geom_vline(xintercept = timing_microarray_metabolite*0.001, color = "grey", lty = 2) +
        labs(x = "Time difference (h)", y = "Distance", title = str_glue("The optimal time difference is {timing_microarray_metabolite*0.001*60} minutes")) +
        theme_classic()
    # ggsave(file.path(dir_res, "Distance_dim2.pdf"), p, height = 7, width = 10)
    ggsave(file.path(dir_res, "Distance_dim2.svg"), p, height = 5, width = 7, device = "svg")
    
} else if (align_method == "distance2") {
    
    
    dis_microarray_metabolite <- distance2_ts(interp_metabolite$value, interp_microarray$value)
    timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]
    
    p <- ggplot(dis_microarray_metabolite, aes(x = Time_Difference*0.001, y = Distance)) +
        geom_line() +
        geom_vline(xintercept = timing_microarray_metabolite*0.001, color = "grey", lty = 2) +
        labs(x = "Time Difference (h)", y = "Distance", title = str_glue("Distance between curves; the timing is {timing_microarray_metabolite*0.001*60} minutes")) +
        theme_classic()
    ggsave(file.path(dir_res, "Distance_dim2.pdf"), p, height = 7, width = 10)
    
}


ts1 <- interp_metabolite$value
ts2 <- interp_microarray$value
i <- timing_microarray_metabolite
shifted_ts1 <- if (i < 0) {
    c(rep(NA, abs(i)), ts1[1:(length(ts1) + i)])
} else {
    c(ts1[(i + 1):length(ts1)], rep(NA, i))
}
p <- ggplot() +
    geom_line(aes(x = 1:(length(ts1))*0.001, y = shifted_ts1, color = "Metabolome")) +
    geom_line(aes(x = 1:(length(ts1))*0.001, y = ts2, color = "Transcriptome of Tu et al. 2005")) +
    labs(x = "Time (h)", y = "Loading", title = "The Optimal alignment of curves") +
    theme_minimal()
ggsave(file.path(dir_res, "Optimal_alignment_dim2.pdf"), p, height = 5, width = 10)



# next phase dim 1-2 -------------------------------------------------------------------

data_microarray_1 <- data_microarray %>% 
    select(time = Time, value = `1`)
data_metabolite_1 <- data_metabolite %>% 
    select(time = Time, value = `2`) 

# 插值计算
if (interpolation_method == "Linear") {
    # linear
    interp_microarray <- tibble(
        time = timepoints,
        value = coredata(na.approx(zoo(data_microarray_1$value, order.by = data_microarray_1$time), xout = timepoints))
    )
    interp_metabolite <- tibble(
        time = timepoints,
        value = coredata(na.approx(zoo(data_metabolite_1$value, order.by = data_metabolite_1$time), xout = timepoints))
    )
} else if (interpolation_method == "Spline") {
    # spline
    interp_microarray <- tibble(
        time = timepoints,
        value = spline(data_microarray_1$time, data_microarray_1$value, xout = timepoints)$y 
    )
    interp_metabolite <- tibble( 
        time = timepoints,
        value = spline(data_metabolite_1$time, data_metabolite_1$value, xout = timepoints)$y
    )
}

# 数据标准化
if (is_norm == "Norm") {
    interp_microarray$value <- normalize(interp_microarray$value)
    interp_metabolite$value <- normalize(interp_metabolite$value)
} else if (is_norm == "Scale") {
    interp_microarray$value <- scalize(interp_microarray$value)
    interp_metabolite$value <- scalize(interp_metabolite$value)
}

# 画出插值结果
p <- ggplot() +
    geom_line(data = interp_microarray, aes(x = time, y = value, color = "Transcriptome of Tu et al. 2005")) +
    geom_line(data = interp_metabolite, aes(x = time, y = value, color = "Metabolome")) +
    labs(title = str_glue("{interpolation_method} interpolation"), x = "Time (h)", y = "Loading") +
    theme_minimal()
ggsave(file.path(dir_res, "Interpolation_dim1-2.pdf"), p, height = 5, width = 10)

# 计算最优对齐
if (align_method == "ccf") {
    
    # 计算ccf，并找到相位差
    ccf_microarray_metabolite <- ccf(interp_metabolite$value, interp_microarray$value, lag.max = 1000, plot = F)
    timing_microarray_metabolite <- ccf_microarray_metabolite$lag[which.max(ccf_microarray_metabolite$acf)]
    # plot
    # pdf(file.path(dir_res, "CCF_dim1-2.pdf"), height = 7, width = 10)
    svg(file.path(dir_res, "CCF_dim1-2.svg"), height = 5, width = 7)
    plot(ccf_microarray_metabolite, main = str_glue("The optimal time difference is {timing_microarray_metabolite*0.001*60} minutes"), xlab = "Lag (0.06 minutes)"  )
    abline(v = timing_microarray_metabolite, col = "grey", lty = 2)
    dev.off()
    
} else if (align_method == "distance1") {
    
    # 计算距离，并找到时间差
    dis_microarray_metabolite <- distance_ts(interp_metabolite$value, interp_microarray$value, lag = 3000)
    timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]
    # 绘制距离随时间差变化的图
    p <- ggplot(dis_microarray_metabolite, aes(x = Time_Difference*0.001, y = Distance)) +
        geom_line() +
        geom_vline(xintercept = timing_microarray_metabolite*0.001, color = "grey", lty = 2) +
        labs(x = "Time difference (h)", y = "Distance", title = str_glue("The optimal time difference is {timing_microarray_metabolite*0.001*60} minutes")) +
        theme_classic()
    # ggsave(file.path(dir_res, "Distance_dim1.pdf"), p, height = 7, width = 10)
    ggsave(file.path(dir_res, "Distance_dim1-2.svg"), p, height = 5, width = 7, device = "svg")
    
} else if (align_method == "distance2") {
    
    # 计算距离，并找到时间差
    dis_microarray_metabolite <- distance2_ts(interp_metabolite$value, interp_microarray$value)
    timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]
    # 绘制距离随时间差变化的图
    p <- ggplot(dis_microarray_metabolite, aes(x = Time_Difference*0.001, y = Distance)) +
        geom_line() +
        geom_vline(xintercept = timing_microarray_metabolite*0.001, color = "grey", lty = 2) +
        labs(x = "Time difference (h)", y = "Distance", title = str_glue("The optimal time difference is {timing_microarray_metabolite*0.001*60} minutes")) +
        theme_classic()
    ggsave(file.path(dir_res, "Distance_dim1-2.pdf"), p, height = 7, width = 10)
    
}

# 绘制最优对齐曲线
ts1 <- interp_metabolite$value
ts2 <- interp_microarray$value
i <- timing_microarray_metabolite
shifted_ts1 <- if (i < 0) {
    c(rep(NA, abs(i)), ts1[1:(length(ts1) + i)])
} else {
    c(ts1[(i + 1):length(ts1)], rep(NA, i))
}
p <- ggplot() +
    geom_line(aes(x = 1:(length(ts1))*0.001, y = shifted_ts1, color = "Metabolome")) +
    geom_line(aes(x = 1:(length(ts1))*0.001, y = ts2, color = "Transcriptome of Tu et al. 2005")) +
    labs(x = "Time (h)", y = "Loading", title = "The Optimal alignment of curves") +
    theme_minimal()
ggsave(file.path(dir_res, "Optimal_alignment_dim1-2.pdf"), p, height = 5, width = 10)



# next phase dim 2-1 -------------------------------------------------------------------

data_microarray_1 <- data_microarray %>% 
    select(time = Time, value = `2`)
data_metabolite_1 <- data_metabolite %>% 
    select(time = Time, value = `1`) %>% 
    mutate(value = -value)

# 插值计算
if (interpolation_method == "Linear") {
    # linear
    interp_microarray <- tibble(
        time = timepoints,
        value = coredata(na.approx(zoo(data_microarray_1$value, order.by = data_microarray_1$time), xout = timepoints))
    )
    interp_metabolite <- tibble(
        time = timepoints,
        value = coredata(na.approx(zoo(data_metabolite_1$value, order.by = data_metabolite_1$time), xout = timepoints))
    )
} else if (interpolation_method == "Spline") {
    # spline
    interp_microarray <- tibble(
        time = timepoints,
        value = spline(data_microarray_1$time, data_microarray_1$value, xout = timepoints)$y 
    )
    interp_metabolite <- tibble( 
        time = timepoints,
        value = spline(data_metabolite_1$time, data_metabolite_1$value, xout = timepoints)$y
    )
}

# 数据标准化
if (is_norm == "Norm") {
    interp_microarray$value <- normalize(interp_microarray$value)
    interp_metabolite$value <- normalize(interp_metabolite$value)
} else if (is_norm == "Scale") {
    interp_microarray$value <- scalize(interp_microarray$value)
    interp_metabolite$value <- scalize(interp_metabolite$value)
}

# 画出插值结果
p <- ggplot() +
    geom_line(data = interp_microarray, aes(x = time, y = value, color = "Transcriptome of Tu et al. 2005")) +
    geom_line(data = interp_metabolite, aes(x = time, y = value, color = "Metabolome")) +
    labs(title = str_glue("{interpolation_method} interpolation"), x = "Time (h)", y = "Loading") +
    theme_minimal()
ggsave(file.path(dir_res, "Interpolation_dim2-1.pdf"), p, height = 5, width = 10)

# 计算最优对齐
if (align_method == "ccf") {
    
    # 计算ccf，并找到相位差
    ccf_microarray_metabolite <- ccf(interp_metabolite$value, interp_microarray$value, lag.max = 1000, plot = F)
    timing_microarray_metabolite <- ccf_microarray_metabolite$lag[which.max(ccf_microarray_metabolite$acf)]
    # plot
    # pdf(file.path(dir_res, "CCF_dim2-1.pdf"), height = 7, width = 10)
    svg(file.path(dir_res, "CCF_dim2-1.svg"), height = 5, width = 7)
    plot(ccf_microarray_metabolite, main = str_glue("The optimal time difference is {timing_microarray_metabolite*0.001*60} minutes"), xlab = "Lag (0.06 minutes)"  )
    abline(v = timing_microarray_metabolite, col = "grey", lty = 2)
    dev.off()
    
} else if (align_method == "distance1") {
    
    # 计算距离，并找到时间差
    dis_microarray_metabolite <- distance_ts(interp_metabolite$value, interp_microarray$value)
    timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]
    # 绘制距离随时间差变化的图
    p <- ggplot(dis_microarray_metabolite, aes(x = Time_Difference*0.001, y = Distance)) +
        geom_line() +
        geom_vline(xintercept = timing_microarray_metabolite*0.001, color = "grey", lty = 2) +
        labs(x = "Time difference (h)", y = "Distance", title = str_glue("The optimal time difference is {timing_microarray_metabolite*0.001*60} minutes")) +
        theme_classic()
    # ggsave(file.path(dir_res, "Distance_dim2-1.pdf"), p, height = 7, width = 10)
    ggsave(file.path(dir_res, "Distance_dim2-1.svg"), p, height = 5, width = 7, device = "svg")
    
} else if (align_method == "distance2") {
    
    # 计算距离，并找到时间差
    dis_microarray_metabolite <- distance2_ts(interp_metabolite$value, interp_microarray$value)
    timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]
    # 绘制距离随时间差变化的图
    p <- ggplot(dis_microarray_metabolite, aes(x = Time_Difference*0.001, y = Distance)) +
        geom_line() +
        geom_vline(xintercept = timing_microarray_metabolite*0.001, color = "grey", lty = 2) +
        labs(x = "Time difference (h)", y = "Distance", title = str_glue("The optimal time difference is {timing_microarray_metabolite*0.001*60} minutes")) +
        theme_classic()
    ggsave(file.path(dir_res, "Distance_dim2-1.pdf"), p, height = 7, width = 10)
    
}

# 绘制最优对齐曲线
ts1 <- interp_metabolite$value
ts2 <- interp_microarray$value
i <- timing_microarray_metabolite
shifted_ts1 <- if (i < 0) {
    c(rep(NA, abs(i)), ts1[1:(length(ts1) + i)])
} else {
    c(ts1[(i + 1):length(ts1)], rep(NA, i))
}
p <- ggplot() +
    geom_line(aes(x = 1:(length(ts1))*0.001, y = shifted_ts1, color = "Metabolome")) +
    geom_line(aes(x = 1:(length(ts1))*0.001, y = ts2, color = "Transcriptome of Tu et al. 2005")) +
    labs(x = "Time (h)", y = "Loading", title = "The Optimal alignment of curves") +
    theme_minimal()
ggsave(file.path(dir_res, "Optimal_alignment_dim2-1.pdf"), p, height = 5, width = 10)


# subsampling -------------------------------------------------------------
# 
# timing_subsampling <- function(i){
# 
#     # data_microarray_subsample <- data_microarray[sample(nrow(data_microarray), 25, replace = F), ]
#     data_microarray_subsample <- data_microarray
#     data_metabolite_subsample <- data_metabolite[sample(nrow(data_metabolite), 17, replace = F), ]
# 
#     # 插值时间点
#     time_min_subsample <- max(min(data_microarray_subsample$Time), min(data_metabolite_subsample$Time), time_min)
#     time_max_subsample <- min(max(data_microarray_subsample$Time), max(data_metabolite_subsample$Time), time_max)
#     timepoints <- seq(from = time_min_subsample, to = time_max_subsample, by = 0.001)
# 
#     
#     # Dim 1
#     data_microarray_1 <- data_microarray_subsample %>%
#         select(time = Time, value = `1`)
#     data_metabolite_1 <- data_metabolite_subsample %>%
#         select(time = Time, value = `1`)
#     
#     # 数据标准化
#     if (is_norm == "Norm") {
#         data_microarray_1$value <- normalize(data_microarray_1$value)
#         data_metabolite_1$value <- normalize(data_metabolite_1$value)
#     } else if (is_norm == "Scale") {
#         data_microarray_1$value <- scalize(data_microarray_1$value)
#         data_metabolite_1$value <- scalize(data_metabolite_1$value)
#     }
#     
#     # 插值计算
#     if (interpolation_method == "Linear") {
#         # linear
#         interp_microarray <- tibble(
#             time = timepoints,
#             value = coredata(na.approx(zoo(data_microarray_1$value, order.by = data_microarray_1$time), xout = timepoints))
#         )
#         interp_metabolite <- tibble(
#             time = timepoints,
#             value = coredata(na.approx(zoo(data_metabolite_1$value, order.by = data_metabolite_1$time), xout = timepoints))
#         )
#     } else if (interpolation_method == "Spline") {
#         # spline
#         interp_microarray <- tibble(
#             time = timepoints,
#             value = spline(data_microarray_1$time, data_microarray_1$value, xout = timepoints)$y
#         )
#         interp_metabolite <- tibble(
#             time = timepoints,
#             value = spline(data_metabolite_1$time, data_metabolite_1$value, xout = timepoints)$y
#         )
#     }
# 
#     # 计算最优对齐
#     if (align_method == "ccf") {
# 
#         # 计算ccf，并找到相位差
#         ccf_microarray_metabolite <- ccf(interp_metabolite$value, interp_microarray$value, lag.max = 1000, plot = F)
#         timing_microarray_metabolite <- ccf_microarray_metabolite$lag[which.max(ccf_microarray_metabolite$acf)]
# 
#     } else if (align_method == "distance1") {
# 
#         # 计算距离，并找到时间差
#         dis_microarray_metabolite <- distance_ts(interp_metabolite$value, interp_microarray$value, lag = 1000)
#         timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]
# 
#     } else if (align_method == "distance2") {
# 
#         # 计算距离，并找到时间差
#         dis_microarray_metabolite <- distance2_ts(interp_metabolite$value, interp_microarray$value, lag = 1000)
#         timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]
# 
#     }
# 
#     res <- c(
#         "dim1" = timing_microarray_metabolite*0.001*60
#     )
# 
# 
#     # Dim 2
#     data_microarray_1 <- data_microarray_subsample %>%
#         select(time = Time, value = `2`)
#     data_metabolite_1 <- data_metabolite_subsample %>%
#         select(time = Time, value = `2`)
#     
#     # 数据标准化
#     if (is_norm == "Norm") {
#         data_microarray_1$value <- normalize(data_microarray_1$value)
#         data_metabolite_1$value <- normalize(data_metabolite_1$value)
#     } else if (is_norm == "Scale") {
#         data_microarray_1$value <- scalize(data_microarray_1$value)
#         data_metabolite_1$value <- scalize(data_metabolite_1$value)
#     }
#     
#     # 插值计算
#     if (interpolation_method == "Linear") {
#         # linear
#         interp_microarray <- tibble(
#             time = timepoints,
#             value = coredata(na.approx(zoo(data_microarray_1$value, order.by = data_microarray_1$time), xout = timepoints))
#         )
#         interp_metabolite <- tibble(
#             time = timepoints,
#             value = coredata(na.approx(zoo(data_metabolite_1$value, order.by = data_metabolite_1$time), xout = timepoints))
#         )
#     } else if (interpolation_method == "Spline") {
#         # spline
#         interp_microarray <- tibble(
#             time = timepoints,
#             value = spline(data_microarray_1$time, data_microarray_1$value, xout = timepoints)$y
#         )
#         interp_metabolite <- tibble(
#             time = timepoints,
#             value = spline(data_metabolite_1$time, data_metabolite_1$value, xout = timepoints)$y
#         )
#     }
# 
#     # 计算最优对齐
#     if (align_method == "ccf") {
# 
#         # 计算ccf，并找到相位差
#         ccf_microarray_metabolite <- ccf(interp_metabolite$value, interp_microarray$value, lag.max = 1000, plot = F)
#         timing_microarray_metabolite <- ccf_microarray_metabolite$lag[which.max(ccf_microarray_metabolite$acf)]
# 
#     } else if (align_method == "distance1") {
# 
#         # 计算距离，并找到时间差
#         dis_microarray_metabolite <- distance_ts(interp_metabolite$value, interp_microarray$value, lag = 1000)
#         timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]
# 
#     } else if (align_method == "distance2") {
# 
#         # 计算距离，并找到时间差
#         dis_microarray_metabolite <- distance2_ts(interp_metabolite$value, interp_microarray$value, lag = 1000)
#         timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]
# 
#     }
# 
#     res <- c(
#         res,
#         "dim2" = timing_microarray_metabolite*0.001*60
#     )
# 
#     return(res)
# }
# # a <- t(sapply(1:100, timing_subsampling))
# sfInit(parallel = TRUE, cpus = CPUs)
# sfExportAll()
# sfLibrary(tidyverse)
# sfLibrary(zoo)
# timing_distribution <- sfSapply(1:10000, timing_subsampling) %>% t()
# sfStop()
# 
# 
# # 置信区间
# 
# # dim1
# x <- timing_distribution[, "dim1"]
# # 计算95%置信区间
# ci_95 <- quantile(x, probs = c(0.025, 0.975))  # 2.5百分位和97.5百分位
# # 绘制密度图
# pdf(file.path(dir_res, "Density_subsample_dim1.pdf"), height = 6, width = 10)
# plot(density(x), main = str_glue("Density of timing; the 95% CI lower and upper: {ci_95[1]} and {ci_95[2]}"), xlab = "Timing (minutes)")
# # 标记95%置信区间的上下界
# abline(v = ci_95[1], col = "red", lty = 2, lwd = 2)  # 下界
# abline(v = ci_95[2], col = "red", lty = 2, lwd = 2)  # 上界
# text(ci_95[1], 0.1, str_glue("95% CI Lower: {ci_95[1]}"), pos = 2, col = "red")  # 下界注释
# text(ci_95[2], 0.1, str_glue("95% CI Upper: {ci_95[2]}"), pos = 4, col = "red")  # 上界注释
# dev.off()
# 
# # dim2
# x <- timing_distribution[, "dim2"]
# # 计算95%置信区间
# ci_95 <- quantile(x, probs = c(0.025, 0.975))  # 2.5百分位和97.5百分位
# # 绘制密度图
# pdf(file.path(dir_res, "Density_subsample_dim2.pdf"), height = 6, width = 10)
# plot(density(x), main = str_glue("Density of timing; the 95% CI lower and upper: {ci_95[1]} and {ci_95[2]}"), xlab = "Timing (minutes)")
# # 标记95%置信区间的上下界
# abline(v = ci_95[1], col = "red", lty = 2, lwd = 2)  # 下界
# abline(v = ci_95[2], col = "red", lty = 2, lwd = 2)  # 上界
# text(ci_95[1], 0.1, str_glue("95% CI Lower: {ci_95[1]}"), pos = 2, col = "red")  # 下界注释
# text(ci_95[2], 0.1, str_glue("95% CI Upper: {ci_95[2]}"), pos = 4, col = "red")  # 上界注释
# dev.off()
# 
# 
# 
# perturbation ------------------------------------------------------------

timing_perturbation <- function(i){

    # 插值时间点
    timepoints <- seq(from = time_min, to = time_max, by = 0.001)

    # Dim 1
    data_microarray_1 <- data_microarray %>%
        select(time = Time, value = `1`)
    data_metabolite_1 <- data_metabolite %>%
        select(time = Time, value = `1`)

    # 添加随机扰动
    # 正态分布扰动，扰动标准差为0.1标准差
    data_microarray_1$value <- data_microarray_1$value +
        rnorm(nrow(data_microarray_1), mean = 0, sd = 0.1 * sd(data_microarray_1$value))
    data_metabolite_1$value <- data_metabolite_1$value +
        rnorm(nrow(data_metabolite_1), mean = 0, sd = 0.1 * sd(data_metabolite_1$value))

    # 插值计算
    if (interpolation_method == "Linear") {
        # linear
        interp_microarray <- tibble(
            time = timepoints,
            value = coredata(na.approx(zoo(data_microarray_1$value, order.by = data_microarray_1$time), xout = timepoints))
        )
        interp_metabolite <- tibble(
            time = timepoints,
            value = coredata(na.approx(zoo(data_metabolite_1$value, order.by = data_metabolite_1$time), xout = timepoints))
        )
    } else if (interpolation_method == "Spline") {
        # spline
        interp_microarray <- tibble(
            time = timepoints,
            value = spline(data_microarray_1$time, data_microarray_1$value, xout = timepoints)$y
        )
        interp_metabolite <- tibble(
            time = timepoints,
            value = spline(data_metabolite_1$time, data_metabolite_1$value, xout = timepoints)$y
        )
    }
    
    # 数据标准化
    if (is_norm == "Norm") {
        interp_microarray$value <- normalize(interp_microarray$value)
        interp_metabolite$value <- normalize(interp_metabolite$value)
    } else if (is_norm == "Scale") {
        interp_microarray$value <- scalize(interp_microarray$value)
        interp_metabolite$value <- scalize(interp_metabolite$value)
    }

    # 计算最优对齐
    if (align_method == "ccf") {

        # 计算ccf，并找到相位差
        ccf_microarray_metabolite <- ccf(interp_metabolite$value, interp_microarray$value, lag.max = 1500, plot = F)
        timing_microarray_metabolite <- ccf_microarray_metabolite$lag[which.max(ccf_microarray_metabolite$acf)]

    } else if (align_method == "distance1") {

        # 计算距离，并找到时间差
        dis_microarray_metabolite <- distance_ts(interp_metabolite$value, interp_microarray$value)
        timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]

    } else if (align_method == "distance2") {

        # 计算距离，并找到时间差
        dis_microarray_metabolite <- distance2_ts(interp_metabolite$value, interp_microarray$value)
        timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]

    }

    res <- c(
        "dim1" = timing_microarray_metabolite*0.001*60
    )


    # Dim 2
    data_microarray_1 <- data_microarray %>%
        select(time = Time, value = `2`)
    data_metabolite_1 <- data_metabolite %>%
        select(time = Time, value = `2`)

    # 添加随机扰动
    # 正态分布扰动，扰动标准差为0.1标准差
    data_microarray_1$value <- data_microarray_1$value +
        rnorm(nrow(data_microarray_1), mean = 0, sd = 0.1 * sd(data_microarray_1$value))
    data_metabolite_1$value <- data_metabolite_1$value +
        rnorm(nrow(data_metabolite_1), mean = 0, sd = 0.1 * sd(data_metabolite_1$value))
    
    # 插值计算
    if (interpolation_method == "Linear") {
        # linear
        interp_microarray <- tibble(
            time = timepoints,
            value = coredata(na.approx(zoo(data_microarray_1$value, order.by = data_microarray_1$time), xout = timepoints))
        )
        interp_metabolite <- tibble(
            time = timepoints,
            value = coredata(na.approx(zoo(data_metabolite_1$value, order.by = data_metabolite_1$time), xout = timepoints))
        )
    } else if (interpolation_method == "Spline") {
        # spline
        interp_microarray <- tibble(
            time = timepoints,
            value = spline(data_microarray_1$time, data_microarray_1$value, xout = timepoints)$y
        )
        interp_metabolite <- tibble(
            time = timepoints,
            value = spline(data_metabolite_1$time, data_metabolite_1$value, xout = timepoints)$y
        )
    }

    # 数据标准化
    if (is_norm == "Norm") {
        interp_microarray$value <- normalize(interp_microarray$value)
        interp_metabolite$value <- normalize(interp_metabolite$value)
    } else if (is_norm == "Scale") {
        interp_microarray$value <- scalize(interp_microarray$value)
        interp_metabolite$value <- scalize(interp_metabolite$value)
    }
    
    # 计算最优对齐
    if (align_method == "ccf") {

        # 计算ccf，并找到相位差
        ccf_microarray_metabolite <- ccf(interp_metabolite$value, interp_microarray$value, lag.max = 1500, plot = F)
        timing_microarray_metabolite <- ccf_microarray_metabolite$lag[which.max(ccf_microarray_metabolite$acf)]

    } else if (align_method == "distance1") {

        # 计算距离，并找到时间差
        dis_microarray_metabolite <- distance_ts(interp_metabolite$value, interp_microarray$value)
        timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]

    } else if (align_method == "distance2") {

        # 计算距离，并找到时间差
        dis_microarray_metabolite <- distance2_ts(interp_metabolite$value, interp_microarray$value)
        timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]

    }

    res <- c(
        res,
        "dim2" = timing_microarray_metabolite*0.001*60
    )

    
    
    # Dim 1-2
    data_microarray_1 <- data_microarray %>% 
        select(time = Time, value = `1`)
    data_metabolite_1 <- data_metabolite %>% 
        select(time = Time, value = `2`) 
    
    # 添加随机扰动
    # 正态分布扰动，扰动标准差为0.1标准差
    data_microarray_1$value <- data_microarray_1$value +
        rnorm(nrow(data_microarray_1), mean = 0, sd = 0.1 * sd(data_microarray_1$value))
    data_metabolite_1$value <- data_metabolite_1$value +
        rnorm(nrow(data_metabolite_1), mean = 0, sd = 0.1 * sd(data_metabolite_1$value))
    
    # 插值计算
    if (interpolation_method == "Linear") {
        # linear
        interp_microarray <- tibble(
            time = timepoints,
            value = coredata(na.approx(zoo(data_microarray_1$value, order.by = data_microarray_1$time), xout = timepoints))
        )
        interp_metabolite <- tibble(
            time = timepoints,
            value = coredata(na.approx(zoo(data_metabolite_1$value, order.by = data_metabolite_1$time), xout = timepoints))
        )
    } else if (interpolation_method == "Spline") {
        # spline
        interp_microarray <- tibble(
            time = timepoints,
            value = spline(data_microarray_1$time, data_microarray_1$value, xout = timepoints)$y
        )
        interp_metabolite <- tibble(
            time = timepoints,
            value = spline(data_metabolite_1$time, data_metabolite_1$value, xout = timepoints)$y
        )
    }
    
    # 数据标准化
    if (is_norm == "Norm") {
        interp_microarray$value <- normalize(interp_microarray$value)
        interp_metabolite$value <- normalize(interp_metabolite$value)
    } else if (is_norm == "Scale") {
        interp_microarray$value <- scalize(interp_microarray$value)
        interp_metabolite$value <- scalize(interp_metabolite$value)
    }
    
    # 计算最优对齐
    if (align_method == "ccf") {
        
        # 计算ccf，并找到相位差
        ccf_microarray_metabolite <- ccf(interp_metabolite$value, interp_microarray$value, lag.max = 1500, plot = F)
        timing_microarray_metabolite <- ccf_microarray_metabolite$lag[which.max(ccf_microarray_metabolite$acf)]
        
    } else if (align_method == "distance1") {
        
        # 计算距离，并找到时间差
        dis_microarray_metabolite <- distance_ts(interp_metabolite$value, interp_microarray$value, lag = 3000)
        timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]
        
    } else if (align_method == "distance2") {
        
        # 计算距离，并找到时间差
        dis_microarray_metabolite <- distance2_ts(interp_metabolite$value, interp_microarray$value)
        timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]
        
    }
    
    res <- c(
        res,
        "dim1-2" = timing_microarray_metabolite*0.001*60
    )
    
    
    # Dim 2-1
    data_microarray_1 <- data_microarray %>% 
        select(time = Time, value = `2`)
    data_metabolite_1 <- data_metabolite %>% 
        select(time = Time, value = `1`) %>% 
        mutate(value = -value)
    
    # 添加随机扰动
    # 正态分布扰动，扰动标准差为0.1标准差
    data_microarray_1$value <- data_microarray_1$value +
        rnorm(nrow(data_microarray_1), mean = 0, sd = 0.1 * sd(data_microarray_1$value))
    data_metabolite_1$value <- data_metabolite_1$value +
        rnorm(nrow(data_metabolite_1), mean = 0, sd = 0.1 * sd(data_metabolite_1$value))
    
    # 插值计算
    if (interpolation_method == "Linear") {
        # linear
        interp_microarray <- tibble(
            time = timepoints,
            value = coredata(na.approx(zoo(data_microarray_1$value, order.by = data_microarray_1$time), xout = timepoints))
        )
        interp_metabolite <- tibble(
            time = timepoints,
            value = coredata(na.approx(zoo(data_metabolite_1$value, order.by = data_metabolite_1$time), xout = timepoints))
        )
    } else if (interpolation_method == "Spline") {
        # spline
        interp_microarray <- tibble(
            time = timepoints,
            value = spline(data_microarray_1$time, data_microarray_1$value, xout = timepoints)$y
        )
        interp_metabolite <- tibble(
            time = timepoints,
            value = spline(data_metabolite_1$time, data_metabolite_1$value, xout = timepoints)$y
        )
    }
    
    # 数据标准化
    if (is_norm == "Norm") {
        interp_microarray$value <- normalize(interp_microarray$value)
        interp_metabolite$value <- normalize(interp_metabolite$value)
    } else if (is_norm == "Scale") {
        interp_microarray$value <- scalize(interp_microarray$value)
        interp_metabolite$value <- scalize(interp_metabolite$value)
    }
    
    # 计算最优对齐
    if (align_method == "ccf") {
        
        # 计算ccf，并找到相位差
        ccf_microarray_metabolite <- ccf(interp_metabolite$value, interp_microarray$value, lag.max = 1500, plot = F)
        timing_microarray_metabolite <- ccf_microarray_metabolite$lag[which.max(ccf_microarray_metabolite$acf)]
        
    } else if (align_method == "distance1") {
        
        # 计算距离，并找到时间差
        dis_microarray_metabolite <- distance_ts(interp_metabolite$value, interp_microarray$value)
        timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]
        
    } else if (align_method == "distance2") {
        
        # 计算距离，并找到时间差
        dis_microarray_metabolite <- distance2_ts(interp_metabolite$value, interp_microarray$value)
        timing_microarray_metabolite <- dis_microarray_metabolite$Time_Difference[which.min(dis_microarray_metabolite$Distance)]
        
    }
    
    res <- c(
        res,
        "dim2-1" = timing_microarray_metabolite*0.001*60
    )
    
    
    return(res)
}

sfInit(parallel = TRUE, cpus = CPUs)
sfExportAll()
sfLibrary(tidyverse)
sfLibrary(zoo)
timing_distribution <- sfSapply(1:10000, timing_perturbation) %>% t()
sfStop()


# 置信区间

# dim1
x <- timing_distribution[, "dim1"]
# 计算95%置信区间
ci_95 <- quantile(x, probs = c(0.025, 0.975))  # 2.5百分位和97.5百分位
# 绘制密度图
svg(file.path(dir_res, "Density_perturb_dim1.svg"), height = 6, width = 10)
plot(density(x), main = str_glue("Density of time difference; the 95% CI lower and upper: {ci_95[1]} and {ci_95[2]}"), xlab = "Timing (minutes)")
# 标记95%置信区间的上下界
abline(v = ci_95[1], col = "red", lty = 2, lwd = 2)  # 下界
abline(v = ci_95[2], col = "red", lty = 2, lwd = 2)  # 上界
text(ci_95[1], 0.1, str_glue("95% CI Lower"), pos = 2, col = "red")  # 下界注释
text(ci_95[2], 0.1, str_glue("95% CI Upper"), pos = 4, col = "red")  # 上界注释
dev.off()

# dim2
x <- timing_distribution[, "dim2"]
# 计算95%置信区间
ci_95 <- quantile(x, probs = c(0.025, 0.975))  # 2.5百分位和97.5百分位
# 绘制密度图
svg(file.path(dir_res, "Density_perturb_dim2.svg"), height = 6, width = 10)
plot(density(x), main = str_glue("Density of time difference; the 95% CI lower and upper: {ci_95[1]} and {ci_95[2]}"), xlab = "Timing (minutes)")
# 标记95%置信区间的上下界
abline(v = ci_95[1], col = "red", lty = 2, lwd = 2)  # 下界
abline(v = ci_95[2], col = "red", lty = 2, lwd = 2)  # 上界
text(ci_95[1], 0.1, str_glue("95% CI Lower"), pos = 2, col = "red")  # 下界注释
text(ci_95[2], 0.1, str_glue("95% CI Upper"), pos = 4, col = "red")  # 上界注释
dev.off()


# dim1-2
x <- timing_distribution[, "dim1-2"]
# 计算95%置信区间
ci_95 <- quantile(x, probs = c(0.025, 0.975))  # 2.5百分位和97.5百分位
# 绘制密度图
svg(file.path(dir_res, "Density_perturb_dim1-2.svg"), height = 6, width = 10)
plot(density(x), main = str_glue("Density of time difference; the 95% CI lower and upper: {ci_95[1]} and {ci_95[2]}"), xlab = "Timing (minutes)")
# 标记95%置信区间的上下界
abline(v = ci_95[1], col = "red", lty = 2, lwd = 2)  # 下界
abline(v = ci_95[2], col = "red", lty = 2, lwd = 2)  # 上界
text(ci_95[1], 0.1, str_glue("95% CI Lower"), pos = 2, col = "red")  # 下界注释
text(ci_95[2], 0.1, str_glue("95% CI Upper"), pos = 4, col = "red")  # 上界注释
dev.off()

# dim2-1
x <- timing_distribution[, "dim2-1"]
# 计算95%置信区间
ci_95 <- quantile(x, probs = c(0.025, 0.975))  # 2.5百分位和97.5百分位
# 绘制密度图
svg(file.path(dir_res, "Density_perturb_dim2-1.svg"), height = 6, width = 10)
plot(density(x), main = str_glue("Density of time difference; the 95% CI lower and upper: {ci_95[1]} and {ci_95[2]}"), xlab = "Timing (minutes)")
# 标记95%置信区间的上下界
abline(v = ci_95[1], col = "red", lty = 2, lwd = 2)  # 下界
abline(v = ci_95[2], col = "red", lty = 2, lwd = 2)  # 上界
text(ci_95[1], 0.1, str_glue("95% CI Lower"), pos = 2, col = "red")  # 下界注释
text(ci_95[2], 0.1, str_glue("95% CI Upper"), pos = 4, col = "red")  # 上界注释
dev.off()
