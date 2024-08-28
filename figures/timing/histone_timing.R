suppressMessages({
    library(snowfall)
    library(zoo)
    library(svglite)
    library(rlist)
    library(tidyverse)
})

setwd("Z:/home/wanglinting/Yeast/CodeOcean/figures/timing")

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

dir_res <- str_glue("./histone_{time_align}_{is_norm}_{interpolation_method}_{align_method}_{time_range}")
if (!dir.exists(dir_res)) dir.create(dir_res, recursive = TRUE)


# functions ---------------------------------------------------------------

if (time_align == "dtw") {
    # dtw
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
            T1 = 0.345,
            T2 = 0.610,
            T3 = 0.745,
            T4 = 0.888,
            T5 = 1.106,
            T6 = 1.362,
            T7 = 1.616,
            T8 = 1.780,
            T9 = 1.912,
            T10 = 2.031,
            T11 = 2.258,
            T12 = 2.400,
            T13 = 2.886,
            T14 = 3.251,
            T15 = 3.647,
            T16 = 4.017
        )
    }
} else if (time_align == "LinearTransform") {
    # LinearTransform
    getRNASeqTime <- function(x) {
        switch(
            x,
            T1 = 50.7556,
            T2 = 51.0284,
            T3 = 51.134,
            T4 = 51.2044,
            T5 = 51.3364,
            T6 = 51.486,
            T7 = 51.6356,
            T8 = 51.7764,
            T9 = 51.926,
            T10 = 52.0756,
            T11 = 52.2164,
            T12 = 52.3748,
            T13 = 52.5244,
            T14 = 52.9996,
            T15 = 53.466,
            T16 = 53.95
        )
    }
    getChipSeqTime <- function(x) {
        switch(
            x,
            T1 = 50.5,
            T2 = 50.9,
            T3 = 51.08,
            T4 = 51.18,
            T5 = 51.34,
            T6 = 51.51,
            T7 = 51.68,
            T8 = 51.84,
            T9 = 52.01,
            T10 = 52.18,
            T11 = 52.38,
            T12 = 52.54,
            T13 = 52.97,
            T14 = 53.42,
            T15 = 53.83,
            T16 = 54.25
        )
    }
    
} else if (time_align == "ddtw") {
    
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
    
}


# 数据正规化
scalize <- function(x) {
    # 标准化
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)  
}
normalize <- function(x) {
    # 归一化
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}


# 计算序列的距离
# 使用mapply函数优化计算
distance_ts <- function(ts1, ts2, lag = 1000) {
    
    calculate_distance <- function(i, ts1, ts2) {
        # 根据时间差i移动时间序列ts1
        shifted_ts1 <- if (i < 0) {
            c(rep(NA, abs(i)), ts1[1:(length(ts1) + i)])
        } else {
            c(ts1[(i + 1):length(ts1)], rep(NA, i))
        }

        # 计算并返回移动后的距离
        mean(abs(shifted_ts1 - ts2), na.rm = TRUE)
    }

    # 确定实际的最大时间差
    l <- min(lag, length(ts1))
    
    # 使用mapply对时间差范围内的每一个值应用calculate_distance函数
    distances <- mapply(calculate_distance, -l:l, MoreArgs = list(ts1 = ts1, ts2 = ts2))
    
    # 创建一个数据框存储结果
    distance_df <- data.frame(Time_Difference = -l:l, Distance = distances)
    
    return(distance_df)
}

distance2_ts <- function(ts1, ts2, lag = 1000) {
    
    calculate_distance <- function(i, ts1, ts2) {
        # 根据时间差i移动时间序列ts1
        shifted_ts1 <- if (i < 0) {
            c(rep(NA, abs(i)), ts1[1:(length(ts1) + i)])
        } else {
            c(ts1[(i + 1):length(ts1)], rep(NA, i))
        }
        
        # ts2只保留固定时间点
        # 特定的时间点，包括超出范围的值
        specific_times <- c(
            0.356,
            0.666,
            0.782,
            0.866,
            1.017,
            1.187,
            1.357,
            1.516,
            1.686,
            1.856,
            2.031,
            2.196,
            2.366,
            2.908,
            3.437,
            3.987
        )  # 包括超出timepoints范围的时间点
        
        # 创建结果向量
        selected_ts2 <- rep(NA, length(ts2))  # 初始化结果向量，全部值设为NA
        # 找出特定时间点在timepoints中的位置
        indices <- match(specific_times, timepoints)
        # 从a中提取对应特定时间点的值，仅当索引有效时（不是NA）
        valid_indices <- !is.na(indices)
        selected_ts2[indices[valid_indices]] <- ts2[indices[valid_indices]]
        
        # 计算并返回移动后的距离
        mean(abs(shifted_ts1 - selected_ts2), na.rm = TRUE)
    }
    
    # 确定实际的最大时间差
    l <- min(lag, length(ts1))
    
    # 使用mapply对时间差范围内的每一个值应用calculate_distance函数
    distances <- mapply(calculate_distance, -l:l, MoreArgs = list(ts1 = ts1, ts2 = ts2))
    
    # 创建一个数据框存储结果
    distance_df <- data.frame(Time_Difference = -l:l, Distance = distances)
    
    return(distance_df)
}

# read data ---------------------------------------------------------------

# rnaseq data
file_rnaseq_sample <- "../../results/temp/rnaseq14/sample_name.txt"
file_rnaseq_v <- "../../results/temp/rnaseq14/svd/V.txt"
data_rnaseq_v <- read_table(file_rnaseq_v,
                            col_names = as.character(0:(number_dim - 1)),
                            show_col_types = F)
data_rnaseq <- read_table(file_rnaseq_sample,
                          col_names = "Sample",
                          show_col_types = F) %>%
    mutate(Time = sapply(Sample, getRNASeqTime)) %>%
    bind_cols(data_rnaseq_v)

# H3K9ac
file_H3K9ac_sample <- "../../results/temp/H3K9ac/sample_name.txt"
file_H3K9ac_v <- "../../results/temp/H3K9ac/svd/V.txt"
data_H3K9ac_v <-  read_table(file_H3K9ac_v,
                             col_names = as.character(0:(number_dim - 1)),
                             show_col_types = F)
data_H3K9ac <- read_table(file_H3K9ac_sample,
                          col_names = "Sample",
                          show_col_types = F) %>%
    mutate(Time = sapply(Sample, getChipSeqTime)) %>%
    bind_cols(data_H3K9ac_v)


# H3K18ac
file_H3K18ac_sample <- "../../results/temp/H3K18ac/sample_name.txt"
file_H3K18ac_v <- "../../results/temp/H3K18ac/svd/V.txt"
data_H3K18ac_v <-  read_table(file_H3K18ac_v,
                              col_names = as.character(0:(number_dim - 1)),
                              show_col_types = F)
data_H3K18ac <- read_table(file_H3K18ac_sample,
                           col_names = "Sample",
                           show_col_types = F) %>%
    mutate(Time = sapply(Sample, getChipSeqTime)) %>%
    bind_cols(data_H3K18ac_v)


# 插值时间点
if (time_range == "all") {
    time_min <- max(min(data_rnaseq$Time), min(data_H3K9ac$Time), min(data_H3K18ac$Time))
    time_max <- min(max(data_rnaseq$Time), max(data_H3K9ac$Time), max(data_H3K18ac$Time))
} else if (time_range == "sub1") {
    time_min <- max(min(data_rnaseq$Time), min(data_H3K9ac$Time), min(data_H3K18ac$Time), 0.5)
    time_max <- min(max(data_rnaseq$Time), max(data_H3K9ac$Time), max(data_H3K18ac$Time), 2.5)
} else if (time_range == "sub2") {
    time_min <- max(min(data_rnaseq$Time), min(data_H3K9ac$Time), min(data_H3K18ac$Time), 0.5)
    time_max <- min(max(data_rnaseq$Time), max(data_H3K9ac$Time), max(data_H3K18ac$Time), 2)
} else if (time_range == "sub3") {
    time_min <- max(min(data_rnaseq$Time), min(data_H3K9ac$Time), min(data_H3K18ac$Time), 0.8)
    time_max <- min(max(data_rnaseq$Time), max(data_H3K9ac$Time), max(data_H3K18ac$Time), 1.3)
}

timepoints <- seq(from = time_min, to = time_max, by = 0.001)


# dim 1 -------------------------------------------------------------------

data_rnaseq_1 <- data_rnaseq %>% 
    select(time = Time, value = `1`)
data_H3K9ac_1 <- data_H3K9ac %>% 
    select(time = Time, value = `1`)
data_H3K18ac_1 <- data_H3K18ac %>% 
    select(time = Time, value = `1`)

# 插值计算
if (interpolation_method == "Linear") {
    # linear
    interp_rnaseq <- tibble(
        time = timepoints,
        value = coredata(na.approx(zoo(data_rnaseq_1$value, order.by = data_rnaseq_1$time), xout = timepoints))
    )
    interp_H3K9ac <- tibble(
        time = timepoints,
        value = coredata(na.approx(zoo(data_H3K9ac_1$value, order.by = data_H3K9ac_1$time), xout = timepoints))
    )
    interp_H3K18ac <- tibble(
        time = timepoints,
        value = coredata(na.approx(zoo(data_H3K18ac_1$value, order.by = data_H3K18ac_1$time), xout = timepoints))
    )
} else if (interpolation_method == "Spline") {
    # spline
    interp_rnaseq <- tibble(
        time = timepoints,
        value = spline(data_rnaseq_1$time, data_rnaseq_1$value, xout = timepoints)$y 
    )
    interp_H3K9ac <- tibble( 
        time = timepoints,
        value = spline(data_H3K9ac_1$time, data_H3K9ac_1$value, xout = timepoints)$y
    )
    interp_H3K18ac <- tibble(
        time = timepoints,
        value = spline(data_H3K18ac_1$time, data_H3K18ac_1$value, xout = timepoints)$y
    )
}

# 数据标准化
if (is_norm == "Norm") {
    interp_rnaseq$value <- normalize(interp_rnaseq$value)
    interp_H3K9ac$value <- normalize(interp_H3K9ac$value)
    interp_H3K18ac$value <- normalize(interp_H3K18ac$value)
}  else if (is_norm == "Scale") {
    interp_rnaseq$value <- normalize(interp_rnaseq$value)
    interp_H3K9ac$value <- normalize(interp_H3K9ac$value)
    interp_H3K18ac$value <- normalize(interp_H3K18ac$value)
}

# 画出插值结果
p <- ggplot() +
    geom_line(data = interp_rnaseq, aes(x = time, y = value, color = "Transcriptome of Zheng et al. 2014")) +
    geom_line(data = interp_H3K9ac, aes(x = time, y = value, color = "Epigenome of histone modification H3K9ac")) +
    geom_line(data = interp_H3K18ac, aes(x = time, y = value, color = "Epigenome of histone modification H3K18ac")) +
    labs(title = str_glue("{interpolation_method} interpolation"), x = "Time (h)", y = "Loading") +
    theme_minimal()
ggsave(file.path(dir_res, "Interpolation_dim1.pdf"), p, height = 5, width = 10)

# 计算最优对齐
if (align_method == "ccf") {

    # H3K9ac
    # 计算ccf，并找到相位差
    ccf_rnaseq_H3K9ac <- ccf(interp_H3K9ac$value, interp_rnaseq$value, lag.max = 1000, plot = F)
    timing_rnaseq_H3K9ac <- ccf_rnaseq_H3K9ac$lag[which.max(ccf_rnaseq_H3K9ac$acf)]
    # plot
    # pdf(file.path(dir_res, "CCF_H3K9ac_dim1.pdf"), height = 7, width = 10)
    svg(file.path(dir_res, "CCF_H3K9ac_dim1.svg"), height = 5, width = 7)
    plot(ccf_rnaseq_H3K9ac, main = str_glue("The optimal time difference is {timing_rnaseq_H3K9ac*0.001*60} minutes"), xlab = "Lag (0.06 minutes)"  )
    abline(v = timing_rnaseq_H3K9ac, col = "grey", lty = 2)
    dev.off()

    # H3K18ac
    ccf_rnaseq_H3K18ac <- ccf(interp_H3K18ac$value, interp_rnaseq$value, lag.max = 1000, plot = F)
    timing_rnaseq_H3K18ac <- ccf_rnaseq_H3K18ac$lag[which.max(ccf_rnaseq_H3K18ac$acf)]
    # plot
    # pdf(file.path(dir_res, "CCF_H3K18ac_dim1.pdf"), height = 7, width = 10)
    svg(file.path(dir_res, "CCF_H3K18ac_dim1.svg"), height = 5, width = 7)
    plot(ccf_rnaseq_H3K18ac, main = str_glue("The optimal time difference is {timing_rnaseq_H3K18ac*0.001*60} minutes"), xlab = "Lag (0.06 minutes)"  )
    abline(v = timing_rnaseq_H3K9ac, col = "grey", lty = 2)
    dev.off()

} else if (align_method == "distance1") {

    # H3K9ac
    # 计算距离，并找到时间差
    dis_rnaseq_H3K9ac <- distance_ts(interp_H3K9ac$value, interp_rnaseq$value, lag = 1000)
    timing_rnaseq_H3K9ac <- dis_rnaseq_H3K9ac$Time_Difference[which.min(dis_rnaseq_H3K9ac$Distance)]
    # 绘制距离随时间差变化的图
    p <- ggplot(dis_rnaseq_H3K9ac, aes(x = Time_Difference*0.001, y = Distance)) +
        geom_line() +
        geom_vline(xintercept = timing_rnaseq_H3K9ac*0.001, color = "grey", lty = 2) +
        labs(x = "Time difference (h)", y = "Distance", title = str_glue("The optimal time difference is {timing_rnaseq_H3K9ac*0.001*60} minutes")) +
        theme_classic()
    # ggsave(file.path(dir_res, "Distance_H3K9ac_dim1.pdf"), p, height = 7, width = 10)
    ggsave(file.path(dir_res, "Distance_H3K9ac_dim1.svg"), p, height = 5, width = 7, device = "svg")
    
    # H3K18ac
    dis_rnaseq_H3K18ac <- distance_ts(interp_H3K18ac$value, interp_rnaseq$value, lag = 1000)
    timing_rnaseq_H3K18ac <- dis_rnaseq_H3K18ac$Time_Difference[which.min(dis_rnaseq_H3K18ac$Distance)]
    p <- ggplot(dis_rnaseq_H3K18ac, aes(x = Time_Difference*0.001, y = Distance)) +
        geom_line() +
        geom_vline(xintercept = timing_rnaseq_H3K18ac*0.001, color = "grey", lty = 2) +
        labs(x = "Time difference (h)", y = "Distance", title = str_glue("The optimal time difference is {timing_rnaseq_H3K18ac*0.001*60} minutes")) +
        theme_classic()
    # ggsave(file.path(dir_res, "Distance_H3K18ac_dim1.pdf"), p, height = 7, width = 10)
    ggsave(file.path(dir_res, "Distance_H3K18ac_dim1.svg"), p, height = 5, width = 7, device = "svg")

} else if (align_method == "distance2") {
    
    # H3K9ac
    # 计算距离，并找到时间差
    dis_rnaseq_H3K9ac <- distance2_ts(interp_H3K9ac$value, interp_rnaseq$value, lag = 1000)
    timing_rnaseq_H3K9ac <- dis_rnaseq_H3K9ac$Time_Difference[which.min(dis_rnaseq_H3K9ac$Distance)]
    # 绘制距离随时间差变化的图
    p <- ggplot(dis_rnaseq_H3K9ac, aes(x = Time_Difference*0.001, y = Distance)) +
        geom_line() +
        geom_vline(xintercept = timing_rnaseq_H3K9ac*0.001, color = "grey", lty = 2) +
        labs(x = "Time Difference (h)", y = "Distance", title = str_glue("The timing is {timing_rnaseq_H3K9ac*0.001*60} minutes")) +
        theme_classic()
    ggsave(file.path(dir_res, "Distance_H3K9ac_dim1.pdf"), p, height = 7, width = 10)
    
    # H3K18ac
    dis_rnaseq_H3K18ac <- distance2_ts(interp_H3K18ac$value, interp_rnaseq$value, lag = 1000)
    timing_rnaseq_H3K18ac <- dis_rnaseq_H3K18ac$Time_Difference[which.min(dis_rnaseq_H3K18ac$Distance)]
    p <- ggplot(dis_rnaseq_H3K18ac, aes(x = Time_Difference*0.001, y = Distance)) +
        geom_line() +
        geom_vline(xintercept = timing_rnaseq_H3K18ac*0.001, color = "grey", lty = 2) +
        labs(x = "Time Difference (h)", y = "Distance", title = str_glue("The timing is {timing_rnaseq_H3K18ac*0.001*60} minutes")) +
        theme_classic()
    ggsave(file.path(dir_res, "Distance_H3K18ac_dim1.pdf"), p, height = 7, width = 10)
    
}

# 绘制最小距离时的曲线 H3K9ac
ts1 <- interp_H3K9ac$value
ts2 <- interp_rnaseq$value
i <- timing_rnaseq_H3K9ac
shifted_ts1 <- if (i < 0) {
    c(rep(NA, abs(i)), ts1[1:(length(ts1) + i)])
} else {
    c(ts1[(i + 1):length(ts1)], rep(NA, i))
}
p <- ggplot() +
    geom_line(aes(x = 1:(length(ts1))*0.001, y = shifted_ts1, color = "Epigenome of histone modification H3K9ac")) +
    geom_line(aes(x = 1:(length(ts1))*0.001, y = ts2, color = "Transcriptome of Zheng et al. 2014")) +
    labs(x = "Time (h)", y = "Loading", title = "The Optimal alignment of curves") +
    theme_classic()
ggsave(file.path(dir_res, "Optimal_alignment_H3K9ac_dim1.pdf"), p, height = 5, width = 10)


# 绘制最小距离时的曲线 H3K18ac
ts1 <- interp_H3K18ac$value
ts2 <- interp_rnaseq$value
i <- timing_rnaseq_H3K18ac
shifted_ts1 <- if (i < 0) {
    c(rep(NA, abs(i)), ts1[1:(length(ts1) + i)])
} else {
    c(ts1[(i + 1):length(ts1)], rep(NA, i))
}
p <- ggplot() +
    geom_line(aes(x = 1:(length(ts1))*0.001, y = shifted_ts1, color = "Epigenome of histone modification H3K18ac")) +
    geom_line(aes(x = 1:(length(ts1))*0.001, y = ts2, color = "Transcriptome of Zheng et al. 2014")) +
    labs(x = "Time (h)", y = "Loading", title = "The Optimal alignment of curves") +
    theme_minimal()
ggsave(file.path(dir_res, "Optimal_alignment_H3K18ac_dim1.pdf"), p, height = 5, width = 10)


# dim 2 -------------------------------------------------------------------

data_rnaseq_1 <- data_rnaseq %>% 
    select(time = Time, value = `2`)
data_H3K9ac_1 <- data_H3K9ac %>% 
    select(time = Time, value = `2`)
data_H3K18ac_1 <- data_H3K18ac %>% 
    select(time = Time, value = `2`)

# 插值计算
if (interpolation_method == "Linear") {
    # linear
    interp_rnaseq <- tibble(
        time = timepoints,
        value = coredata(na.approx(zoo(data_rnaseq_1$value, order.by = data_rnaseq_1$time), xout = timepoints))
    )
    interp_H3K9ac <- tibble(
        time = timepoints,
        value = coredata(na.approx(zoo(data_H3K9ac_1$value, order.by = data_H3K9ac_1$time), xout = timepoints))
    )
    interp_H3K18ac <- tibble(
        time = timepoints,
        value = coredata(na.approx(zoo(data_H3K18ac_1$value, order.by = data_H3K18ac_1$time), xout = timepoints))
    )
} else if (interpolation_method == "Spline") {
    # spline
    interp_rnaseq <- tibble(
        time = timepoints,
        value = spline(data_rnaseq_1$time, data_rnaseq_1$value, xout = timepoints)$y 
    )
    interp_H3K9ac <- tibble( 
        time = timepoints,
        value = spline(data_H3K9ac_1$time, data_H3K9ac_1$value, xout = timepoints)$y
    )
    interp_H3K18ac <- tibble(
        time = timepoints,
        value = spline(data_H3K18ac_1$time, data_H3K18ac_1$value, xout = timepoints)$y
    )
}

# 数据标准化
if (is_norm == "Norm") {
    interp_rnaseq$value <- normalize(interp_rnaseq$value)
    interp_H3K9ac$value <- normalize(interp_H3K9ac$value)
    interp_H3K18ac$value <- normalize(interp_H3K18ac$value)
}  else if (is_norm == "Scale") {
    interp_rnaseq$value <- normalize(interp_rnaseq$value)
    interp_H3K9ac$value <- normalize(interp_H3K9ac$value)
    interp_H3K18ac$value <- normalize(interp_H3K18ac$value)
}

# 画出插值结果
p <- ggplot() +
    geom_line(data = interp_rnaseq, aes(x = time, y = value, color = "Transcriptome of Zheng et al. 2014")) +
    geom_line(data = interp_H3K9ac, aes(x = time, y = value, color = "Epigenome of histone modification H3K9ac")) +
    geom_line(data = interp_H3K18ac, aes(x = time, y = value, color = "Epigenome of histone modification H3K18ac")) +
    labs(title = str_glue("{interpolation_method} interpolation"), x = "Time (h)", y = "Loading") +
    theme_minimal()
ggsave(file.path(dir_res, "Interpolation_dim2.pdf"), p, height = 5, width = 10)

# 计算最优对齐
if (align_method == "ccf") {

    # H3K9ac
    # 计算ccf，并找到相位差
    ccf_rnaseq_H3K9ac <- ccf(interp_H3K9ac$value, interp_rnaseq$value, lag.max = 1000, plot = F)
    timing_rnaseq_H3K9ac <- ccf_rnaseq_H3K9ac$lag[which.max(ccf_rnaseq_H3K9ac$acf)]
    # plot
    # pdf(file.path(dir_res, "CCF_H3K9ac_dim2.pdf"), height = 7, width = 10)
    svg(file.path(dir_res, "CCF_H3K9ac_dim2.svg"), height = 5, width = 7)
    plot(ccf_rnaseq_H3K9ac, main = str_glue("The optimal time difference is {timing_rnaseq_H3K9ac*0.001*60} minutes"), xlab = "Lag (0.06 minutes)"  )
    abline(v = timing_rnaseq_H3K9ac, col = "grey", lty = 2)
    dev.off()

    # H3K18ac
    ccf_rnaseq_H3K18ac <- ccf(interp_H3K18ac$value, interp_rnaseq$value, lag.max = 1000, plot = F)
    timing_rnaseq_H3K18ac <- ccf_rnaseq_H3K18ac$lag[which.max(ccf_rnaseq_H3K18ac$acf)]
    # plot
    # pdf(file.path(dir_res, "CCF_H3K18ac_dim2.pdf"), height = 7, width = 10)
    svg(file.path(dir_res, "CCF_H3K18ac_dim2.svg"), height = 5, width = 7)
    plot(ccf_rnaseq_H3K18ac, main = str_glue("The optimal time difference is {timing_rnaseq_H3K18ac*0.001*60} minutes"), xlab = "Lag (0.06 minutes)"  )
    abline(v = timing_rnaseq_H3K9ac, col = "grey", lty = 2)
    dev.off()
    
} else if (align_method == "distance1") {
    
    # H3K9ac
    # 计算距离，并找到时间差
    dis_rnaseq_H3K9ac <- distance_ts(interp_H3K9ac$value, interp_rnaseq$value, lag = 1000)
    timing_rnaseq_H3K9ac <- dis_rnaseq_H3K9ac$Time_Difference[which.min(dis_rnaseq_H3K9ac$Distance)]
    # 绘制距离随时间差变化的图
    p <- ggplot(dis_rnaseq_H3K9ac, aes(x = Time_Difference*0.001, y = Distance)) +
        geom_line() +
        geom_vline(xintercept = timing_rnaseq_H3K9ac*0.001, color = "grey", lty = 2) +
        labs(x = "Time difference (h)", y = "Distance", title = str_glue("The optimal time difference is {timing_rnaseq_H3K9ac*0.001*60} minutes")) +
        theme_classic()
    # ggsave(file.path(dir_res, "Distance_H3K9ac_dim2.pdf"), p, height = 7, width = 10)
    ggsave(file.path(dir_res, "Distance_H3K9ac_dim2.svg"), p, height = 5, width = 7, device = "svg")
    
    # H3K18ac
    dis_rnaseq_H3K18ac <- distance_ts(interp_H3K18ac$value, interp_rnaseq$value, lag = 1000)
    timing_rnaseq_H3K18ac <- dis_rnaseq_H3K18ac$Time_Difference[which.min(dis_rnaseq_H3K18ac$Distance)]
    p <- ggplot(dis_rnaseq_H3K18ac, aes(x = Time_Difference*0.001, y = Distance)) +
        geom_line() +
        geom_vline(xintercept = timing_rnaseq_H3K18ac*0.001, color = "grey", lty = 2) +
        labs(x = "Time difference (h)", y = "Distance", title = str_glue("The optimal time difference is {timing_rnaseq_H3K18ac*0.001*60} minutes")) +
        theme_classic()
    # ggsave(file.path(dir_res, "Distance_H3K18ac_dim2.pdf"), p, height = 7, width = 10)
    ggsave(file.path(dir_res, "Distance_H3K18ac_dim2.svg"), p, height = 5, width = 7, device = "svg")
    
    
} else if (align_method == "distance2") {
    
    # H3K9ac
    # 计算距离，并找到时间差
    dis_rnaseq_H3K9ac <- distance2_ts(interp_H3K9ac$value, interp_rnaseq$value, lag = 1000)
    timing_rnaseq_H3K9ac <- dis_rnaseq_H3K9ac$Time_Difference[which.min(dis_rnaseq_H3K9ac$Distance)]
    # 绘制距离随时间差变化的图
    p <- ggplot(dis_rnaseq_H3K9ac, aes(x = Time_Difference*0.001, y = Distance)) +
        geom_line() +
        geom_vline(xintercept = timing_rnaseq_H3K9ac*0.001, color = "grey", lty = 2) +
        labs(x = "Time Difference (h)", y = "Distance", title = str_glue("The timing is {timing_rnaseq_H3K9ac*0.001*60} minutes")) +
        theme_classic()
    ggsave(file.path(dir_res, "Distance_H3K9ac_dim2.pdf"), p, height = 7, width = 10)
    
    # H3K18ac
    dis_rnaseq_H3K18ac <- distance2_ts(interp_H3K18ac$value, interp_rnaseq$value, lag = 1000)
    timing_rnaseq_H3K18ac <- dis_rnaseq_H3K18ac$Time_Difference[which.min(dis_rnaseq_H3K18ac$Distance)]
    p <- ggplot(dis_rnaseq_H3K18ac, aes(x = Time_Difference*0.001, y = Distance)) +
        geom_line() +
        geom_vline(xintercept = timing_rnaseq_H3K18ac*0.001, color = "grey", lty = 2) +
        labs(x = "Time Difference (h)", y = "Distance", title = str_glue("The timing is {timing_rnaseq_H3K18ac*0.001*60} minutes")) +
        theme_classic()
    ggsave(file.path(dir_res, "Distance_H3K18ac_dim2.pdf"), p, height = 7, width = 10)
    
}

# 绘制最小距离时的曲线 H3K9ac
ts1 <- interp_H3K9ac$value
ts2 <- interp_rnaseq$value
i <- timing_rnaseq_H3K9ac
shifted_ts1 <- if (i < 0) {
    c(rep(NA, abs(i)), ts1[1:(length(ts1) + i)])
} else {
    c(ts1[(i + 1):length(ts1)], rep(NA, i))
}
p <- ggplot() +
    geom_line(aes(x = 1:(length(ts1))*0.001, y = shifted_ts1, color = "Epigenome of histone modification H3K9ac")) +
    geom_line(aes(x = 1:(length(ts1))*0.001, y = ts2, color = "Transcriptome of Zheng et al. 2014")) +
    labs(x = "Time (h)", y = "Loading", title = "The Optimal alignment of curves") +
    theme_classic()
ggsave(file.path(dir_res, "Optimal_alignment_H3K9ac_dim2.pdf"), p, height = 5, width = 10)

# 绘制最小距离时的曲线 H3K18ac
ts1 <- interp_H3K18ac$value
ts2 <- interp_rnaseq$value
i <- timing_rnaseq_H3K18ac
shifted_ts1 <- if (i < 0) {
    c(rep(NA, abs(i)), ts1[1:(length(ts1) + i)])
} else {
    c(ts1[(i + 1):length(ts1)], rep(NA, i))
}
p <- ggplot() +
    geom_line(aes(x = 1:(length(ts1))*0.001, y = shifted_ts1, color = "Epigenome of histone modification H3K18ac")) +
    geom_line(aes(x = 1:(length(ts1))*0.001, y = ts2, color = "Transcriptome of Zheng et al. 2014")) +
    labs(x = "Time (h)", y = "Loading", title = "The Optimal alignment of curves") +
    theme_minimal()
ggsave(file.path(dir_res, "Optimal_alignment_H3K18ac_dim2.pdf"), p, height = 5, width = 10)


# subsampling -------------------------------------------------------------
# 
# timing_subsampling <- function(i){
#     
#     data_rnaseq_subsample <- data_rnaseq[sample(nrow(data_rnaseq), 13, replace = F), ]
#     data_H3K9ac_subsample <- data_H3K9ac[sample(nrow(data_H3K9ac), 13, replace = F), ]
#     data_H3K18ac_subsample <- data_H3K18ac[sample(nrow(data_H3K18ac), 13, replace = F), ]
#     
#     # 插值时间点
#     time_min_subsample <- max(min(data_rnaseq_subsample$Time), min(data_H3K9ac_subsample$Time), min(data_H3K18ac_subsample$Time), time_min)
#     time_max_subsample <- min(max(data_rnaseq_subsample$Time), max(data_H3K9ac_subsample$Time), max(data_H3K18ac_subsample$Time), time_max)
#     timepoints <- seq(from = time_min_subsample, to = time_max_subsample, by = 0.001)
#     
#     # Dim 1 
#     data_rnaseq_1 <- data_rnaseq_subsample %>% 
#         select(time = Time, value = `1`)
#     data_H3K9ac_1 <- data_H3K9ac_subsample %>% 
#         select(time = Time, value = `1`)
#     data_H3K18ac_1 <- data_H3K18ac_subsample %>% 
#         select(time = Time, value = `1`)
#     
#     # 数据标准化
#     if (is_norm == "Norm") {
#         data_rnaseq_1$value <- normalize(data_rnaseq_1$value)
#         data_H3K9ac_1$value <- normalize(data_H3K9ac_1$value)
#         data_H3K18ac_1$value <- normalize(data_H3K18ac_1$value)
#     }  else if (is_norm == "Scale") {
#         data_rnaseq_1$value <- scalize(data_rnaseq_1$value)
#         data_H3K9ac_1$value <- scalize(data_H3K9ac_1$value)
#         data_H3K18ac_1$value <- scalize(data_H3K18ac_1$value)
#     }
#     
#     # 插值计算
#     if (interpolation_method == "Linear") {
#         # linear
#         interp_rnaseq <- tibble(
#             time = timepoints,
#             value = coredata(na.approx(zoo(data_rnaseq_1$value, order.by = data_rnaseq_1$time), xout = timepoints))
#         )
#         interp_H3K9ac <- tibble(
#             time = timepoints,
#             value = coredata(na.approx(zoo(data_H3K9ac_1$value, order.by = data_H3K9ac_1$time), xout = timepoints))
#         )
#         interp_H3K18ac <- tibble(
#             time = timepoints,
#             value = coredata(na.approx(zoo(data_H3K18ac_1$value, order.by = data_H3K18ac_1$time), xout = timepoints))
#         )
#     } else if (interpolation_method == "Spline") {
#         # spline
#         interp_rnaseq <- tibble(
#             time = timepoints,
#             value = spline(data_rnaseq_1$time, data_rnaseq_1$value, xout = timepoints)$y 
#         )
#         interp_H3K9ac <- tibble( 
#             time = timepoints,
#             value = spline(data_H3K9ac_1$time, data_H3K9ac_1$value, xout = timepoints)$y
#         )
#         interp_H3K18ac <- tibble(
#             time = timepoints,
#             value = spline(data_H3K18ac_1$time, data_H3K18ac_1$value, xout = timepoints)$y
#         )
#     }
#     
#     # 计算最优对齐
#     if (align_method == "ccf") {
#         
#         # H3K9ac
#         # 计算ccf，并找到相位差
#         ccf_rnaseq_H3K9ac <- ccf(interp_H3K9ac$value, interp_rnaseq$value, lag.max = 1000, plot = F)
#         timing_rnaseq_H3K9ac <- ccf_rnaseq_H3K9ac$lag[which.max(ccf_rnaseq_H3K9ac$acf)]
#         
#         # H3K18ac
#         ccf_rnaseq_H3K18ac <- ccf(interp_H3K18ac$value, interp_rnaseq$value, lag.max = 1000, plot = F)
#         timing_rnaseq_H3K18ac <- ccf_rnaseq_H3K18ac$lag[which.max(ccf_rnaseq_H3K18ac$acf)]
#         
#     } else if (align_method == "distance1") {
#         
#         # H3K9ac
#         # 计算距离，并找到时间差
#         dis_rnaseq_H3K9ac <- distance_ts(interp_H3K9ac$value, interp_rnaseq$value, lag = 1000)
#         timing_rnaseq_H3K9ac <- dis_rnaseq_H3K9ac$Time_Difference[which.min(dis_rnaseq_H3K9ac$Distance)]
# 
#         # H3K18ac
#         dis_rnaseq_H3K18ac <- distance_ts(interp_H3K18ac$value, interp_rnaseq$value, lag = 1000)
#         timing_rnaseq_H3K18ac <- dis_rnaseq_H3K18ac$Time_Difference[which.min(dis_rnaseq_H3K18ac$Distance)]
# 
#     } else if (align_method == "distance2") {
#         
#         # H3K9ac
#         # 计算距离，并找到时间差
#         dis_rnaseq_H3K9ac <- distance2_ts(interp_H3K9ac$value, interp_rnaseq$value, lag = 1000)
#         timing_rnaseq_H3K9ac <- dis_rnaseq_H3K9ac$Time_Difference[which.min(dis_rnaseq_H3K9ac$Distance)]
# 
#         # H3K18ac
#         dis_rnaseq_H3K18ac <- distance2_ts(interp_H3K18ac$value, interp_rnaseq$value, lag = 1000)
#         timing_rnaseq_H3K18ac <- dis_rnaseq_H3K18ac$Time_Difference[which.min(dis_rnaseq_H3K18ac$Distance)]
# 
#     }
#     res <- c(
#         "H3K9ac_dim1" = timing_rnaseq_H3K9ac*0.001*60,
#         "H3K18ac_dim1" = timing_rnaseq_H3K18ac*0.001*60
#     )
#     
#     
#     # Dim 2 
#     data_rnaseq_1 <- data_rnaseq_subsample %>% 
#         select(time = Time, value = `2`)
#     data_H3K9ac_1 <- data_H3K9ac_subsample %>% 
#         select(time = Time, value = `2`)
#     data_H3K18ac_1 <- data_H3K18ac_subsample %>% 
#         select(time = Time, value = `2`)
#     
#     # 数据标准化
#     if (is_norm == "Norm") {
#         data_rnaseq_1$value <- normalize(data_rnaseq_1$value)
#         data_H3K9ac_1$value <- normalize(data_H3K9ac_1$value)
#         data_H3K18ac_1$value <- normalize(data_H3K18ac_1$value)
#     }  else if (is_norm == "Scale") {
#         data_rnaseq_1$value <- scalize(data_rnaseq_1$value)
#         data_H3K9ac_1$value <- scalize(data_H3K9ac_1$value)
#         data_H3K18ac_1$value <- scalize(data_H3K18ac_1$value)
#     }
# 
#     # 插值计算
#     if (interpolation_method == "Linear") {
#         # linear
#         interp_rnaseq <- tibble(
#             time = timepoints,
#             value = coredata(na.approx(zoo(data_rnaseq_1$value, order.by = data_rnaseq_1$time), xout = timepoints))
#         )
#         interp_H3K9ac <- tibble(
#             time = timepoints,
#             value = coredata(na.approx(zoo(data_H3K9ac_1$value, order.by = data_H3K9ac_1$time), xout = timepoints))
#         )
#         interp_H3K18ac <- tibble(
#             time = timepoints,
#             value = coredata(na.approx(zoo(data_H3K18ac_1$value, order.by = data_H3K18ac_1$time), xout = timepoints))
#         )
#     } else if (interpolation_method == "Spline") {
#         # spline
#         interp_rnaseq <- tibble(
#             time = timepoints,
#             value = spline(data_rnaseq_1$time, data_rnaseq_1$value, xout = timepoints)$y 
#         )
#         interp_H3K9ac <- tibble( 
#             time = timepoints,
#             value = spline(data_H3K9ac_1$time, data_H3K9ac_1$value, xout = timepoints)$y
#         )
#         interp_H3K18ac <- tibble(
#             time = timepoints,
#             value = spline(data_H3K18ac_1$time, data_H3K18ac_1$value, xout = timepoints)$y
#         )
#     }
#     
#     # 计算最优对齐
#     if (align_method == "ccf") {
#         
#         # H3K9ac
#         # 计算ccf，并找到相位差
#         ccf_rnaseq_H3K9ac <- ccf(interp_H3K9ac$value, interp_rnaseq$value, lag.max = 1000, plot = F)
#         timing_rnaseq_H3K9ac <- ccf_rnaseq_H3K9ac$lag[which.max(ccf_rnaseq_H3K9ac$acf)]
#         
#         # H3K18ac
#         ccf_rnaseq_H3K18ac <- ccf(interp_H3K18ac$value, interp_rnaseq$value, lag.max = 1000, plot = F)
#         timing_rnaseq_H3K18ac <- ccf_rnaseq_H3K18ac$lag[which.max(ccf_rnaseq_H3K18ac$acf)]
#         
#     } else if (align_method == "distance1") {
#         
#         # H3K9ac
#         # 计算距离，并找到时间差
#         dis_rnaseq_H3K9ac <- distance_ts(interp_H3K9ac$value, interp_rnaseq$value, lag = 1000)
#         timing_rnaseq_H3K9ac <- dis_rnaseq_H3K9ac$Time_Difference[which.min(dis_rnaseq_H3K9ac$Distance)]
#         
#         # H3K18ac
#         dis_rnaseq_H3K18ac <- distance_ts(interp_H3K18ac$value, interp_rnaseq$value, lag = 1000)
#         timing_rnaseq_H3K18ac <- dis_rnaseq_H3K18ac$Time_Difference[which.min(dis_rnaseq_H3K18ac$Distance)]
#         
#     } else if (align_method == "distance2") {
#         
#         # H3K9ac
#         # 计算距离，并找到时间差
#         dis_rnaseq_H3K9ac <- distance2_ts(interp_H3K9ac$value, interp_rnaseq$value, lag = 1000)
#         timing_rnaseq_H3K9ac <- dis_rnaseq_H3K9ac$Time_Difference[which.min(dis_rnaseq_H3K9ac$Distance)]
#         
#         # H3K18ac
#         dis_rnaseq_H3K18ac <- distance2_ts(interp_H3K18ac$value, interp_rnaseq$value, lag = 1000)
#         timing_rnaseq_H3K18ac <- dis_rnaseq_H3K18ac$Time_Difference[which.min(dis_rnaseq_H3K18ac$Distance)]
#         
#     }
# 
#     res <- c(
#         res,
#         "H3K9ac_dim2" = timing_rnaseq_H3K9ac*0.001*60,
#         "H3K18ac_dim2" = timing_rnaseq_H3K18ac*0.001*60
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
# # H3K9ac_dim1
# x <- timing_distribution[, "H3K9ac_dim1"]
# # 计算95%置信区间
# ci_95 <- quantile(x, probs = c(0.025, 0.975))  # 2.5百分位和97.5百分位
# # 绘制密度图
# pdf(file.path(dir_res, "Density_subsample_H3K9ac_dim1.pdf"), height = 6, width = 10)
# plot(density(x), main = str_glue("Density of timing; the 95% CI lower and upper: {ci_95[1]} and {ci_95[2]}"), xlab = "Timing (minutes)")
# # 标记95%置信区间的上下界
# abline(v = ci_95[1], col = "red", lty = 2, lwd = 2)  # 下界
# abline(v = ci_95[2], col = "red", lty = 2, lwd = 2)  # 上界
# text(ci_95[1], 0.1, str_glue("95% CI Lower: {ci_95[1]}"), pos = 2, col = "red")  # 下界注释
# text(ci_95[2], 0.1, str_glue("95% CI Upper: {ci_95[2]}"), pos = 4, col = "red")  # 上界注释
# dev.off()
# 
# # H3K9ac_dim2
# x <- timing_distribution[, "H3K9ac_dim2"]
# # 计算95%置信区间
# ci_95 <- quantile(x, probs = c(0.025, 0.975))  # 2.5百分位和97.5百分位
# # 绘制密度图
# pdf(file.path(dir_res, "Density_subsample_H3K9ac_dim2.pdf"), height = 6, width = 10)
# plot(density(x), main = str_glue("Density of timing; the 95% CI lower and upper: {ci_95[1]} and {ci_95[2]}"), xlab = "Timing (minutes)")
# # 标记95%置信区间的上下界
# abline(v = ci_95[1], col = "red", lty = 2, lwd = 2)  # 下界
# abline(v = ci_95[2], col = "red", lty = 2, lwd = 2)  # 上界
# text(ci_95[1], 0.1, str_glue("95% CI Lower: {ci_95[1]}"), pos = 2, col = "red")  # 下界注释
# text(ci_95[2], 0.1, str_glue("95% CI Upper: {ci_95[2]}"), pos = 4, col = "red")  # 上界注释
# dev.off()
# 
# # H3K18ac_dim1
# x <- timing_distribution[, "H3K18ac_dim1"]
# # 计算95%置信区间
# ci_95 <- quantile(x, probs = c(0.025, 0.975))  # 2.5百分位和97.5百分位
# # 绘制密度图
# pdf(file.path(dir_res, "Density_subsample_H3K18ac_dim1.pdf"), height = 6, width = 10)
# plot(density(x), main = str_glue("Density of timing; the 95% CI lower and upper: {ci_95[1]} and {ci_95[2]}"), xlab = "Timing (minutes)")
# # 标记95%置信区间的上下界
# abline(v = ci_95[1], col = "red", lty = 2, lwd = 2)  # 下界
# abline(v = ci_95[2], col = "red", lty = 2, lwd = 2)  # 上界
# text(ci_95[1], 0.1, str_glue("95% CI Lower: {ci_95[1]}"), pos = 2, col = "red")  # 下界注释
# text(ci_95[2], 0.1, str_glue("95% CI Upper: {ci_95[2]}"), pos = 4, col = "red")  # 上界注释
# dev.off()
# 
# # H3K9ac_dim2
# x <- timing_distribution[, "H3K18ac_dim2"]
# # 计算95%置信区间
# ci_95 <- quantile(x, probs = c(0.025, 0.975))  # 2.5百分位和97.5百分位
# # 绘制密度图
# pdf(file.path(dir_res, "Density_subsample_H3K18ac_dim2.pdf"), height = 6, width = 10)
# plot(density(x), main = str_glue("Density of timing; the 95% CI lower and upper: {ci_95[1]} and {ci_95[2]}"), xlab = "Timing (minutes)")
# # 标记95%置信区间的上下界
# abline(v = ci_95[1], col = "red", lty = 2, lwd = 2)  # 下界
# abline(v = ci_95[2], col = "red", lty = 2, lwd = 2)  # 上界
# text(ci_95[1], 0.1, str_glue("95% CI Lower: {ci_95[1]}"), pos = 2, col = "red")  # 下界注释
# text(ci_95[2], 0.1, str_glue("95% CI Upper: {ci_95[2]}"), pos = 4, col = "red")  # 上界注释
# dev.off()


# perturbation ------------------------------------------------------------

timing_perturbation <- function(i){

    # 插值时间点
    timepoints <- seq(from = time_min, to = time_max, by = 0.001)

    # Dim 1 
    data_rnaseq_1 <- data_rnaseq %>% 
        select(time = Time, value = `1`)
    data_H3K9ac_1 <- data_H3K9ac %>% 
        select(time = Time, value = `1`)
    data_H3K18ac_1 <- data_H3K18ac %>% 
        select(time = Time, value = `1`)
    
    # 添加随机扰动
    # 正态分布扰动，扰动标准差为0.1标准差
    data_rnaseq_1$value <- data_rnaseq_1$value + 
        rnorm(nrow(data_rnaseq_1), mean = 0, sd = 0.1 * sd(data_rnaseq_1$value))  
    data_H3K9ac_1$value <- data_H3K9ac_1$value + 
        rnorm(nrow(data_H3K9ac_1), mean = 0, sd = 0.1 * sd(data_H3K9ac_1$value))  
    data_H3K18ac_1$value <- data_H3K18ac_1$value + 
        rnorm(nrow(data_H3K18ac_1), mean = 0, sd = 0.1 * sd(data_H3K18ac_1$value))  

    # 插值时间点
    timepoints <- seq(from = time_min, to = time_max, by = 0.001)
    
    # 插值计算
    if (interpolation_method == "Linear") {
        # linear
        interp_rnaseq <- tibble(
            time = timepoints,
            value = coredata(na.approx(zoo(data_rnaseq_1$value, order.by = data_rnaseq_1$time), xout = timepoints))
        )
        interp_H3K9ac <- tibble(
            time = timepoints,
            value = coredata(na.approx(zoo(data_H3K9ac_1$value, order.by = data_H3K9ac_1$time), xout = timepoints))
        )
        interp_H3K18ac <- tibble(
            time = timepoints,
            value = coredata(na.approx(zoo(data_H3K18ac_1$value, order.by = data_H3K18ac_1$time), xout = timepoints))
        )
    } else if (interpolation_method == "Spline") {
        # spline
        interp_rnaseq <- tibble(
            time = timepoints,
            value = spline(data_rnaseq_1$time, data_rnaseq_1$value, xout = timepoints)$y 
        )
        interp_H3K9ac <- tibble( 
            time = timepoints,
            value = spline(data_H3K9ac_1$time, data_H3K9ac_1$value, xout = timepoints)$y
        )
        interp_H3K18ac <- tibble(
            time = timepoints,
            value = spline(data_H3K18ac_1$time, data_H3K18ac_1$value, xout = timepoints)$y
        )
    }
    
    # 数据标准化
    if (is_norm == "Norm") {
        interp_rnaseq$value <- normalize(interp_rnaseq$value)
        interp_H3K9ac$value <- normalize(interp_H3K9ac$value)
        interp_H3K18ac$value <- normalize(interp_H3K18ac$value)
    }  else if (is_norm == "Scale") {
        interp_rnaseq$value <- normalize(interp_rnaseq$value)
        interp_H3K9ac$value <- normalize(interp_H3K9ac$value)
        interp_H3K18ac$value <- normalize(interp_H3K18ac$value)
    }
    
    # 计算最优对齐
    if (align_method == "ccf") {
        
        # H3K9ac
        # 计算ccf，并找到相位差
        ccf_rnaseq_H3K9ac <- ccf(interp_H3K9ac$value, interp_rnaseq$value, lag.max = 1000, plot = F)
        timing_rnaseq_H3K9ac <- ccf_rnaseq_H3K9ac$lag[which.max(ccf_rnaseq_H3K9ac$acf)]
        
        # H3K18ac
        ccf_rnaseq_H3K18ac <- ccf(interp_H3K18ac$value, interp_rnaseq$value, lag.max = 1000, plot = F)
        timing_rnaseq_H3K18ac <- ccf_rnaseq_H3K18ac$lag[which.max(ccf_rnaseq_H3K18ac$acf)]
        
    } else if (align_method == "distance1") {
        
        # H3K9ac
        # 计算距离，并找到时间差
        dis_rnaseq_H3K9ac <- distance_ts(interp_H3K9ac$value, interp_rnaseq$value, lag = 1000)
        timing_rnaseq_H3K9ac <- dis_rnaseq_H3K9ac$Time_Difference[which.min(dis_rnaseq_H3K9ac$Distance)]
        
        # H3K18ac
        dis_rnaseq_H3K18ac <- distance_ts(interp_H3K18ac$value, interp_rnaseq$value, lag = 1000)
        timing_rnaseq_H3K18ac <- dis_rnaseq_H3K18ac$Time_Difference[which.min(dis_rnaseq_H3K18ac$Distance)]
        
    } else if (align_method == "distance2") {
        
        # H3K9ac
        # 计算距离，并找到时间差
        dis_rnaseq_H3K9ac <- distance2_ts(interp_H3K9ac$value, interp_rnaseq$value, lag = 1000)
        timing_rnaseq_H3K9ac <- dis_rnaseq_H3K9ac$Time_Difference[which.min(dis_rnaseq_H3K9ac$Distance)]
        
        # H3K18ac
        dis_rnaseq_H3K18ac <- distance2_ts(interp_H3K18ac$value, interp_rnaseq$value, lag = 1000)
        timing_rnaseq_H3K18ac <- dis_rnaseq_H3K18ac$Time_Difference[which.min(dis_rnaseq_H3K18ac$Distance)]
        
    }
    res <- c(
        "H3K9ac_dim1" = timing_rnaseq_H3K9ac*0.001*60,
        "H3K18ac_dim1" = timing_rnaseq_H3K18ac*0.001*60
    )
    
    
    # Dim 2 
    data_rnaseq_1 <- data_rnaseq %>% 
        select(time = Time, value = `2`)
    data_H3K9ac_1 <- data_H3K9ac %>% 
        select(time = Time, value = `2`)
    data_H3K18ac_1 <- data_H3K18ac %>% 
        select(time = Time, value = `2`)
    
    # 添加随机扰动
    # 正态分布扰动，扰动标准差为0.1标准差
    data_rnaseq_1$value <- data_rnaseq_1$value + 
        rnorm(nrow(data_rnaseq_1), mean = 0, sd = 0.1 * sd(data_rnaseq_1$value))  
    data_H3K9ac_1$value <- data_H3K9ac_1$value + 
        rnorm(nrow(data_H3K9ac_1), mean = 0, sd = 0.1 * sd(data_H3K9ac_1$value))  
    data_H3K18ac_1$value <- data_H3K18ac_1$value + 
        rnorm(nrow(data_H3K18ac_1), mean = 0, sd = 0.1 * sd(data_H3K18ac_1$value))  
    
    # 插值计算
    if (interpolation_method == "Linear") {
        # linear
        interp_rnaseq <- tibble(
            time = timepoints,
            value = coredata(na.approx(zoo(data_rnaseq_1$value, order.by = data_rnaseq_1$time), xout = timepoints))
        )
        interp_H3K9ac <- tibble(
            time = timepoints,
            value = coredata(na.approx(zoo(data_H3K9ac_1$value, order.by = data_H3K9ac_1$time), xout = timepoints))
        )
        interp_H3K18ac <- tibble(
            time = timepoints,
            value = coredata(na.approx(zoo(data_H3K18ac_1$value, order.by = data_H3K18ac_1$time), xout = timepoints))
        )
    } else if (interpolation_method == "Spline") {
        # spline
        interp_rnaseq <- tibble(
            time = timepoints,
            value = spline(data_rnaseq_1$time, data_rnaseq_1$value, xout = timepoints)$y 
        )
        interp_H3K9ac <- tibble( 
            time = timepoints,
            value = spline(data_H3K9ac_1$time, data_H3K9ac_1$value, xout = timepoints)$y
        )
        interp_H3K18ac <- tibble(
            time = timepoints,
            value = spline(data_H3K18ac_1$time, data_H3K18ac_1$value, xout = timepoints)$y
        )
    }
    
    # 数据标准化
    if (is_norm == "Norm") {
        interp_rnaseq$value <- normalize(interp_rnaseq$value)
        interp_H3K9ac$value <- normalize(interp_H3K9ac$value)
        interp_H3K18ac$value <- normalize(interp_H3K18ac$value)
    }  else if (is_norm == "Scale") {
        interp_rnaseq$value <- normalize(interp_rnaseq$value)
        interp_H3K9ac$value <- normalize(interp_H3K9ac$value)
        interp_H3K18ac$value <- normalize(interp_H3K18ac$value)
    }
    
    # 计算最优对齐
    if (align_method == "ccf") {
        
        # H3K9ac
        # 计算ccf，并找到相位差
        ccf_rnaseq_H3K9ac <- ccf(interp_H3K9ac$value, interp_rnaseq$value, lag.max = 1000, plot = F)
        timing_rnaseq_H3K9ac <- ccf_rnaseq_H3K9ac$lag[which.max(ccf_rnaseq_H3K9ac$acf)]
        
        # H3K18ac
        ccf_rnaseq_H3K18ac <- ccf(interp_H3K18ac$value, interp_rnaseq$value, lag.max = 1000, plot = F)
        timing_rnaseq_H3K18ac <- ccf_rnaseq_H3K18ac$lag[which.max(ccf_rnaseq_H3K18ac$acf)]
        
    } else if (align_method == "distance1") {
        
        # H3K9ac
        # 计算距离，并找到时间差
        dis_rnaseq_H3K9ac <- distance_ts(interp_H3K9ac$value, interp_rnaseq$value, lag = 1000)
        timing_rnaseq_H3K9ac <- dis_rnaseq_H3K9ac$Time_Difference[which.min(dis_rnaseq_H3K9ac$Distance)]
        
        # H3K18ac
        dis_rnaseq_H3K18ac <- distance_ts(interp_H3K18ac$value, interp_rnaseq$value, lag = 1000)
        timing_rnaseq_H3K18ac <- dis_rnaseq_H3K18ac$Time_Difference[which.min(dis_rnaseq_H3K18ac$Distance)]
        
    } else if (align_method == "distance2") {
        
        # H3K9ac
        # 计算距离，并找到时间差
        dis_rnaseq_H3K9ac <- distance2_ts(interp_H3K9ac$value, interp_rnaseq$value, lag = 1000)
        timing_rnaseq_H3K9ac <- dis_rnaseq_H3K9ac$Time_Difference[which.min(dis_rnaseq_H3K9ac$Distance)]
        
        # H3K18ac
        dis_rnaseq_H3K18ac <- distance2_ts(interp_H3K18ac$value, interp_rnaseq$value, lag = 1000)
        timing_rnaseq_H3K18ac <- dis_rnaseq_H3K18ac$Time_Difference[which.min(dis_rnaseq_H3K18ac$Distance)]
        
    }
    
    res <- c(
        res,
        "H3K9ac_dim2" = timing_rnaseq_H3K9ac*0.001*60,
        "H3K18ac_dim2" = timing_rnaseq_H3K18ac*0.001*60
    )

    
    return(res)
}

# a <- t(sapply(1:100, timing_subsampling))
sfInit(parallel = TRUE, cpus = CPUs)
sfExportAll()
sfLibrary(tidyverse)
sfLibrary(zoo)
timing_distribution <- sfSapply(1:10000, timing_perturbation) %>% t()
sfStop()


# 置信区间

# H3K9ac_dim1
x <- timing_distribution[, "H3K9ac_dim1"]
# 计算95%置信区间
ci_95 <- quantile(x, probs = c(0.025, 0.975))  # 2.5百分位和97.5百分位
# 绘制密度图
# pdf(file.path(dir_res, "Density_perturb_H3K9ac_dim1.pdf"), height = 6, width = 10)
svg(file.path(dir_res, "Density_perturb_H3K9ac_dim1.svg"), height = 6, width = 10)
plot(density(x), main = str_glue("Density of time difference; the 95% CI lower and upper: {ci_95[1]} and {ci_95[2]}"), xlab = "Time difference (minutes)")
# 标记95%置信区间的上下界
abline(v = ci_95[1], col = "red", lty = 2, lwd = 2)  # 下界
abline(v = ci_95[2], col = "red", lty = 2, lwd = 2)  # 上界
text(ci_95[1], 0.1, str_glue("95% CI Lower"), pos = 2, col = "red")  # 下界注释
text(ci_95[2], 0.1, str_glue("95% CI Upper"), pos = 4, col = "red")  # 上界注释
dev.off()

# H3K9ac_dim2
x <- timing_distribution[, "H3K9ac_dim2"]
# 计算95%置信区间
ci_95 <- quantile(x, probs = c(0.025, 0.975))  # 2.5百分位和97.5百分位
# 绘制密度图
# pdf(file.path(dir_res, "Density_perturb_H3K9ac_dim2.pdf"), height = 6, width = 10)
svg(file.path(dir_res, "Density_perturb_H3K9ac_dim2.svg"), height = 6, width = 10)
plot(density(x), main = str_glue("Density of time difference; the 95% CI lower and upper: {ci_95[1]} and {ci_95[2]}"), xlab = "Time difference (minutes)")
# 标记95%置信区间的上下界
abline(v = ci_95[1], col = "red", lty = 2, lwd = 2)  # 下界
abline(v = ci_95[2], col = "red", lty = 2, lwd = 2)  # 上界
text(ci_95[1], 0.1, str_glue("95% CI Lower"), pos = 2, col = "red")  # 下界注释
text(ci_95[2], 0.1, str_glue("95% CI Upper"), pos = 4, col = "red")  # 上界注释
dev.off()

# H3K18ac_dim1
x <- timing_distribution[, "H3K18ac_dim1"]
# 计算95%置信区间
ci_95 <- quantile(x, probs = c(0.025, 0.975))  # 2.5百分位和97.5百分位
# 绘制密度图
# pdf(file.path(dir_res, "Density_perturb_H3K18ac_dim1.pdf"), height = 6, width = 10)
svg(file.path(dir_res, "Density_perturb_H3K18ac_dim1.svg"), height = 6, width = 10)
plot(density(x), main = str_glue("Density of time difference; the 95% CI lower and upper: {ci_95[1]} and {ci_95[2]}"), xlab = "Time difference (minutes)")
# 标记95%置信区间的上下界
abline(v = ci_95[1], col = "red", lty = 2, lwd = 2)  # 下界
abline(v = ci_95[2], col = "red", lty = 2, lwd = 2)  # 上界
text(ci_95[1], 0.1, str_glue("95% CI Lower"), pos = 2, col = "red")  # 下界注释
text(ci_95[2], 0.1, str_glue("95% CI Upper"), pos = 4, col = "red")  # 上界注释
dev.off()

# H3K9ac_dim2
x <- timing_distribution[, "H3K18ac_dim2"]
# 计算95%置信区间
ci_95 <- quantile(x, probs = c(0.025, 0.975))  # 2.5百分位和97.5百分位
# 绘制密度图
# pdf(file.path(dir_res, "Density_perturb_H3K18ac_dim2.pdf"), height = 6, width = 10)
svg(file.path(dir_res, "Density_perturb_H3K18ac_dim2.svg"), height = 6, width = 10)
plot(density(x), main = str_glue("Density of time difference; the 95% CI lower and upper: {ci_95[1]} and {ci_95[2]}"), xlab = "Time difference (minutes)")
# 标记95%置信区间的上下界
abline(v = ci_95[1], col = "red", lty = 2, lwd = 2)  # 下界
abline(v = ci_95[2], col = "red", lty = 2, lwd = 2)  # 上界
text(ci_95[1], 0.1, str_glue("95% CI Lower"), pos = 2, col = "red")  # 下界注释
text(ci_95[2], 0.1, str_glue("95% CI Upper"), pos = 4, col = "red")  # 上界注释
dev.off()
