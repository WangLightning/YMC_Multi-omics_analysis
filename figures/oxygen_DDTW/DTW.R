if (!require("dtw")) {
    install.packages("dtw")
}
library(dtw)
library(DTWBI)
library(tidyverse)



# read data
data_raw_YMC2005 <- read_csv("YMC2005_oxygen.csv", show_col_types = FALSE) %>% 
    mutate(
        n = row_number(),
        Time_old = Time,
        dO2_old = dO2,
        Time_dtw = Time_old * 5/12 ,
        dO2_dtw = (dO2_old - min(dO2_old)) / (max(dO2_old) - min(dO2_old)),
        Time = Time_dtw,
        dO2 = dO2_dtw*40 + 20
    )
data_raw_YMC2014 <- read_csv("YMC2014_oxygen_rnaseq.csv", show_col_types = FALSE) %>% 
    mutate(
        n = row_number(),
        Time_old = Time,
        dO2_old = dO2,
        Time_dtw = Time_old - min(Time_old),
        dO2_dtw = (dO2_old - min(dO2_old)) / (max(dO2_old) - min(dO2_old)),
        Time = Time_dtw,
        dO2 = dO2_dtw*40 + 20
    )
data_raw_histone <- read_csv("YMC2014_oxygen_chipseq.csv", show_col_types = FALSE) %>% 
    mutate(
        n = row_number(),
        Time_dtw = Time - min(Time),
        dO2_dtw = (dO2 - min(dO2)) / (max(dO2) - min(dO2))
    )
data_raw_meta_lc <- read_csv("LC-MS_oxygen.csv", show_col_types = FALSE) %>% 
    filter(Time > 42.7) %>% 
    mutate(
        n = row_number(),
        Time_dtw = Time - min(Time),
        dO2_dtw = (dO2 - min(dO2)) / (max(dO2) - min(dO2))
    )
data_raw_meta_gc <- read_csv("TOFMS_oxygen.csv", show_col_types = FALSE) %>% 
    filter(Time > 65.5) %>% 
    mutate(
        n = row_number(),
        Time_dtw = Time - min(Time),
        dO2_dtw = (dO2 - min(dO2)) / (max(dO2) - min(dO2))
    )

times_YMC2005 <- c(
    1,
    30,
    60,
    89,
    118,
    147,
    177,
    206,
    235,
    264,
    294,
    323,
    352,
    381,
    411,
    440,
    469,
    498,
    528,
    557,
    586,
    615,
    645,
    674,
    703,
    732,
    762,
    791,
    820,
    849,
    879,
    908,
    937,
    966,
    996,
    1025
)


times_meta_lc <- c(
    4,
    16,
    28,
    37,
    46,
    56,
    68,
    79,
    90,
    100,
    110,
    121,
    127,
    137,
    146,
    156,
    167,
    176,
    184,
    194,
    205,
    216,
    227,
    238
)


times_meta_gc <- c(
    2,
    14,
    24,
    36,
    48,
    60,
    72,
    83,
    95,
    107,
    118,
    127,
    136,
    147,
    157,
    167,
    178,
    190,
    202,
    212,
    224,
    236,
    248,
    255
)

# times_YMC2014 <- c(
#     138,
#     253,
#     298,
#     327,
#     383,
#     446,
#     509,
#     568,
#     631,
#     694,
#     754,
#     820,
#     883,
#     1084,
#     1280,
#     1484
# )
# 由于原坐标轴与时间有所偏差，所以对这里的时间进行了微调
times_YMC2014 <- c(
    133,
    248,
    291,
    322,
    378,
    441,
    504,
    563,
    626,
    689,
    754,
    815,
    878,
    1079,
    1275,
    1479
)
times_histone <- c(
    191,
    339,
    406,
    443,
    503,
    566,
    629,
    688,
    752,
    815,
    889,
    948,
    1108,
    1275,
    1428,
    1584
)

data_point_YMC2005 <- data_raw_YMC2005 %>% 
    filter(n %in% times_YMC2005) %>% 
    mutate(t = 1:36)
data_point_YMC2014 <- data_raw_YMC2014 %>% 
    filter(n %in% times_YMC2014) %>% 
    mutate(t = 1:16)


# 保存ref的氧气曲线
data_ref_05 <- data_raw_YMC2005 %>% 
    select(Time, dO2)
write_csv(data_ref_05, "oxygen_ref_YMC2005.csv")

data_ref_14 <- data_raw_YMC2014 %>% 
    select(Time, dO2)
write_csv(data_ref_14, "oxygen_ref_YMC2014.csv")



# lc to YMC2005 -----------------------------------------------------------

# ggplot() +
#     geom_line(data = data_raw_meta_lc, aes(x = Time_dtw, y = dO2_dtw, colour = 'LC-MS'), linewidth = 1) +
#     geom_line(data = data_raw_YMC2005, aes(x = Time_dtw, y = dO2_dtw, colour = 'Tu et al. 2005'), linewidth = 1)


# DTW
# dtw_result <- dtw(data_raw_meta_lc$dO2_dtw, data_raw_YMC2005$dO2_dtw, step.pattern = symmetric1)
# dtw_result <- dtw(data_raw_meta_lc$dO2_dtw, data_raw_YMC2005$dO2_dtw, step.pattern = symmetric2)
# dtw_result <- dtw(data_raw_meta_lc$dO2_dtw, data_raw_YMC2005$dO2_dtw, step.pattern = asymmetric, open.end = TRUE, open.begin = TRUE)
# dtw_result <- dtw(data_raw_meta_lc$dO2_dtw, data_raw_YMC2005$dO2_dtw, step.pattern = asymmetricP0, open.end = TRUE, open.begin = TRUE)
# dtw_result <- dtw(data_raw_meta_lc$dO2_dtw, data_raw_YMC2005$dO2_dtw, step.pattern = asymmetricP05, open.end = TRUE, open.begin = TRUE)
dtw_result <- dtw(local.derivative.ddtw(data_raw_meta_lc$dO2_dtw), local.derivative.ddtw(data_raw_YMC2005$dO2_dtw), step.pattern = asymmetricP05, open.end = TRUE, open.begin = TRUE)
# dtw_result <- dtw(data_raw_meta_lc$dO2_dtw, data_raw_YMC2005$dO2_dtw, step.pattern = asymmetricP1, open.end = TRUE, open.begin = TRUE)
# dtw_result <- dtw(data_raw_meta_lc$dO2_dtw, data_raw_YMC2005$dO2_dtw, step.pattern = asymmetricP2, open.end = TRUE, open.begin = TRUE)


data_mapping_plot1 <- tibble(
    t1 = data_raw_meta_lc$Time_dtw[dtw_result$index1]+0.5,
    y1 = data_raw_meta_lc$dO2_dtw[dtw_result$index1]+1,
    t2 = data_raw_YMC2005$Time_dtw[dtw_result$index2],
    y2 = data_raw_YMC2005$dO2_dtw[dtw_result$index2]
)
data_mapping_plot1 <- data_mapping_plot1[seq(1, nrow(data_mapping_plot1), by = 5),]
p <- ggplot() +
    geom_line(data = data_raw_meta_lc, aes(x = Time_dtw+0.5, y = dO2_dtw+1, colour = 'Metabolome of LC-MS'), linewidth = 1) +
    geom_line(data = data_raw_YMC2005, aes(x = Time_dtw, y = dO2_dtw, colour = 'Transcriptome of Tu et al. 2005'), linewidth = 1) +
    geom_segment(data = data_mapping_plot1, aes(x = t1, y = y1, xend = t2, yend = y2), colour = "gray") +
    labs(title = "DDTW alignment results", x = "Time (h)", y = "dO2", colour = "data") +
    xlim(c(0,12)) +
    theme_classic() +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "top"
    )
# print(p)
# ggsave("DTW_lc_YMC2005_asymmetricP05_plot1.pdf", p, height = 5, width = 10)
ggsave("DDTW_lc_YMC2005_asymmetricP05_plot1.pdf", p, height = 5, width = 7)


data_mapping_plot2 <- tibble(
    Time = data_raw_YMC2005$Time_dtw[dtw_result$index2],
    dO2 = data_raw_meta_lc$dO2_dtw[dtw_result$index1]+1,
    n = data_raw_meta_lc$n[dtw_result$index1]
    ) %>%
    distinct(n, .keep_all = T)
# p <- ggplot() +
#     geom_line(data = data_mapping_plot2, aes(x = Time, y = dO2, colour = 'Metabolome LC-MS'), linewidth = 1) +
#     geom_line(data = data_raw_YMC2005, aes(x = Time_dtw, y = dO2_dtw, colour = 'Tu et al. 2005'), linewidth = 1) +
#     labs(title = "DTW alignment between LC-MS and Tu et al. 2005", x = "Time (h)", y = "dO2", colour = "data") +
#     theme_light() +
#     theme(
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#     )
# print(p)
# ggsave("DTW_lc_YMC2005_asymmetricP05_plot2.pdf", p, height = 5, width = 10)
# ggsave("DDTW_lc_YMC2005_asymmetricP05_plot2.pdf", p, height = 5, width = 10)


data_point_lc <- data_mapping_plot2 %>% 
    filter(n %in% times_meta_lc) %>% 
    mutate(t = 1:24)
p <- ggplot() +
    geom_line(data = data_mapping_plot2, aes(x = Time, y = dO2, colour = 'Metabolome of LC-MS'), linewidth = 1) +
    geom_line(data = data_raw_YMC2005, aes(x = Time_dtw, y = dO2_dtw, colour = 'Transcriptome of Tu et al. 2005'), linewidth = 1) +
    geom_point(data = data_point_lc, aes(x = Time, y = dO2), color = "grey30") +
    # geom_text(data = data_point_lc, aes(x = Time, y = dO2, label = t), vjust = 1.5, hjust = -0.2, size = 3) +
    geom_point(data = data_point_YMC2005, aes(x = Time, y = dO2_dtw), color = "grey30") +
    # geom_text(data = data_point_YMC2005, aes(x = Time, y = dO2_dtw, label = t), vjust = 1.5, hjust = -0.2, size = 3) +
    labs(title = "The aligned oxygen concentration curves", x = "Time (h)", y = "dO2", colour = "data") +
    xlim(c(0,12)) +
    theme_light() +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "top"
    )
# print(p)
# ggsave("DTW_lc_YMC2005_asymmetricP05_plot3.pdf", p, height = 5, width = 10)
ggsave("DDTW_lc_YMC2005_asymmetricP05_plot3.pdf", p, height = 5, width = 7)


# gc to YMC2005 -----------------------------------------------------------

# ggplot() +
#     geom_line(data = data_raw_meta_gc, aes(x = Time_dtw, y = dO2_dtw, colour = 'GC-TFOMS'), linewidth = 1) +
#     geom_line(data = data_raw_YMC2005, aes(x = Time_dtw, y = dO2_dtw, colour = 'Tu et al. 2005'), linewidth = 1)


# DTW
# dtw_result <- dtw(data_raw_meta_gc$dO2_dtw, data_raw_YMC2005$dO2_dtw, step.pattern = symmetric1)
# dtw_result <- dtw(data_raw_meta_gc$dO2_dtw, data_raw_YMC2005$dO2_dtw, step.pattern = symmetric2)
# dtw_result <- dtw(data_raw_meta_gc$dO2_dtw, data_raw_YMC2005$dO2_dtw, step.pattern = asymmetric, open.end = TRUE, open.begin = TRUE)
# dtw_result <- dtw(data_raw_meta_gc$dO2_dtw, data_raw_YMC2005$dO2_dtw, step.pattern = asymmetricP0, open.end = TRUE, open.begin = TRUE)
# dtw_result <- dtw(data_raw_meta_gc$dO2_dtw, data_raw_YMC2005$dO2_dtw, step.pattern = asymmetricP05, open.end = TRUE, open.begin = TRUE)
dtw_result <- dtw(local.derivative.ddtw(data_raw_meta_gc$dO2_dtw), local.derivative.ddtw(data_raw_YMC2005$dO2_dtw), step.pattern = asymmetricP05, open.end = TRUE, open.begin = TRUE)
# dtw_result <- dtw(data_raw_meta_gc$dO2_dtw, data_raw_YMC2005$dO2_dtw, step.pattern = asymmetricP1, open.end = TRUE, open.begin = TRUE)
# dtw_result <- dtw(data_raw_meta_gc$dO2_dtw, data_raw_YMC2005$dO2_dtw, step.pattern = asymmetricP2, open.end = TRUE, open.begin = TRUE)

data_mapping_plot1 <- tibble(
    t1 = data_raw_meta_gc$Time_dtw[dtw_result$index1]+0.5,
    y1 = data_raw_meta_gc$dO2_dtw[dtw_result$index1]+1,
    t2 = data_raw_YMC2005$Time_dtw[dtw_result$index2],
    y2 = data_raw_YMC2005$dO2_dtw[dtw_result$index2]
)
data_mapping_plot1 <- data_mapping_plot1[seq(1, nrow(data_mapping_plot1), by = 5),]
p <- ggplot() +
    geom_line(data = data_raw_meta_gc, aes(x = Time_dtw+0.5, y = dO2_dtw+1, colour = 'Metabolome of GC-TFOMS'), linewidth = 1) +
    geom_line(data = data_raw_YMC2005, aes(x = Time_dtw, y = dO2_dtw, colour = 'Transcriptome of Tu et al. 2005'), linewidth = 1) +
    geom_segment(data = data_mapping_plot1, aes(x = t1, y = y1, xend = t2, yend = y2), colour = "gray") +
    labs(title = "DDTW alignment results", x = "Time (h)", y = "dO2", colour = "data") +
    xlim(c(0,12)) +
    theme_classic() +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "top"
    )
# print(p)
# ggsave("DTW_gc_YMC2005_asymmetricP05_plot1.pdf", p, height = 5, width = 10)
ggsave("DDTW_gc_YMC2005_asymmetricP05_plot1.pdf", p, height = 5, width = 7)


data_mapping_plot2 <- tibble(
    Time = data_raw_YMC2005$Time_dtw[dtw_result$index2],
    dO2 = data_raw_meta_gc$dO2_dtw[dtw_result$index1]+1,
    n = data_raw_meta_gc$n[dtw_result$index1]
) %>% 
    distinct(n, .keep_all = T)
# p <- ggplot() +
#     geom_line(data = data_mapping_plot2, aes(x = Time, y = dO2, colour = 'Metabolome GC-TFOMS'), linewidth = 1) +
#     geom_line(data = data_raw_YMC2005, aes(x = Time_dtw, y = dO2_dtw, colour = 'Tu et al. 2005'), linewidth = 1) +
#     labs(title = "DTW alignment between GC-TFOMS and Tu et al. 2005", x = "Time (h)", y = "dO2", colour = "data") +
#     theme_light() +
#     theme(
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#     )
# # print(p)
# # ggsave("DTW_gc_YMC2005_asymmetricP05_plot2.pdf", p, height = 5, width = 10)
# ggsave("DDTW_gc_YMC2005_asymmetricP05_plot2.pdf", p, height = 5, width = 10)


data_point_gc <- data_mapping_plot2 %>% 
    filter(n %in% times_meta_gc) %>% 
    mutate(t = 1:24)
p <- ggplot() +
    geom_line(data = data_mapping_plot2, aes(x = Time, y = dO2, colour = 'Metabolome of GC-TFOMS'), linewidth = 1) +
    geom_line(data = data_raw_YMC2005, aes(x = Time_dtw, y = dO2_dtw, colour = 'Transcriptome of Tu et al. 2005'), linewidth = 1) +
    geom_point(data = data_point_gc, aes(x = Time, y = dO2), color = "grey30") +
    # geom_text(data = data_point_gc, aes(x = Time, y = dO2, label = t), vjust = 1.5, hjust = -0.2, size = 3) +
    geom_point(data = data_point_YMC2005, aes(x = Time, y = dO2_dtw), color = "grey30") +
    # geom_text(data = data_point_YMC2005, aes(x = Time, y = dO2_dtw, label = t), vjust = 1.5, hjust = -0.2, size = 3) +
    labs(title = "The aligned oxygen concentration curves", x = "Time (h)", y = "dO2", colour = "data") +
    xlim(c(0,12)) +
    theme_light() +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "top"
    )
# print(p)
# ggsave("DTW_gc_YMC2005_asymmetricP05_plot3.pdf", p, height = 5, width = 10)
ggsave("DDTW_gc_YMC2005_asymmetricP05_plot3.pdf", p, height = 5, width = 7)


# histone to YMC2014 -----------------------------------------------------------

# ggplot() +
#     geom_line(data = data_raw_histone, aes(x = Time_dtw, y = dO2_dtw, colour = 'Histone modification'), linewidth = 1) +
#     geom_line(data = data_raw_YMC2014, aes(x = Time_dtw, y = dO2_dtw, colour = 'Zheng et al. 2014'), linewidth = 1)


# DTW
# dtw_result <- dtw(data_raw_histone$dO2_dtw, data_raw_YMC2014$dO2_dtw, step.pattern = asymmetric, open.end = TRUE, open.begin = TRUE)
# dtw_result <- dtw(data_raw_histone$dO2_dtw, data_raw_YMC2014$dO2_dtw, step.pattern = asymmetricP05, open.end = TRUE, open.begin = TRUE)
# dtw_result <- dtw(data_raw_histone$dO2_dtw, data_raw_YMC2014$dO2_dtw, step.pattern = symmetric1)
# dtw_result <- dtw(data_raw_histone$dO2_dtw, data_raw_YMC2014$dO2_dtw, step.pattern = symmetric2)
# dtw_result <- dtw(data_raw_histone$dO2_dtw, data_raw_YMC2014$dO2_dtw, step.pattern = asymmetricP0, open.end = TRUE, open.begin = TRUE)
# dtw_result <- dtw(data_raw_histone$dO2_dtw, data_raw_YMC2014$dO2_dtw, step.pattern = asymmetricP1, open.end = TRUE, open.begin = TRUE)
# dtw_result <- dtw(data_raw_histone$dO2_dtw, data_raw_YMC2014$dO2_dtw, step.pattern = asymmetricP2, open.end = TRUE, open.begin = TRUE)
dtw_result <- dtw(local.derivative.ddtw(data_raw_histone$dO2_dtw), local.derivative.ddtw(data_raw_YMC2014$dO2_dtw), step.pattern = asymmetricP2, open.end = TRUE, open.begin = TRUE)
# dtw_result <- dtw(data_raw_histone$dO2_dtw, data_raw_YMC2014$dO2_dtw, step.pattern = asymmetricP2, window.type = "itakura", open.end = TRUE, open.begin = TRUE)
# dtw_result <- dtw(data_raw_histone$dO2_dtw, data_raw_YMC2014$dO2_dtw, step.pattern = asymmetricP2, window.type = "sakoechiba", window.size = 100, open.end = TRUE, open.begin = TRUE)


data_mapping_plot1 <- tibble(
    t1 = data_raw_histone$Time_dtw[dtw_result$index1],
    y1 = data_raw_histone$dO2_dtw[dtw_result$index1]+1,
    t2 = data_raw_YMC2014$Time_dtw[dtw_result$index2]+0.2,
    y2 = data_raw_YMC2014$dO2_dtw[dtw_result$index2]
)
data_mapping_plot1 <- data_mapping_plot1[seq(1, nrow(data_mapping_plot1), by = 8),]
p <- ggplot() +
    geom_line(data = data_raw_histone, aes(x = Time_dtw, y = dO2_dtw+1, colour = 'Epigenome of histone modification'), linewidth = 1) +
    geom_line(data = data_raw_YMC2014, aes(x = Time_dtw+0.2, y = dO2_dtw, colour = 'Transcriptome of Zheng et al. 2014'), linewidth = 1) +
    geom_segment(data = data_mapping_plot1, aes(x = t1, y = y1, xend = t2, yend = y2), colour = "gray") +
    labs(title = "DDTW alignment results", x = "Time (h)", y = "dO2", colour = "data") +
    theme_classic() +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "top"
    )
# print(p)
# ggsave("DTW_histone_YMC2014_asymmetricP2_plot1.pdf", p, height = 5, width = 10)
ggsave("DDTW_histone_YMC2014_asymmetricP2_plot1.pdf", p, height = 5, width = 7)


data_mapping_plot2 <- tibble(
    Time = data_raw_YMC2014$Time_dtw[dtw_result$index2],
    dO2 = data_raw_histone$dO2_dtw[dtw_result$index1]+1,
    n = data_raw_histone$n[dtw_result$index1]
    ) %>% 
    distinct(n, .keep_all = T)
# p <- ggplot() +
#     geom_line(data = data_mapping_plot2, aes(x = Time, y = dO2, colour = 'Histone modification'), linewidth = 1) +
#     geom_line(data = data_raw_YMC2014, aes(x = Time_dtw, y = dO2_dtw, colour = 'Zheng et al. 2014'), linewidth = 1) +
#     labs(title = "DTW alignment between Histone modification and Zheng et al. 2014", x = "Time (h)", y = "dO2", colour = "data") +
#     theme_light() +
#     theme(
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#     )
# # print(p)
# # ggsave("DTW_histone_YMC2014_asymmetricP2_plot2.pdf", p, height = 5, width = 10)
# ggsave("DDTW_histone_YMC2014_asymmetricP2_plot2.pdf", p, height = 5, width = 10)

data_point_histone <- data_mapping_plot2 %>% 
    filter(n %in% times_histone) %>% 
    mutate(t = 1:16)
p <- ggplot() +
    geom_line(data = data_mapping_plot2, aes(x = Time, y = dO2, colour = 'Epigenome of histone modification'), linewidth = 1) +
    geom_line(data = data_raw_YMC2014, aes(x = Time_dtw, y = dO2_dtw, colour = 'Transcriptome of Zheng et al. 2014'), linewidth = 1) +
    geom_point(data = data_point_histone, aes(x = Time, y = dO2), color = "grey30") +
    # geom_text(data = data_point_histone, aes(x = Time, y = dO2, label = t), vjust = 1.5, hjust = -0.2, size = 3) +
    geom_point(data = data_point_YMC2014, aes(x = Time, y = dO2_dtw), color = "grey30") +
    # geom_text(data = data_point_YMC2014, aes(x = Time, y = dO2_dtw, label = t), vjust = 1.5, hjust = -0.2, size = 3) +
    labs(title = "The aligned oxygen concentration curves", x = "Time (h)", y = "dO2", colour = "data") +
    theme_light() +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "top"
    )
# print(p)
# ggsave("DTW_histone_YMC2014_asymmetricP2_plot3.pdf", p, height = 5, width = 10)
ggsave("DDTW_histone_YMC2014_asymmetricP2_plot3.pdf", p, height = 5, width = 7)



