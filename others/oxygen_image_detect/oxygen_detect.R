library(tidyverse)


# LC-MS -------------------------------------------------------------------

data <- read_csv('LC-MS_all.csv') %>%
    filter(x >= 50, y <= 270 )

data_res <- data %>%
    mutate(x = (x - 55)/30 + 40,
           y = max(y) - y,
           y = (53/max(y)) * y + 20) %>%
    group_by(x) %>%
    summarise(y_mean = mean(y)) %>%
    rename(Time = x, dO2 = y_mean)

write_csv(data_res, 'LC-MS_oxygen.csv')




# TOFMS -------------------------------------------------------------------

data <- read_csv('TOFMS_all.csv') %>%
    filter(x >= 50, y <= 270 )

data_res <- data %>%
    mutate(x = (x - 51)/30 + 64,
           y = max(y) - y,
           y = (36/max(y)) * y + 22) %>%
    group_by(x) %>%
    summarise(y_mean = mean(y)) %>%
    rename(Time = x, dO2 = y_mean)

write_csv(data_res, 'TOFMS_oxygen.csv')


# YMC2014 -------------------------------------------------------------------

data <- read_csv('YMC2014_oxygen_rnaseq_all.csv') 

data_res <- data %>%
    mutate(x = ((x - 63)*4.3)/1594 + 47.95,
           y = max(y) - y,
           y = (42/max(y)) * y + 11) %>%
    group_by(x) %>%
    summarise(y_mean = mean(y)) %>%
    rename(Time = x, dO2 = y_mean)

write_csv(data_res, 'YMC2014_oxygen_rnaseq.csv')



data <- read_csv('YMC2014_oxygen_chipseq_all.csv') 


data_res <- data %>%
    mutate(x = ((x - 13)*4.46)/1657 + 49.99,
           y = max(y) - y,
           y = (48/max(y)) * y + 6) %>%
    group_by(x) %>%
    summarise(y_mean = mean(y)) %>%
    rename(Time = x, dO2 = y_mean)

write_csv(data_res, 'YMC2014_oxygen_chipseq.csv')



# YMC2005 -----------------------------------------------------------------


data <- read_csv('YMC2005_oxygen_all.csv') 

ggplot(data, mapping = aes(x = x, y = -y)) +
    geom_point()


data_res <- data %>%
    mutate(x = (x*35)/1024 + 1,
           y = max(y) - y,
           y = (45/max(y)) * y + 5) %>%
    group_by(x) %>%
    summarise(y_mean = mean(y)) %>%
    rename(Time = x, dO2 = y_mean)

ggplot(data_res, mapping = aes(x = Time, y = dO2)) +
    geom_point()


write_csv(data_res, 'YMC2005_oxygen.csv')

