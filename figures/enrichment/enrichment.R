suppressMessages({
    library(tidyverse)
    library(openxlsx)
})


data <- read.xlsx("../../results/concatenate/enrichment/yeast_all.xlsx") %>%
    as_tibble() %>% 
    select(-description)

pathway_data_1 <- read.xlsx("pathway_ids.xlsx", sheet = 1) %>%
    as_tibble()  %>%
    left_join(data, by = "id") %>% 
    mutate(
        pathwayGeneNum = as.integer(pathwayGeneNum),
        pvalue = pvalue.gene1,
        `-log10P` = -sign(pvalue) * log10(abs(pvalue)),
        database = database %>%
            case_match(
                "go.bp" ~ "GO.bp",
                "go.cc" ~ "GO.cc",
                "go.mf" ~ "GO.mf",
                "kegg" ~ "KEGG",
                "reactome" ~ "Reactome"
            ) %>%
            factor(
                levels = c("GO.bp", "GO.cc", "GO.mf", "KEGG", "Reactome")
            ),
        pathway = str_c(description, database, sep = ", ") %>%
            fct_inorder() %>%
            fct_rev()
    ) %>%
    select(
        id,
        pathway,
        pathwayGeneNum,
        `-log10P`,
    )

pathway_data_2 <- read.xlsx("pathway_ids.xlsx", sheet = 2) %>%
    as_tibble()  %>%
    left_join(data, by = "id") %>% 
    mutate(
        pathwayGeneNum = as.integer(pathwayGeneNum),
        pvalue = pvalue.gene2,
        `-log10P` = -sign(pvalue) * log10(abs(pvalue)),
        database = database %>%
            case_match(
                "go.bp" ~ "GO.bp",
                "go.cc" ~ "GO.cc",
                "go.mf" ~ "GO.mf",
                "kegg" ~ "KEGG",
                "reactome" ~ "Reactome"
            ) %>%
            factor(
                levels = c("GO.bp", "GO.cc", "GO.mf", "KEGG", "Reactome")
            ),
        pathway = str_c(description, database, sep = ", ") %>%
            fct_inorder() %>%
            fct_rev()
    ) %>%
    select(
        id,
        pathway,
        pathwayGeneNum,
        `-log10P`,
    )


# plot --------------------------------------------------------------------

# set p-value color gradient
colorPalette.pos <- c("white",  "orange2")
colorList.pos <- colorRampPalette(colorPalette.pos)(20L)[c(10L, 17L, 18L, 19L, 20L, 20L, 20L, 20L, 20L, 20L)]
colorPalette.neg <- c("white", "royalblue")
colorList.neg <- colorRampPalette(colorPalette.neg)(10L)[c(4L, 5L, 6L, 7L, 8L, 9L, 10L, 10L, 10L)]
colorList <- c(rev(colorList.neg), "white", colorList.pos)
# scales::show_col(colorList)

# # 500*600
# # 100=2.645cm
# ggplot(
#     data = pathway_data_1,
#     mapping = aes(
#         x = "yeast 1st",
#         y = pathway,
#         color = `-log10P`,
#         size = pathwayGeneNum
#     )
# ) +
#     geom_point() +
#     labs(
#         x = NULL,
#         y = NULL,
#         color = expression(-Log[10] * P),
#         size = "Gene set size"
#     ) +
#     scale_x_discrete(
#         labels = ~ str_wrap(., width = 5L)
#     ) +
#     scale_y_discrete(
#         labels = ~ str_wrap(., width = 55L)
#     ) +
#     scale_color_gradientn(
#         colors = colorList,
#         limits = c(-95, 15),
#         values = c(0, 95/110, 1),
#         breaks = c(-90, -60, -30, 0, 10),
#         labels = c(90, 60, 30, 0, 10) %>% as.character()
#     ) +
#     theme_classic() +
#     theme(
#         legend.title = element_text(size = 7),
#         legend.text = element_text(size = 6),
#         legend.key.size = unit(0.2, "inches"),
#         axis.text.x = element_blank(),
#         axis.line.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.y = element_text(
#             size = 10,
#             color = "black"
#         ),
#         axis.line.y = element_line(
#             linetype = 1,
#             color = "black",
#             linewidth = 0.3
#         ),
#         axis.ticks.y = element_line(
#             color = "black",
#             linewidth = 0.3,
#             lineend = 1
#         )
#     ) 

ggplot(
    data = pathway_data_1,
    mapping = aes(
        x = "yeast 1st",
        y = pathway,
        color = `-log10P`,
        size = pathwayGeneNum
    )
    ) +
    geom_point() +
    labs(
        x = NULL,
        y = NULL,
        color = expression(-Log[10] * P),
        size = "Gene set size"
    ) +
    scale_size(range = c(1, 3)) +
    scale_x_discrete(labels = ~ str_wrap(., width = 5L)) +
    scale_y_discrete(labels = ~ str_wrap(., width = 55L)) +
    scale_color_gradientn(
        colors = colorList,
        limits = c(-15, 95),
        values = c(0, 15/110, 1),
        breaks = c(-10, 0, 10, 30, 60, 90),
        labels = c(10, 0, 10, 30, 60, 90) %>% as.character()
    ) +
    theme_classic() +
    theme(
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "inches"),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.line.y = element_line(linetype = 1, color = "black", linewidth = 0.25),
        axis.ticks.y = element_line(color = "black", linewidth = 0.25, lineend = 1),
    ) 
ggsave("Enrichment1.pdf", width = 9.5, height = 9.5, units = "cm" ) 



# set p-value color gradient
colorPalette.pos <- c("white", "#FF4848")
colorList.pos <- colorRampPalette(colorPalette.pos)(20L)[c(14L, 15L, 16L, 17L, 18L, 19L, 20L, 20L, 20L, 20L)]
colorPalette.neg <- c("white", "#36A336")
colorList.neg <- colorRampPalette(colorPalette.neg)(20L)[c(16L, 17L, 18L, 19L, 20L, 20L, 20L, 20L, 20L, 20L)]
colorList <- c(rev(colorList.neg), "white", colorList.pos)
# scales::show_col(colorList)

# # 500*600
# ggplot(
#     data = pathway_data_2,
#     mapping = aes(
#         x = "yeast 2rd",
#         y = pathway,
#         color = `-log10P`,
#         size = pathwayGeneNum
#     )
# ) +
#     geom_point() +
#     labs(
#         x = NULL,
#         y = NULL,
#         color = expression(-Log[10] * P),
#         size = "Gene set size"
#     ) +
#     scale_x_discrete(
#         labels = ~ str_wrap(., width = 5L)
#     ) +
#     scale_y_discrete(
#         labels = ~ str_wrap(., width = 55L)
#     ) +
#     scale_color_gradientn(
#         colors = colorList,
#         limits = c(-55, 85),
#         values = c(0, 55/140, 1),
#         breaks = c(-50, -20, 0, 20, 50, 80),
#         labels = c(50, 20, 0, 20, 50, 80) %>% as.character()
#     ) +
#     theme_classic() +
#     theme(
#         legend.title = element_text(size = 7),
#         legend.text = element_text(size = 6),
#         legend.key.size = unit(0.2, "inches"),
#         axis.text.x = element_blank(),
#         axis.line.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.y = element_text(
#             size = 10,
#             color = "black"
#         ),
#         axis.line.y = element_line(
#             linetype = 1,
#             color = "black",
#             linewidth = 0.3
#         ),
#         axis.ticks.y = element_line(
#             color = "black",
#             linewidth = 0.3,
#             lineend = 1
#         )
#     )

ggplot(
    data = pathway_data_2,
    mapping = aes(
        x = "yeast 2rd",
        y = pathway,
        color = `-log10P`,
        size = pathwayGeneNum
    )
    ) +
    geom_point() +
    labs(
        x = NULL,
        y = NULL,
        color = expression(-Log[10] * P),
        size = "Gene set size"
    ) +
    scale_size(range = c(1, 3)) +
    scale_x_discrete(labels = ~ str_wrap(., width = 5L)) +
    scale_y_discrete(labels = ~ str_wrap(., width = 55L)) +
    scale_color_gradientn(
        colors = colorList,
        limits = c(-55, 85),
        values = c(0, 55/140, 1),
        breaks = c(-50, -20, 0, 20, 50, 80),
        labels = c(50, 20, 0, 20, 50, 80) %>% as.character()
    ) +
    theme_classic() +
    theme(
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "inches"),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.line.y = element_line(linetype = 1, color = "black", linewidth = 0.25),
        axis.ticks.y = element_line(color = "black", linewidth = 0.25, lineend = 1),
    ) 
ggsave("Enrichment2.pdf", width = 9.5, height = 9.5, units = "cm" ) 
