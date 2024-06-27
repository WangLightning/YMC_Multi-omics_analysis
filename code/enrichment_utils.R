
enrichment_wilcoxon_test.1 <- function(pathway.data, loading.data) {

    # 富集计算方法1：
    # 以如下两组基因的载荷作为样本做 Wilcoxon 两样本秩和检验：
    # （1）基因特征向量全部基因和某个通路的交；
    # （2）基因特征向量全部基因和“所有通路并集减该通路“的交；
    # 注意：每个通路正极的 p-value 和负极的 p-value 和为1。

    pathway.data <- pathway.data %>%
        select(
            id,
            geneNum,
            description,
            # all_of(transAPIList),
            ends_with("Trans"),
            gene
        ) %>%
        mutate(
            gene = gene %>% str_split(",")
        )

    all_pathway_gene <- pathway.data %>%
        pull(gene) %>%
        flatten_chr() %>%
        unique() %>%
        sort()
    loading.data <- loading.data %>%
        # filter(Symbol %in% all_pathway_gene)
        filter(toupper(Symbol) %in% toupper(all_pathway_gene))
    vector_gene_num <- nrow(loading.data)

    result <- pathway.data %>%
        rename(
            pathwayGeneNum = geneNum,
            pathwayGene = gene
        ) %>%
        rowwise() %>%
        mutate(
            # index = list(loading.data$Symbol %in% pathwayGene),
            index = list(toupper(loading.data$Symbol) %in% toupper(pathwayGene)),
            num.of.targets = sum(index),
            num.of.non.targets = vector_gene_num - num.of.targets,
            pathwayGene_in_geneEigenvectors = loading.data$Symbol[index] %>%
                # sort() %>%
                str_c(collapse = ",") %>%
                list(),
            pathwayGene = pathwayGene %>%
                str_c(collapse = ",") %>%
                list(),
            # 不可用if_else，ifelse短路但if_else不短路，
            # 在condition为FALSE的位置依然会计算true参数
            pvalue.pos = ifelse(
                num.of.targets >= 1L,
                wilcox.test(
                    loading.data$Loading[index],
                    loading.data$Loading[!index],
                    "greater"
                )$p.val,
                NA_real_
            ),
            pvalue.neg = ifelse(
                num.of.targets >= 1L,
                wilcox.test(
                    loading.data$Loading[index],
                    loading.data$Loading[!index],
                    "less"
                )$p.val,
                NA_real_
            )
        ) %>%
        ungroup() %>%
        select(all_of(outputCols))

    result

}

enrichment_wilcoxon_test.2 <- function(pathway.data, loading.data) {

    # 富集计算方法2：
    # 以如下两组基因的载荷作为样本做 Wilcoxon 两样本秩和检验：
    # （1）基因特征向量全部基因和某个通路的交；
    # （2）基因特征向量全部基因减该通路；
    # 注意：每个通路正极的 p-value 和负极的 p-value 和为1。

    pathway.data <- pathway.data %>%
        select(
            id,
            geneNum,
            description,
            # all_of(transAPIList),
            ends_with("Trans"),
            gene
        ) %>%
        mutate(
            gene = gene %>% str_split(",")
        )

    vector_gene_num <- nrow(loading.data)

    result <- pathway.data %>%
        rename(
            pathwayGeneNum = geneNum,
            pathwayGene = gene
        ) %>%
        rowwise() %>%
        mutate(
            # index = list(loading.data$Symbol %in% pathwayGene),
            index = list(toupper(loading.data$Symbol) %in% toupper(pathwayGene)),
            num.of.targets = sum(index),
            num.of.non.targets = vector_gene_num - num.of.targets,
            pathwayGene_in_geneEigenvectors = loading.data$Symbol[index] %>%
                # sort() %>%
                str_c(collapse = ",") %>%
                list(),
            pathwayGene = pathwayGene %>%
                str_c(collapse = ",") %>%
                list(),
            # 不可用if_else，ifelse短路但if_else不短路，
            # 在condition为FALSE的位置依然会计算true参数
            pvalue.pos = ifelse(
                num.of.targets >= 1L,
                wilcox.test(
                    loading.data$Loading[index],
                    loading.data$Loading[!index],
                    "greater"
                )$p.val,
                NA_real_
            ),
            pvalue.neg = ifelse(
                num.of.targets >= 1L,
                wilcox.test(
                    loading.data$Loading[index],
                    loading.data$Loading[!index],
                    "less"
                )$p.val,
                NA_real_
            )
        ) %>%
        ungroup() %>%
        select(all_of(outputCols))

    result

}

enrichment_wilcoxon_test.3 <- function(pathway.data, loading.data) {

    # 富集计算方法3:
    # 计算正极的 p-value 时，把基因特征向量中负的载荷都置为0；
    # 计算负极的 p-value 时，把基因特征向量中正的载荷都置为0；
    # 然后以这两个新的载荷向量作为参与检验的载荷向量，
    # 用方法1分别计算正极的 p—value 和负极的 p—value。
    # 注意：每个通路正极的 p-value 和负极的 p-value 和不为1。

    pathway.data <- pathway.data %>%
        select(
            id,
            geneNum,
            description,
            # all_of(transAPIList),
            ends_with("Trans"),
            gene
        ) %>%
        mutate(
            gene = gene %>% str_split(",")
        )

    all_pathway_gene <- pathway.data %>%
        pull(gene) %>%
        flatten_chr() %>%
        unique() %>%
        sort()
    loading.data <- loading.data %>%
        mutate(
            Loading.pos = if_else(Loading >= 0, Loading, 0),
            Loading.neg = if_else(Loading <= 0, Loading, 0)
        ) %>%
        # filter(Symbol %in% all_pathway_gene)
        filter(toupper(Symbol) %in% toupper(all_pathway_gene))
    vector_gene_num <- nrow(loading.data)

    result <- pathway.data %>%
        rename(
            pathwayGeneNum = geneNum,
            pathwayGene = gene
        ) %>%
        rowwise() %>%
        mutate(
            # index = list(loading.data$Symbol %in% pathwayGene),
            index = list(toupper(loading.data$Symbol) %in% toupper(pathwayGene)),
            num.of.targets = sum(index),
            num.of.non.targets = vector_gene_num - num.of.targets,
            pathwayGene_in_geneEigenvectors = loading.data$Symbol[index] %>%
                # sort() %>%
                str_c(collapse = ",") %>%
                list(),
            pathwayGene = pathwayGene %>%
                str_c(collapse = ",") %>%
                list(),
            # 不可用if_else，ifelse短路但if_else不短路，
            # 在condition为FALSE的位置依然会计算true参数
            pvalue.pos = ifelse(
                num.of.targets >= 1L,
                wilcox.test(
                    loading.data$Loading.pos[index],
                    loading.data$Loading.pos[!index],
                    "greater"
                )$p.val,
                NA_real_
            ),
            pvalue.neg = ifelse(
                num.of.targets >= 1L,
                wilcox.test(
                    loading.data$Loading.neg[index],
                    loading.data$Loading.neg[!index],
                    "less"
                )$p.val,
                NA_real_
            )
        ) %>%
        ungroup() %>%
        select(all_of(outputCols))

    result

}

enrichment_wilcoxon_test.4 <- function(pathway.data, loading.data) {

    # 富集计算方法4：
    # 计算正极的 p-value 时，把基因特征向量中负的载荷都置为0；
    # 计算负极的 p-value 时，把基因特征向量中正的载荷都置为0；
    # 然后以这两个新的载荷向量作为参与检验的载荷向量，
    # 用方法2分别计算正极的 p—value 和负极的 p—value。
    # 注意：每个通路正极的 p-value 和负极的 p-value 和不为1。

    pathway.data <- pathway.data %>%
        select(
            id,
            geneNum,
            description,
            # all_of(transAPIList),
            ends_with("Trans"),
            gene
        ) %>%
        mutate(
            gene = gene %>% str_split(",")
        )

    loading.data <- loading.data %>%
        mutate(
            Loading.pos = if_else(Loading >= 0, Loading, 0),
            Loading.neg = if_else(Loading <= 0, Loading, 0)
        )
    vector_gene_num <- nrow(loading.data)

    result <- pathway.data %>%
        rename(
            pathwayGeneNum = geneNum,
            pathwayGene = gene
        ) %>%
        rowwise() %>%
        mutate(
            # index = list(loading.data$Symbol %in% pathwayGene),
            index = list(toupper(loading.data$Symbol) %in% toupper(pathwayGene)),
            num.of.targets = sum(index),
            num.of.non.targets = vector_gene_num - num.of.targets,
            pathwayGene_in_geneEigenvectors = loading.data$Symbol[index] %>%
                # sort() %>%
                str_c(collapse = ",") %>%
                list(),
            pathwayGene = pathwayGene %>%
                str_c(collapse = ",") %>%
                list(),
            # 不可用if_else，ifelse短路但if_else不短路，
            # 在condition为FALSE的位置依然会计算true参数
            pvalue.pos = ifelse(
                num.of.targets >= 1L,
                wilcox.test(
                    loading.data$Loading.pos[index],
                    loading.data$Loading.pos[!index],
                    "greater"
                )$p.val,
                NA_real_
            ),
            pvalue.neg = ifelse(
                num.of.targets >= 1L,
                wilcox.test(
                    loading.data$Loading.neg[index],
                    loading.data$Loading.neg[!index],
                    "less"
                )$p.val,
                NA_real_
            )
        ) %>%
        ungroup() %>%
        select(all_of(outputCols))

    result

}

enrichment_wilcoxon_test <- function(pathway.data, loading.data, method) {
    fun <- str_glue("enrichment_wilcoxon_test.{method}")
    eval(call(fun, pathway.data, loading.data))
}


do_enrichment <- function(list, loading.dir) {

    s <- list$species
    f <- list$family
    g <- list$gene
    m <- 1

    loading.data <- read_tsv(
        file.path(
            loading.dir,
            str_glue("{g}_loadings.txt")
        ),
        col_types = "cd"
    )
    result <- list()
    for (p in pathway_list) {
        result[[p]] <- enrichment_wilcoxon_test(
            pathway.data = pathway_data.list[[f]][[p]],
            loading.data = loading.data,
            method = m
        )
    }

    result

}

write_raw_result <- function(enrichment_result, raw_result.dir) {

    # # 单文件存储
    # list.save(
    #     enrichment_result,
    #     file = file.path(raw_result.dir, "enrichment_result.RData")
    # )

    # 分物种存储
    for (s in species_list) {
        list.save(
            enrichment_result[[s]],
            file = file.path(
                raw_result.dir,
                str_glue("enrichment_result_{s}.RData")
            )
        )
    }

}

read_raw_result <- function(raw_result.dir) {

    # 单文件读取
    # enrichment_result <- list.load(
    #     file = file.path(raw_result.dir, "enrichment_result.RData")
    # )

    # 分物种读取
    enrichment_result <- list()
    for (s in species_list) {
        enrichment_result[[s]] <- list.load(
            file = file.path(
                raw_result.dir,
                str_glue("enrichment_result_{s}.RData")
            )
        )
    }

    enrichment_result

}

mapColNames2Width <- function(colNames) {
    if (colNames == "id") width <- 12L
    else if (colNames == "database") width <- 10L
    else if (colNames == "description") width <- 50L
    else if (colNames %>% str_ends("Trans")) width <- 40L
    else if (colNames %>% str_detect("pvalue")) width <- 13L
    else if (colNames == "pathwayGeneNum") width <- 16L
    else if (colNames == "num.of.targets") width <- 16L
    else if (colNames == "num.of.non.targets") width <- 16L
    else if (colNames == "pathwayGene")  width <- 100L
    else if (colNames == "pathwayGene_in_geneEigenvectors") width <- 100L
    width
}
mapColNames2Width <- Vectorize(mapColNames2Width, USE.NAMES = FALSE)

formatResult1 <- function(enrichment_result, s, g) {
    result <- enrichment_result[[s]][[g]]
    for (p in pathway_list) {
        result[[p]] <- result[[p]] %>%
            arrange(pvalue.pos) %>%
            filter(
                num.of.targets %>% between(minCommonGeneNum, maxCommonGeneNum)
            ) %>%
            select(all_of(outputCols))
    }
    result
}

saveResult1 <- function(result, s, g) {
    wb <- createWorkbook()
    for (p in pathway_list) {

        addWorksheet(wb, p)
        setColWidths(
            wb = wb,
            sheet = p,
            cols = seq_along(outputCols),
            widths = mapColNames2Width(outputCols)
        )
        writeData(wb, p, result[[p]])

        colorCols <- str_which(outputCols, "description|Trans")
        colorRows <- 1L + seq_len(nrow(result[[p]]))
        # 0 < pos <= 0.01
        conditionalFormatting(
            wb = wb,
            sheet = p,
            cols = colorCols,
            rows = colorRows,
            rule = "AND($B2>0,$B2<=0.01)",
            style = createStyle(bgFill = bgFill$pos[1]),
            type = "expression"
        )
        # 0.01 < pos <= 0.05
        conditionalFormatting(
            wb = wb,
            sheet = p,
            cols = colorCols,
            rows = colorRows,
            rule = "AND($B2>0.01,$B2<=0.05)",
            style = createStyle(bgFill = bgFill$pos[2]),
            type = "expression"
        )
        # 0.05 < pos <= 0.10
        conditionalFormatting(
            wb = wb,
            sheet = p,
            cols = colorCols,
            rows = colorRows,
            rule = "AND($B2>0.05,$B2<=0.10)",
            style = createStyle(bgFill = bgFill$pos[3]),
            type = "expression"
        )
        # 0 < neg <= 0.01
        conditionalFormatting(
            wb = wb,
            sheet = p,
            cols = colorCols,
            rows = colorRows,
            rule = "AND($C2>0,$C2<=0.01)",
            style = createStyle(bgFill = bgFill$neg[1]),
            type = "expression"
        )
        # 0.01 < neg <= 0.05
        conditionalFormatting(
            wb = wb,
            sheet = p,
            cols = colorCols,
            rows = colorRows,
            rule = "AND($C2>0.01,$C2<=0.05)",
            style = createStyle(bgFill = bgFill$neg[2]),
            type = "expression"
        )
        # 0.05 < neg <= 0.10
        conditionalFormatting(
            wb = wb,
            sheet = p,
            cols = colorCols,
            rows = colorRows,
            rule = "AND($C2>0.05,$C2<=0.10)",
            style = createStyle(bgFill = bgFill$neg[3]),
            type = "expression"
        )
    }

    saveWorkbook(
        wb = wb,
        file = file.path(
            enrichment.dir,
            str_glue("{s}_{g}.xlsx")
        ),
        overwrite = T
    )
}

formatResult2 <- function(enrichment_result, s, g) {
    result.temp <- enrichment_result[[s]][[g]]
    result <- list()
    for (pole in c("pos", "neg")) {
        pvalue.pole <- sym(str_c("pvalue.", pole, collapse = ""))
        for (p in pathway_list) {
            result[[pole]][[p]] <- result.temp[[p]] %>%
                filter(!!pvalue.pole <= signif_level) %>%
                arrange(!!pvalue.pole) %>%
                add_column(database = p, .before = 1)
        }
        result[[pole]] <- result[[pole]] %>%
            bind_rows() %>%
            rename(pvalue = !!pvalue.pole) %>%
            filter(
                num.of.targets %>% between(minCommonGeneNum, maxCommonGeneNum)
            ) %>%
            select(all_of(outputCols))
    }
    result
}

saveResult2 <- function(result, s, g) {
    wb.pos_neg <- createWorkbook()
    for (pole in c("pos", "neg")) {
        addWorksheet(wb.pos_neg, pole)
        setColWidths(
            wb = wb.pos_neg,
            sheet = pole,
            cols = seq_along(outputCols),
            widths = mapColNames2Width(outputCols)
        )
        writeData(wb.pos_neg, pole, result[[pole]])
    }
    saveWorkbook(
        wb = wb.pos_neg,
        file = file.path(
            enrichment.dir,
            str_glue("{s}_{g}_pos_neg.xlsx")
        ),
        overwrite = T
    )
}

formatResult3 <- function(enrichment_result, s) {
    result.temp <- enrichment_result[[s]]
    for (g in gene_list) {
        for (p in pathway_list) {
            result.temp[[g]][[p]] <- result.temp[[g]][[p]] %>%
                arrange(id) %>%
                add_column(database = p, .before = 1)
        }
        result.temp[[g]] <- result.temp[[g]] %>%
            bind_rows()
    }
    # 检查是否有通路的正极和负极的p值都不高于 signif_level
    # print(s)
    # check_pvalue <- result.temp %>%
    #     list.mapv(min(pmax(pvalue.pos, pvalue.neg), na.rm = TRUE))
    # print(check_pvalue)
    # if (min(check_pvalue) <= signif_level)
    #     cat(str_glue("有通路的正负极p值都不高于{signif_level}!\n"))
    # check_pvalue <- result.temp %>%
    #     list.mapv(sum(pmax(pvalue.pos, pvalue.neg) <= signif_level, na.rm = TRUE)) %>%
    #     unname()
    # print(check_pvalue)
    # if (min(check_pvalue) >= 1)
    #     cat(str_glue("有通路的正负极p值都不高于{signif_level}!\n"))
    result <- result.temp[[gene_list[1]]] %>%
        select(!c(pvalue.pos, pvalue.neg))
    for (g in gene_list) {
        result <- result %>%
            mutate(
                !!str_glue("pvalue.{g}") :=
                    # 假设每个通路只在一极显著（绝大多数情况），
                    # 因此只考虑p值更小的那极
                    with(
                        result.temp[[g]],
                        if_else(
                            pvalue.pos <= pvalue.neg,
                            pvalue.pos,
                            -pvalue.neg
                        )
                    ),
                .before = description
            )
    }
    pvalue_columns <- syms(str_c("pvalue.", gene_list))
    result <- result %>%
        # 将不显著的p-value替换为缺失值
        modify_at(
            .at = str_c("pvalue.", gene_list),
            .f = ~ if_else(between(., -signif_level, signif_level), ., NA_real_)
        ) %>%
        # 过滤掉在每一维都不显著的通路
        filter(!is.na(coalesce(!!!pvalue_columns)))
    result <- result %>%
        # 对 pathwayGene_in_geneEigenvectors 列按基因首字母排序
        mutate(
            pathwayGene_in_geneEigenvectors = pathwayGene_in_geneEigenvectors %>%
                str_split(",") %>%
                map(., sort) %>%
                map(., str_c, collapse = ",")
        )
    result <- result %>%
        filter(
            num.of.targets %>% between(minCommonGeneNum, maxCommonGeneNum)
        ) %>%
        select(all_of(outputCols))
    result
}

saveResult3 <- function(result, s) {
    wb <- createWorkbook()
    addWorksheet(wb, "Sheet 1")
    setColWidths(
        wb = wb,
        sheet = "Sheet 1",
        cols = seq_along(outputCols),
        widths = mapColNames2Width(outputCols)
    )
    writeData(wb, "Sheet 1", result)

    colorCols <- str_which(outputCols, "pvalue.gene")
    colorRows <- 1L + seq_len(nrow(result))
    # 0 < pos <= 0.01
    conditionalFormatting(
        wb = wb,
        sheet = "Sheet 1",
        cols = colorCols,
        rows = colorRows,
        rule = "AND(C2>0,C2<=0.01)",
        style = createStyle(bgFill = bgFill$pos[1]),
        type = "expression"
    )
    # 0.01 < pos <= 0.05
    conditionalFormatting(
        wb = wb,
        sheet = "Sheet 1",
        cols = colorCols,
        rows = colorRows,
        rule = "AND(C2>0.01,C2<=0.05)",
        style = createStyle(bgFill = bgFill$pos[2]),
        type = "expression"
    )
    # 0.05 < pos <= 0.10
    conditionalFormatting(
        wb = wb,
        sheet = "Sheet 1",
        cols = colorCols,
        rows = colorRows,
        rule = "AND(C2>0.05,C2<=0.10)",
        style = createStyle(bgFill = bgFill$pos[3]),
        type = "expression"
    )
    # 0 < neg <= 0.01
    conditionalFormatting(
        wb = wb,
        sheet = "Sheet 1",
        cols = colorCols,
        rows = colorRows,
        rule = "AND(C2<0,C2>=-0.01)",
        style = createStyle(bgFill = bgFill$neg[1]),
        type = "expression"
    )
    # 0.01 < neg <= 0.05
    conditionalFormatting(
        wb = wb,
        sheet = "Sheet 1",
        cols = colorCols,
        rows = colorRows,
        rule = "AND(C2<-0.01,C2>=-0.05)",
        style = createStyle(bgFill = bgFill$neg[2]),
        type = "expression"
    )
    # 0.05 < neg <= 0.10
    conditionalFormatting(
        wb = wb,
        sheet = "Sheet 1",
        cols = colorCols,
        rows = colorRows,
        rule = "AND(C2<-0.05,C2>=-0.10)",
        style = createStyle(bgFill = bgFill$neg[3]),
        type = "expression"
    )
    saveWorkbook(
        wb = wb,
        file = file.path(
            enrichment.dir,
            str_glue("{s}_all.xlsx")
        ),
        overwrite = T
    )
}

compIntersection <- function(species_x, gene_x, species_y, gene_y, pathway, idSet_list) {
    x <- chuck(idSet_list, species_x, gene_x, pathway)
    y <- chuck(idSet_list, species_y, gene_y, pathway)
    length(intersect(x, y)) / min(length(x), length(y))
}

compGoSemanticSimilariy <- function(species_x, gene_x, species_y, gene_y, pathway, idSet_list) {
    go1 <- chuck(idSet_list, species_x, gene_x, pathway)
    go2 <- chuck(idSet_list, species_y, gene_y, pathway)
    d <- GOSemSim::godata(ont = str_sub(pathway, start = 4L) %>% toupper())
    GOSemSim::mgoSim(go1, go2, semData = d, measure = "Wang", combine = "BMA")
}

do_compGOSemSim <- function(s1_s2_g) {
    map(
        pathway_list,
        ~ compGoSemanticSimilariy(
            species_x = s1_s2_g$species_x,
            species_y = s1_s2_g$species_y,
            gene_x = s1_s2_g$gene,
            gene_y = s1_s2_g$gene,
            pathway = .,
            idSet_list = idSet_list
        )
    )
}
