

replaceH <- function(vec) {
  vec <- stringr::str_replace(vec, "H", vec[1])
  return(vec)
}

#' 生成GAPIT的geno文件格式
#'
#'@param geno 数据为Soybase网站下载的snp50k_Williams_82.a24数据，数据列为CHROM POS ID ...
#'@param threshhold 删除稀有变异预值，默认为0.05。计算时不包含H和缺失
#'@return 返回用于gapit 计算的标准数据，其中材料已按名称升序排列
#'

gapit_geno_format <-
  function(geno = snp50k_Williams_82.a24,  threshold = 0.05) {
    ##去除Scaffold
    geno <-
      geno[stringr::str_detect(geno$CHROM, "Gm"), ]

    factor_cols <- sapply(geno, is.factor)

    # 将所有因子变量转换为字符串
    geno[factor_cols] <- lapply(geno[factor_cols], as.character)
    # g2，将H去除，用于统计ATGC个数,并生成alleles列
    g2 <- geno[names(geno)[-1:-3]]
    g2 <- g2[order(colnames(g2))]
    ## 不去除H，用于生成最后文件
    g3 <- g2


    g2[g2 == "U"] <- NA
    g2[g2 == "H"] <- NA

    g2_table <- apply(g2, 1, table)
    g2_table_names <- lapply(g2_table, names)


    alleles <- unlist(lapply(g2_table_names, paste, collapse = "/"))

    rs <- paste(geno$CHROM, geno$POS, geno$ID, sep = "_")
    chrom <- as.numeric(stringr::str_remove(geno$CHROM, "Gm"))

    info_head <- data.frame(rs, alleles, chrom)
    info_head$pos <- geno$POS
    info_head$strand <- NA
    info_head$assembly <- NA
    info_head$center <- NA
    info_head$protLSID <- NA
    info_head$assayLSID <- NA
    info_head$panel <- NA
    info_head$QCcode <- NA


    #处理H
    g3[g3 == "A"] <- "AA"
    g3[g3 == "T"] <- "TT"
    g3[g3 == "G"] <- "GG"
    g3[g3 == "C"] <- "CC"
    g3[g3 == "U"] <- NA
    g3[is.na(g3)] <- "NN"

    hg3 <- data.frame(alleles = alleles, g3)

    myapply <- apply(hg3, 1, replaceH)
    myapply <- t(myapply)
    colnames(myapply) <- colnames(hg3)
    myapply <- data.frame(myapply)
    myapply <- myapply[, -1]

    standard_geno <- cbind(info_head, myapply)


    #去除稀有变异，<-10%
    minsnp <- unlist(lapply(g2_table, min))
    maxsnp <- unlist(lapply(g2_table, max))
    minsnp[minsnp == maxsnp] <- 0

    minorsele <- minsnp / (minsnp + maxsnp)
    standard_geno_sele <- standard_geno[minorsele > threshold, ]


    return(standard_geno_sele)
  }
