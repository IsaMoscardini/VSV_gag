# VSV_gag

rm(list = ls())
options(stringsAsFactors = F)

sapply(c('ggplot2','ggrepel','dplyr','Glimma', "tidyr", "edgeR", "tmod"), require,character.only=T)


#### Quality Control ####

descript_gag <- read.csv2("data/VSV_gag_desc.csv")
count_gag <- read.csv2("data/VSV_gag_COUNTS.csv")
View(count_gag)

rownames(descript_gag) <- paste0("X", descript_gag$sample)
rownames(count_gag) <- count_gag$Gene
count_gag$Gene <- NULL
identical(colnames(count_gag), rownames(descript_gag))

<