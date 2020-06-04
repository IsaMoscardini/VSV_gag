###### 
rm(list = ls())
options(stringsAsFactors = F)

sapply(c('ggplot2','ggrepel','dplyr','Glimma', "tidyr", "edgeR", "tmod"), require,character.only=T)


#### EdgeR using Placebo samples ####

descript_gag <- read.csv2("data/VSV_gag_desc.csv")
count_gag <- read.csv2("data/VSV_gag_COUNTS.csv")

rownames(descript_gag) <- paste0("X", descript_gag$sample)
rownames(count_gag) <- count_gag$Gene
count_gag$Gene <- NULL

descript_gag <- descript_gag[descript_gag$Visit == 8,]
View(descript_gag)
common <- intersect(colnames(count_gag), rownames(descript_gag))
count_gag <- count_gag[, common]

identical(colnames(count_gag), rownames(descript_gag))

# create DEGlist
y_ebovac <- DGEList(counts = count_gag, genes = row.names(count_gag), group= descript_gag$Dose)

keep <- rowSums(cpm(y_ebovac)>1)>=10
y_ebovac <- y_ebovac[keep,]
dim(y_ebovac) # 13124    52

# normalize
y.1 <- calcNormFactors(y_ebovac)

# design matrix
batch <- factor(descript_gag$batch)
descript_gag$Dose <- relevel(as.factor(descript_gag$Dose), ref = "P")
design <- model.matrix(~batch+descript_gag$Dose)
y.2 <- estimateDisp(y.1, design)

# Fitting general linear models to each feature
qlfit_new <- glmQLFit(y.2, design)
colnames(design)

####
####
#### Find DEGs Day 7
# T1
DEA_D7_T1 <- glmTreat(qlfit_new, coef =4, lfc = log2(1.2))
summary(decideTestsDGE(DEA_D7_T1))
# descript_gag$DoseT1
# Down                     0
# NotSig               13123
# Up                       1
toptags07 <- topTags(DEA_D7_T1, n = Inf)
View(toptags07$table)
write.csv2(toptags07$table, "Results/EdgeR/Placebo/DEA_D7_T1", row.names = FALSE)

# T2
DEA_D7_T2 <- glmTreat(qlfit_new, coef =5, lfc = log2(1.2))
summary(decideTestsDGE(DEA_D7_T2))
# descript_gag$DoseT2
# Down                     0
# NotSig               13124
# Up                       0
toptagst2 <- topTags(DEA_D7_T2, n = Inf)
View(toptagst2$table)
write.csv2(toptagst2$table, "Results/EdgeR/Placebo/DEA_D7_T2", row.names = FALSE)

# T3
DEA_D7_T3 <- glmTreat(qlfit_new, coef =6, lfc = log2(1.2))
summary(decideTestsDGE(DEA_D7_T3))
# descript_gag$DoseT3
# Down                     0
# NotSig               13124
# Up                       0
toptagst3 <- topTags(DEA_D7_T3, n = Inf)
View(toptagst3$table)
write.csv2(toptagst3$table, "Results/EdgeR/Placebo/DEA_D7_T3", row.names = FALSE)

# T4
DEA_D7_T4 <- glmTreat(qlfit_new, coef =7, lfc = log2(1.2))
summary(decideTestsDGE(DEA_D7_T4))
# descript_gag$DoseT4
# Down                     0
# NotSig               13122
# Up                       2
toptagst4 <- topTags(DEA_D7_T4, n = Inf)
View(toptagst4$table)
write.csv2(toptagst4$table, "Results/EdgeR/Placebo/DEA_D7_T4", row.names = FALSE)

# T5
DEA_D7_T5 <- glmTreat(qlfit_new, coef =8, lfc = log2(1.2))
summary(decideTestsDGE(DEA_D7_T5))
# descript_gag$DoseT5
# Down                     0
# NotSig               13124
# Up                       0
toptagst5 <- topTags(DEA_D7_T5, n = Inf)
View(toptagst5$table)
write.csv2(toptagst5$table, "Results/EdgeR/Placebo/DEA_D7_T5", row.names = FALSE)

####
####
#### Find DEGs Day 14
# T1
DEA_D14_T1 <- glmTreat(qlfit_new, coef =4, lfc = log2(1.2))
summary(decideTestsDGE(DEA_D14_T1))
# descript_gag$DoseT1
# Down                     0
# NotSig               13185
# Up                       0
toptagst1 <- topTags(DEA_D14_T1, n = Inf)
View(toptagst1$table)
write.csv2(toptagst1$table, "Results/EdgeR/Placebo/DEA_D14_T1", row.names = FALSE)

# T2
DEA_D14_T2 <- glmTreat(qlfit_new, coef =5, lfc = log2(1.2))
summary(decideTestsDGE(DEA_D14_T2))
# descript_gag$DoseT2
# Down                     0
# NotSig               13185
# Up                       0
toptagst2 <- topTags(DEA_D14_T2, n = Inf)
View(toptagst2$table)
write.csv2(toptagst2$table, "Results/EdgeR/Placebo/DEA_D14_T2", row.names = FALSE)

# T3
DEA_D14_T3 <- glmTreat(qlfit_new, coef =6, lfc = log2(1.2))
summary(decideTestsDGE(DEA_D14_T3))
# descript_gag$DoseT3
# Down                     0
# NotSig               13124
# Up                       0
toptagst3 <- topTags(DEA_D14_T3, n = Inf)
View(toptagst3$table)
write.csv2(toptagst3$table, "Results/EdgeR/Placebo/DEA_D14_T3", row.names = FALSE)

# T4
DEA_D14_T4 <- glmTreat(qlfit_new, coef =7, lfc = log2(1.2))
summary(decideTestsDGE(DEA_D14_T4))
# descript_gag$DoseT4
# Down                     0
# NotSig               13185
# Up                       0
toptagst4 <- topTags(DEA_D14_T4, n = Inf)
View(toptagst4$table)
write.csv2(toptagst4$table, "Results/EdgeR/Placebo/DEA_D14_T4", row.names = FALSE)

# T5
DEA_D14_T5 <- glmTreat(qlfit_new, coef =8, lfc = log2(1.2))
summary(decideTestsDGE(DEA_D14_T5))
# descript_gag$DoseT5
# Down                     0
# NotSig               13185
# Up                       0
toptagst5 <- topTags(DEA_D14_T5, n = Inf)
View(toptagst5$table)
write.csv2(toptagst5$table, "Results/EdgeR/Placebo/DEA_D14_T5", row.names = FALSE)


####
####
#### Find DEGs Day 63
# T1
DEA_D63_T1 <- glmTreat(qlfit_new, coef =4, lfc = log2(1.2))
summary(decideTestsDGE(DEA_D63_T1))
# descript_gag$DoseT1
# Down                     0
# NotSig               13107
# Up                       0
toptagst1 <- topTags(DEA_D63_T1, n = Inf)
View(toptagst1$table)
write.csv2(toptagst1$table, "Results/EdgeR/Placebo/DEA_D63_T1", row.names = FALSE)

# T2
DEA_D63_T2 <- glmTreat(qlfit_new, coef =5, lfc = log2(1.2))
summary(decideTestsDGE(DEA_D63_T2))
# descript_gag$DoseT2
# Down                     0
# NotSig               13107
# Up                       0
toptagst2 <- topTags(DEA_D63_T2, n = Inf)
View(toptagst2$table)
write.csv2(toptagst2$table, "Results/EdgeR/Placebo/DEA_D63_T2", row.names = FALSE)

# T3
DEA_D63_T3 <- glmTreat(qlfit_new, coef =6, lfc = log2(1.2))
summary(decideTestsDGE(DEA_D63_T3))
# descript_gag$DoseT3
# Down                     0
# NotSig               13107
# Up                       0
toptagst3 <- topTags(DEA_D63_T3, n = Inf)
View(toptagst3$table)
write.csv2(toptagst3$table, "Results/EdgeR/Placebo/DEA_D63_T3", row.names = FALSE)

# T4
DEA_D63_T4 <- glmTreat(qlfit_new, coef =7, lfc = log2(1.2))
summary(decideTestsDGE(DEA_D63_T4))
# descript_gag$DoseT4
# Down                     0
# NotSig               13107
# Up                       0
toptagst4 <- topTags(DEA_D63_T4, n = Inf)
View(toptagst4$table)
write.csv2(toptagst4$table, "Results/EdgeR/Placebo/DEA_D63_T4", row.names = FALSE)

# T5
DEA_D63_T5 <- glmTreat(qlfit_new, coef =8, lfc = log2(1.2))
summary(decideTestsDGE(DEA_D63_T5))
# descript_gag$DoseT5
# Down                     0
# NotSig               13107
# Up                       0
toptagst5 <- topTags(DEA_D63_T5, n = Inf)
View(toptagst5$table)
write.csv2(toptagst5$table, "Results/EdgeR/Placebo/DEA_D63_T5", row.names = FALSE)



##########################################################################
#### EdgeR for all doses ####

descript_gag <- read.csv2("data/VSV_gag_desc.csv")
count_gag <- read.csv2("data/VSV_gag_COUNTS.csv")
View(descript_gag)

rownames(descript_gag) <- paste0("X", descript_gag$sample)
rownames(count_gag) <- count_gag$Gene
count_gag$Gene <- NULL

descript_gag <- descript_gag[descript_gag$Treatment != "P",]
common <- intersect(colnames(count_gag), rownames(descript_gag))
count_gag <- count_gag[, common] 

identical(colnames(count_gag), rownames(descript_gag))
# create DEGlist
y_ebovac <- DGEList(counts = count_gag, genes = row.names(count_gag), group= descript_gag$Visit)

keep <- rowSums(cpm(y_ebovac)>1)>=10
y_ebovac <- y_ebovac[keep,]
dim(y_ebovac) # 13999   252

# normalize
y.1 <- calcNormFactors(y_ebovac)

# design matrix
batch <- factor(descript_gag$batch)
descript_gag$Visit <- relevel(as.factor(descript_gag$Visit), ref = "2")
design <- model.matrix(~batch+descript_gag$Visit)
y.2 <- estimateDisp(y.1, design)

# Fitting general linear models to each feature
qlfit_new <- glmQLFit(y.2, design)
colnames(design)

#### Find DEGs 
dir.create("Results/EdgeR", recursive = T)
# DAY 07
DEA_day07 <- glmTreat(qlfit_new, coef =4, lfc = log2(1))
summary(decideTestsDGE(DEA_day07))
# Down                     0
# NotSig               13647
# Up                      47
toptags07 <- topTags(DEA_day07, n = Inf)
View(toptags07$table)
write.csv2(toptags07$table, "Results/EdgeR/DEGs_07", row.names = FALSE)

# DAY 14
DEA_day14 <- glmTreat(qlfit_new, coef =5, lfc = log2(1))
summary(decideTestsDGE(DEA_day14))
# Down                    18
# NotSig               13610
# Up                      66
toptags14 <- topTags(DEA_day14, n = Inf)
View(toptags14$table)
write.csv2(toptags14$table, "Results/EdgeR/DEGs_14", row.names = FALSE)

# DAY 63
DEA_day63 <- glmTreat(qlfit_new, coef =6, lfc = log2(1))
summary(decideTestsDGE(DEA_day63))
# Down                     5
# NotSig               13674
# Up                      15
toptags63 <- topTags(DEA_day63, n = Inf)
View(toptags63$table)
write.csv2(toptags63$table, "Results/EdgeR/DEGs_63", row.names = FALSE)

####

### tmod
dir.create("Results/tmod")

d7 <- read.csv2("Results/EdgeR/DEGs_07")
View(martd7)
index_d7 <- order(d7$FDR)
martd7 <- d7[index_d7,]
d7 <- tmodCERNOtest(martd7$genes)
View(d7)
write.csv2(d7, file = "Results/tmod/d7.csv", row.names = F)

###
d14 <- read.csv2("Results/EdgeR/DEGs_14")
index_d14 <- order(d14$FDR)
martd14 <- d14[index_d14,]
d14 <- tmodCERNOtest(martd14$genes)
View(d14)
write.csv2(d14, file = "Results/tmod/d14.csv", row.names = F)

###
d63 <- read.csv2("Results/EdgeR/DEGs_63")
index_d63 <- order(d63$FDR)
martd63 <- d63[index_d63,]
d63 <- tmodCERNOtest(martd63$genes)
View(d63)
write.csv2(d63, file = "Results/tmod/d63.csv", row.names = F)




#### EdgeR T1 ####

descript_gag <- read.csv2("data/VSV_gag_desc.csv")
count_gag <- read.csv2("data/VSV_gag_COUNTS.csv")

rownames(descript_gag) <- paste0("X", descript_gag$sample)
rownames(count_gag) <- count_gag$Gene
count_gag$Gene <- NULL

descript_gag <- descript_gag[descript_gag$Treatment != "P",]
descript_gag <- descript_gag[descript_gag$Visit %in% c(2,)]


common <- intersect(colnames(count_gag), rownames(descript_gag))
count_gag <- count_gag[, common] 

identical(colnames(count_gag), rownames(descript_gag))
# create DEGlist
y_ebovac <- DGEList(counts = count_gag, genes = row.names(count_gag), group= descript_gag$Visit)

keep <- rowSums(cpm(y_ebovac)>1)>=10
y_ebovac <- y_ebovac[keep,]
dim(y_ebovac) # 13999   252

# normalize
y.1 <- calcNormFactors(y_ebovac)

# design matrix
batch <- factor(descript_gag$batch)
descript_gag$Visit <- relevel(as.factor(descript_gag$Visit), ref = "2")
design <- model.matrix(~batch+descript_gag$Visit)
y.2 <- estimateDisp(y.1, design)

# Fitting general linear models to each feature
qlfit_new <- glmQLFit(y.2, design)
colnames(design)

#### Find DEGs 
dir.create("Results/EdgeR", recursive = T)
# DAY 07
DEA_day07 <- glmTreat(qlfit_new, coef =4, lfc = log2(1))
summary(decideTestsDGE(DEA_day07))
# Down                     0
# NotSig               13647
# Up                      47
toptags07 <- topTags(DEA_day07, n = Inf)
View(toptags07$table)
write.csv2(toptags07$table, "Results/EdgeR/DEGs_07", row.names = FALSE)

# DAY 14
DEA_day14 <- glmTreat(qlfit_new, coef =5, lfc = log2(1))
summary(decideTestsDGE(DEA_day14))
# Down                    18
# NotSig               13610
# Up                      66
toptags14 <- topTags(DEA_day14, n = Inf)
View(toptags14$table)
write.csv2(toptags14$table, "Results/EdgeR/DEGs_14", row.names = FALSE)

# DAY 63
DEA_day63 <- glmTreat(qlfit_new, coef =6, lfc = log2(1))
summary(decideTestsDGE(DEA_day63))
# Down                     5
# NotSig               13674
# Up                      15
toptags63 <- topTags(DEA_day63, n = Inf)
View(toptags63$table)
write.csv2(toptags63$table, "Results/EdgeR/DEGs_63", row.names = FALSE)
