###### 
rm(list = ls())
options(stringsAsFactors = F)

sapply(c('ggplot2','ggrepel','dplyr','Glimma', "tidyr", "edgeR", "tmod"), require,character.only=T)

#######################################################################
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



#######################################################################
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


d7<-read.csv2("Results/EdgeR/All_treatments/DEGs_07")
d14<-read.csv2("Results/EdgeR/All_treatments/DEGs_14")
d63<-read.csv2("Results/EdgeR/All_treatments/DEGs_63")

D7<-read.csv2("Results/tmod/d7.csv")
D14<-read.csv2("Results/tmod/d14.csv")
D63<-read.csv2("Results/tmod/d63.csv")

creslist<-list(D7=D7, D14=D14, D63=D63)

tmodPanelPlot(creslist,pval.thr = 10^-2, pval.cutoff=10^-7, filter.unknown =  TRUE, 
              text.cex = 0.5,  legend.style = "tall" )

#######################################################################
#### EdgeR for each Dose #### 

descript_gag <- read.csv2("data/VSV_gag_desc.csv")
count_gag <- read.csv2("data/VSV_gag_COUNTS.csv")

rownames(descript_gag) <- paste0("X", descript_gag$sample)
rownames(count_gag) <- count_gag$Gene
count_gag$Gene <- NULL

descript_gag <- descript_gag[descript_gag$Treatment != "P",]
descript_gag <- descript_gag[descript_gag$Dose %in% c("T1", "T2"),]
View(descript_gag)

common <- intersect(colnames(count_gag), rownames(descript_gag))
count_gag <- count_gag[, common] 

identical(colnames(count_gag), rownames(descript_gag))
#

#### EdgeR T1 ####
# create DEGlist
y_ebovac <- DGEList(counts = count_gag, genes = row.names(count_gag), group= descript_gag$Visit)

keep <- rowSums(cpm(y_ebovac)>1)>=10
y_ebovac <- y_ebovac[keep,]
dim(y_ebovac) # 13002    40

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
# DAY 07
DEA_day07 <- glmTreat(qlfit_new, coef =4, lfc = log2(1))
summary(decideTestsDGE(DEA_day07))
# descript_gag$Visit4
# Down                     0
# NotSig               13735
# Up                      48
toptags07 <- topTags(DEA_day07, n = Inf)
View(toptags07$table)
write.csv2(toptags07$table, "Results/EdgeR/No_Placebo/T1/DEGs_T1_D7", row.names = FALSE)

# DAY 14
DEA_day14 <- glmTreat(qlfit_new, coef =5, lfc = log2(1))
summary(decideTestsDGE(DEA_day14))
# descript_gag$Visit5
# Down                    39
# NotSig               13653
# Up                      91
toptags14 <- topTags(DEA_day14, n = Inf)
View(toptags14$table)
write.csv2(toptags14$table, "Results/EdgeR/No_Placebo/T1/DEGs_T1_D14", row.names = FALSE)

# DAY 63
DEA_day63 <- glmTreat(qlfit_new, coef =6, lfc = log2(1))
summary(decideTestsDGE(DEA_day63))
# descript_gag$Visit8
# Down                    57
# NotSig               13636
# Up                      90
toptags63 <- topTags(DEA_day63, n = Inf)
View(toptags63$table)
write.csv2(toptags63$table, "Results/EdgeR/No_Placebo/T1/DEGs_T1_D63", row.names = FALSE)

#### EdgeR T2 ####

y_ebovac <- DGEList(counts = count_gag, genes = row.names(count_gag), group= descript_gag$Visit)

keep <- rowSums(cpm(y_ebovac)>1)>=10
y_ebovac <- y_ebovac[keep,]
dim(y_ebovac) # 12898    37

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
# DAY 07
DEA_day07 <- glmTreat(qlfit_new, coef =4, lfc = log2(1))
summary(decideTestsDGE(DEA_day07))
# descript_gag$Visit4
# Down                     0
# NotSig               12898
# Up                       0
toptags07 <- topTags(DEA_day07, n = Inf)
View(toptags07$table)
write.csv2(toptags07$table, "Results/EdgeR/No_Placebo/T2/DEGs_T2_D7", row.names = FALSE)

# DAY 14
DEA_day14 <- glmTreat(qlfit_new, coef =5, lfc = log2(1))
summary(decideTestsDGE(DEA_day14))
# descript_gag$Visit5
# Down                     0
# NotSig               12898
# Up                       0
toptags14 <- topTags(DEA_day14, n = Inf)
View(toptags14$table)
write.csv2(toptags14$table, "Results/EdgeR/No_Placebo/T2/DEGs_T2_D14", row.names = FALSE)

# DAY 63
DEA_day63 <- glmTreat(qlfit_new, coef =6, lfc = log2(1))
summary(decideTestsDGE(DEA_day63))
# descript_gag$Visit8
# Down                     0
# NotSig               12898
# Up                       0
toptags63 <- topTags(DEA_day63, n = Inf)
View(toptags63$table)
write.csv2(toptags63$table, "Results/EdgeR/No_Placebo/T2/DEGs_T2_D63", row.names = FALSE)

#### EdgeR T3 ####

y_ebovac <- DGEList(counts = count_gag, genes = row.names(count_gag), group= descript_gag$Visit)

keep <- rowSums(cpm(y_ebovac)>1)>=10
y_ebovac <- y_ebovac[keep,]
dim(y_ebovac) # 12907    40

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
# DAY 07
DEA_day07 <- glmTreat(qlfit_new, coef =4, lfc = log2(1))
summary(decideTestsDGE(DEA_day07))
# descript_gag$Visit4
# Down                     0
# NotSig               12907
# Up                       0
toptags07 <- topTags(DEA_day07, n = Inf)
View(toptags07$table)
write.csv2(toptags07$table, "Results/EdgeR/No_Placebo/T3/DEGs_T3_D7", row.names = FALSE)

# DAY 14
DEA_day14 <- glmTreat(qlfit_new, coef =5, lfc = log2(1))
summary(decideTestsDGE(DEA_day14))
# descript_gag$Visit5
# Down                     0
# NotSig               12907
# Up                       0
toptags14 <- topTags(DEA_day14, n = Inf)
View(toptags14$table)
write.csv2(toptags14$table, "Results/EdgeR/No_Placebo/T3/DEGs_T3_D14", row.names = FALSE)

# DAY 63
DEA_day63 <- glmTreat(qlfit_new, coef =6, lfc = log2(1))
summary(decideTestsDGE(DEA_day63))
# descript_gag$Visit8
# Down                     0
# NotSig               12907
# Up                       0
toptags63 <- topTags(DEA_day63, n = Inf)
View(toptags63$table)
write.csv2(toptags63$table, "Results/EdgeR/No_Placebo/T3/DEGs_T3_D63", row.names = FALSE)


#### EdgeR T4 ####

y_ebovac <- DGEList(counts = count_gag, genes = row.names(count_gag), group= descript_gag$Visit)

keep <- rowSums(cpm(y_ebovac)>1)>=10
y_ebovac <- y_ebovac[keep,]
dim(y_ebovac) #  12785    32

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
# DAY 07
DEA_day07 <- glmTreat(qlfit_new, coef =4, lfc = log2(1))
summary(decideTestsDGE(DEA_day07))
# descript_gag$Visit4
# Down                     0
# NotSig               12785
# Up                       0
toptags07 <- topTags(DEA_day07, n = Inf)
View(toptags07$table)
write.csv2(toptags07$table, "Results/EdgeR/No_Placebo/T4/DEGs_T4_D7", row.names = FALSE)

# DAY 14
DEA_day14 <- glmTreat(qlfit_new, coef =5, lfc = log2(1))
summary(decideTestsDGE(DEA_day14))
# descript_gag$Visit5
# Down                     0
# NotSig               12785
# Up                       0
toptags14 <- topTags(DEA_day14, n = Inf)
View(toptags14$table)
write.csv2(toptags14$table, "Results/EdgeR/No_Placebo/T4/DEGs_T4_D14", row.names = FALSE)

# DAY 63
DEA_day63 <- glmTreat(qlfit_new, coef =6, lfc = log2(1))
summary(decideTestsDGE(DEA_day63))
# descript_gag$Visit8
# Down                     0
# NotSig               12785
# Up                       0
toptags63 <- topTags(DEA_day63, n = Inf)
View(toptags63$table)
write.csv2(toptags63$table, "Results/EdgeR/No_Placebo/T4/DEGs_T4_D63", row.names = FALSE)


#### EdgeR T5 ####

y_ebovac <- DGEList(counts = count_gag, genes = row.names(count_gag), group= descript_gag$Visit)

keep <- rowSums(cpm(y_ebovac)>1)>=10
y_ebovac <- y_ebovac[keep,]
dim(y_ebovac) #  12874    35

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
# DAY 07
DEA_day07 <- glmTreat(qlfit_new, coef =4, lfc = log2(1))
summary(decideTestsDGE(DEA_day07))
# descript_gag$Visit4
# Down                     0
# NotSig               12873
# Up                       1
toptags07 <- topTags(DEA_day07, n = Inf)
View(toptags07$table)
write.csv2(toptags07$table, "Results/EdgeR/No_Placebo/T5/DEGs_T5_D7", row.names = FALSE)

# DAY 14
DEA_day14 <- glmTreat(qlfit_new, coef =5, lfc = log2(1))
summary(decideTestsDGE(DEA_day14))
# descript_gag$Visit5
# Down                     0
# NotSig               12874
# Up                       0
toptags14 <- topTags(DEA_day14, n = Inf)
View(toptags14$table)
write.csv2(toptags14$table, "Results/EdgeR/No_Placebo/T5/DEGs_T5_D14", row.names = FALSE)

# DAY 63
DEA_day63 <- glmTreat(qlfit_new, coef =6, lfc = log2(1))
summary(decideTestsDGE(DEA_day63))
# descript_gag$Visit8
# Down                     0
# NotSig               12874
# Up                       0
toptags63 <- topTags(DEA_day63, n = Inf)
View(toptags63$table)
write.csv2(toptags63$table, "Results/EdgeR/No_Placebo/T5/DEGs_T5_D63", row.names = FALSE)


#### EdgeR T4 and T5 ####

y_ebovac <- DGEList(counts = count_gag, genes = row.names(count_gag), group= descript_gag$Visit)

keep <- rowSums(cpm(y_ebovac)>1)>=10
y_ebovac <- y_ebovac[keep,]
dim(y_ebovac) #  13225    67

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
# DAY 07
DEA_day07 <- glmTreat(qlfit_new, coef =4, lfc = log2(1))
summary(decideTestsDGE(DEA_day07))
# descript_gag$Visit4
# Down                     1
# NotSig               13182
# Up                      42
toptags07 <- topTags(DEA_day07, n = Inf)
View(toptags07$table)
write.csv2(toptags07$table, "Results/EdgeR/No_Placebo/T45/DEGs_T45_D7", row.names = FALSE)

# DAY 14
DEA_day14 <- glmTreat(qlfit_new, coef =5, lfc = log2(1))
summary(decideTestsDGE(DEA_day14))
# descript_gag$Visit5
# Down                     0
# NotSig               13225
# Up                       0
toptags14 <- topTags(DEA_day14, n = Inf)
View(toptags14$table)
write.csv2(toptags14$table, "Results/EdgeR/No_Placebo/T45/DEGs_T45_D14", row.names = FALSE)

# DAY 63
DEA_day63 <- glmTreat(qlfit_new, coef =6, lfc = log2(1))
summary(decideTestsDGE(DEA_day63))
# descript_gag$Visit8
# Down                     0
# NotSig               13225
# Up                       0
toptags63 <- topTags(DEA_day63, n = Inf)
View(toptags63$table)
write.csv2(toptags63$table, "Results/EdgeR/No_Placebo/T45/DEGs_T5_D63", row.names = FALSE)



#### EdgeR T3, T4 and T5 ####

y_ebovac <- DGEList(counts = count_gag, genes = row.names(count_gag), group= descript_gag$Visit)

keep <- rowSums(cpm(y_ebovac)>1)>=10
y_ebovac <- y_ebovac[keep,]
dim(y_ebovac) #   13450   107

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
# DAY 07
DEA_day07 <- glmTreat(qlfit_new, coef =4, lfc = log2(1))
summary(decideTestsDGE(DEA_day07))
# descript_gag$Visit4
# Down                     2
# NotSig               13401
# Up                      47
toptags07 <- topTags(DEA_day07, n = Inf)
View(toptags07$table)
write.csv2(toptags07$table, "Results/EdgeR/No_Placebo/T345/DEGs_T345_D7", row.names = FALSE)

# DAY 14
DEA_day14 <- glmTreat(qlfit_new, coef =5, lfc = log2(1))
summary(decideTestsDGE(DEA_day14))
# descript_gag$Visit5
# Down                     0
# NotSig               13446
# Up                       4
toptags14 <- topTags(DEA_day14, n = Inf)
View(toptags14$table)
write.csv2(toptags14$table, "Results/EdgeR/No_Placebo/T345/DEGs_T345_D14", row.names = FALSE)

# DAY 63
DEA_day63 <- glmTreat(qlfit_new, coef =6, lfc = log2(1))
summary(decideTestsDGE(DEA_day63))
# descript_gag$Visit8
# Down                     0
# NotSig               13450
# Up                       0
toptags63 <- topTags(DEA_day63, n = Inf)
View(toptags63$table)
write.csv2(toptags63$table, "Results/EdgeR/No_Placebo/T345/DEGs_T345_D63", row.names = FALSE)


##### TMOD ####
library(tmod)
# T4 and T5
#dir.create("Results/tmod/T4_T5")

d7 <- read.csv2("Results/EdgeR/No_Placebo/T45/DEGs_T45_D7")
index_d7 <- order(d7$FDR)
martd7 <- d7[index_d7,]
d7 <- tmodCERNOtest(martd7$genes)
View(d7)
write.csv2(d7, file = "Results/tmod/T4_T5/T4_T5_d7.csv", row.names = F)

###
d14 <- read.csv2("Results/EdgeR/No_Placebo/T45/DEGs_T45_D14")
index_d14 <- order(d14$FDR)
martd14 <- d14[index_d14,]
d14 <- tmodCERNOtest(martd14$genes)
View(d14)
write.csv2(d14, file = "Results/tmod/T4_T5/T4_T5_d14.csv", row.names = F)

###
d63 <- read.csv2("Results/EdgeR/No_Placebo/T45/DEGs_T45_D63")
index_d63 <- order(d63$FDR)
martd63 <- d63[index_d63,]
d63 <- tmodCERNOtest(martd63$genes)
View(d63)
write.csv2(d63, file = "Results/tmod/d63.csv", row.names = F)


#################
t123_7 <- read.csv2("Results/EdgeR/No_Placebo/T123/DEGs_T123_D7")
t45_7 <- read.csv2("Results/EdgeR/No_Placebo/T45/DEGs_T45_D7")
t123_14 <- read.csv2("Results/EdgeR/No_Placebo/T123/DEGs_T123_D14")
t45_14 <- read.csv2("Results/EdgeR/No_Placebo/T45/DEGs_T45_D14")
t123_63 <- read.csv2("Results/EdgeR/No_Placebo/T123/DEGs_T123_D63")
t45_63 <- read.csv2("Results/EdgeR/No_Placebo/T45/DEGs_T45_D63")
#View(t123_63)

T123_D7 <- tmodCERNOtest(t123_7$genes, mset ="LI")
T45_D7 <- tmodCERNOtest(t45_7$genes, mset ="LI")
T123_D14 <- tmodCERNOtest(t123_14$genes, mset ="LI")
T45_D14 <- tmodCERNOtest(t45_14$genes, mset ="LI")
T123_D63 <- tmodCERNOtest(t123_63$genes, mset ="LI")
T45_D63 <- tmodCERNOtest(t45_63$genes, mset ="LI")

creslist <- list(T123_D7=T123_D7,T45_D7=T45_D7, T123_D14=T123_D14, T45_D14=T45_D14,
               T123_D63=T123_D63, T45_D63=T45_D63)

tmodPanelPlot(creslist,pval.thr = 10^-2, pval.cutoff=10^-7, filter.unknown =  TRUE, 
              text.cex = 0.7,  legend.style = "tall" )


#################

pieD7<- tmodDecideTests(martd7$genes, lfc=martd7$logFC, pval = martd7$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

pieD14<- tmodDecideTests(martd14$genes, lfc=martd14$logFC, pval = martd14$FDR, 
                         mset = "LI", pval.thr = 0.05, lfc.thr = 0)

pieD63<- tmodDecideTests(martd63$genes, lfc=martd63$logFC, pval = martd63$FDR, 
                         mset = "LI", pval.thr = 0.05, lfc.thr = 0)


pieA<- as.data.frame(pieD7)
pieB<- as.data.frame(pieD14)
pieC<- as.data.frame(pieD63)

colnames(pieA) <- gsub("X.*\\.", "", colnames(pieA))
colnames(pieB) <- gsub("X.*\\.", "", colnames(pieB))
colnames(pieC) <- gsub("X.*\\.", "", colnames(pieC))

D7_LI<- list(d7)
D14_LI<- list(d14)
D63_LI<- list(d63)


panel_list_7_14_63 <- c(D7_LI, D14_LI, D63_LI)
is.list(panel_list_7_14_63)
pie_list_7_14_63 <- c(pieA, pieB, pieC)

tmodPanelPlot(panel_list_7_14_63, pval.thr = 10^-3, 
              pval.cutoff= 10^-30, filter.unknown = T, 
              text.cex =0.60, clust = "qval", pie =pie_list_7_14_63, pie.style = "pie")


#### EdgeR T1, T2 and T3 ####

y_ebovac <- DGEList(counts = count_gag, genes = row.names(count_gag), group= descript_gag$Visit)

keep <- rowSums(cpm(y_ebovac)>1)>=10
y_ebovac <- y_ebovac[keep,]
dim(y_ebovac) # 13519   117

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

dir.create("Results/EdgeR/No_Placebo/T123")

#### Find DEGs 
# DAY 07
DEA_day07 <- glmTreat(qlfit_new, coef =4, lfc = log2(1))
summary(decideTestsDGE(DEA_day07))
# descript_gag$Visit4
# Down                     0
# NotSig               13518
# Up                       1
toptags07 <- topTags(DEA_day07, n = Inf)
View(toptags07$table)
write.csv2(toptags07$table, "Results/EdgeR/No_Placebo/T123/DEGs_T123_D7", row.names = FALSE)

# DAY 14
DEA_day14 <- glmTreat(qlfit_new, coef =5, lfc = log2(1))
summary(decideTestsDGE(DEA_day14))
# descript_gag$Visit5
# Down                   104
# NotSig               13345
# Up                      70
toptags14 <- topTags(DEA_day14, n = Inf)
View(toptags14$table)
write.csv2(toptags14$table, "Results/EdgeR/No_Placebo/T123/DEGs_T123_D14", row.names = FALSE)

# DAY 63
DEA_day63 <- glmTreat(qlfit_new, coef =6, lfc = log2(1))
summary(decideTestsDGE(DEA_day63))
# descript_gag$Visit8
# Down                     0
# NotSig               13519
# Up                       0
toptags63 <- topTags(DEA_day63, n = Inf)
View(toptags63$table)
write.csv2(toptags63$table, "Results/EdgeR/No_Placebo/T123/DEGs_T123_D63", row.names = FALSE)





#### EdgeR T1 and T2 ####

y_ebovac <- DGEList(counts = count_gag, genes = row.names(count_gag), group= descript_gag$Visit)

keep <- rowSums(cpm(y_ebovac)>1)>=10
y_ebovac <- y_ebovac[keep,]
dim(y_ebovac) # 13519   117

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

dir.create("Results/EdgeR/No_Placebo/T12")

#### Find DEGs 
# DAY 07
DEA_day07 <- glmTreat(qlfit_new, coef =4, lfc = log2(1))
summary(decideTestsDGE(DEA_day07))
# descript_gag$Visit4
# Down                     0
# NotSig               13329
# Up                       0
toptags07 <- topTags(DEA_day07, n = Inf)
View(toptags07$table)
write.csv2(toptags07$table, "Results/EdgeR/No_Placebo/T12/DEGs_T12_D7", row.names = FALSE)

# DAY 14
DEA_day14 <- glmTreat(qlfit_new, coef =5, lfc = log2(1))
summary(decideTestsDGE(DEA_day14))
# descript_gag$Visit5
# Down                     0
# NotSig               13329
# Up                       0
toptags14 <- topTags(DEA_day14, n = Inf)
View(toptags14$table)
write.csv2(toptags14$table, "Results/EdgeR/No_Placebo/T12/DEGs_T12_D14", row.names = FALSE)

# DAY 63
DEA_day63 <- glmTreat(qlfit_new, coef =6, lfc = log2(1))
summary(decideTestsDGE(DEA_day63))
# descript_gag$Visit8
# Down                     0
# NotSig               13329
# Up                       0
toptags63 <- topTags(DEA_day63, n = Inf)
View(toptags63$table)
write.csv2(toptags63$table, "Results/EdgeR/No_Placebo/T12/DEGs_T12_D63", row.names = FALSE)




