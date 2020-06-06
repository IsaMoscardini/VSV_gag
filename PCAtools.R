rm(list = ls())
options(stringsAsFactors = F)

descript_gag <- read.csv2("data/VSV_gag_desc.csv")
count_gag <- read.csv2("data/VSV_gag_COUNTS.csv")

Anti <- read.csv("data/Antibod.csv")
Anam <- read.csv("data/Anam.csv")
#View(Anti)

descript_gag$Cod <- Anam$Participant.ID[match(descript_gag$nID, Anam$ID)]
descript_gag$Age <- Anam$AGE..at.form.completion.[match(descript_gag$nID, Anam$ID)]
descript_gag$Gender <- Anam$Current.Gender[match(descript_gag$nID, Anam$ID)]
descript_gag$Pctpos_adj <- Anti$pctpos_adj[match(descript_gag$Cod, Anti$Participant.ID)]
id_tcell <- c(126, 275, 319, 435, 25, 225, 175, 209)
descript_gag$Resp_Tcell <- ifelse(descript_gag$nID %in% id_tcell, "Y", "N")

ID_art <- c(122, 38, 52, 71, 203, 36, 126, 275)
descript_gag$Arthralgia <- ifelse(descript_gag$nID %in% ID_art, "Y", "N")

ID_


rownames(descript_gag) <- paste0("X", descript_gag$sample)

rownames(count_gag) <- count_gag$Gene

descript_gag$X <- NULL
count_gag$Gene <- NULL
#View(count_gag)
identical(colnames(count_gag), rownames(descript_gag))


#### PCA ####
descript_gag <- descript_gag[descript_gag$Visit == 8,]
com <- intersect(colnames(count_gag), rownames(descript_gag))
count_gag <- count_gag[,com]

y <- DGEList(counts = count_gag, genes = row.names(count_gag), group= descript_gag$Visit)
keep <- rowSums(cpm(y)>2)>=10
y.1 <- y[keep,]
dim(y.1) # 12774   221

y.1 <- calcNormFactors(y.1)

#library(ggplot2)
#library(vegan)

cp <- log2(cpm(y.1)+1)
sample.var <- apply(cp, 2, var)
sample.var <- sort(sample.var)
pca<-prcomp(t(cp))

pcapanel<- as.data.frame(pca$x)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar<- round(100* percentVar)

#group.colors <- c(NI_NS = "#66C2A5", NI_T = "#8DA0CB", I_NS ="#FC8D62", I_T = "#E78AC3")

ggplot(pcapanel, aes(PC2, PC3, color=as.factor(descript_gag$Resp_Tcell), shape = as.character(descript_gag$Gender)))+
  geom_point(size = 5)+ theme_minimal()+ theme(legend.position = "top")+ #labs(shape="Time Point", group = "group")+
  geom_text(label=descript_gag$nID, nudge_x = 0.9, nudge_y = 1)+
  xlab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% variance"))

load <- abs(pca$rotation)
View(load)

descript_gag$Gender <- gsub("M", "2", descript_gag$Gender)
descript_gag$Gender <- gsub("F", "1", descript_gag$Gender)

plot(descript_gag$Gender, descript_gag$titer_p24_visit9)

#### PCA tools ####

# factors variables
factor_vars <- c('batch','Dose','Treatment','PFU','response_p24_visit9','response_vsv_visit5',
                 'response_vsv_visit9','Gender','Resp_Tcell','Arthralgia')
descript_gag[,factor_vars] <- lapply(descript_gag[,factor_vars], as.factor)


# numeric variables
num_vars <- c('Visit','titer_p24_visit9','titer_vsv_visit5','titer_vsv_visit9','Age',
              'Pctpos_adj')
descript_gag[,num_vars] <- lapply(descript_gag[,num_vars], as.factor)

descript_gag[,c('Sample','nID', 'Cod')] <- NULL

#--- Tidy counts data
# Normalisation w/ DESeq2 (VST)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=count_gag,colData=descript_gag,design=~ 1)
dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds,normalized = T)>1)>=10;table(keep)
norm_data <- assay(vst(dds[keep,]))

identical(colnames(norm_data),rownames(descript_gag))

## PCA plot analysis ##
library(PCAtools)
library(pca3d)
p <- pca(norm_data, metadata = descript_gag, removeVar = 0.1)

# find optimal number of components
horn <- parallelPCA(norm_data)
elbow <- findElbowPoint(p$variance)
screeplot(p, components = getComponents(p, 1:20), vline = c(horn$n, elbow)) +
  geom_text(aes(horn$n + 1, 50, label = "Horn's", vjust = -1)) +
  geom_text(aes(elbow + 1, 50, label = "Elbow", vjust = -1))

screeplot(p)
biplot(p)
plotloadings(p)
eigencorplot(p, metavars=c(num_vars,factor_vars))
