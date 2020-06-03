# VSV_gag

rm(list = ls())
options(stringsAsFactors = F)

sapply(c('ggplot2','ggrepel','dplyr','Glimma', "tidyr", "edgeR", "tmod", "mixOmics"), require,character.only=T)


#### Quality Control ####

descript_gag <- read.csv2("data/VSV_gag_desc.csv")
count_gag <- read.csv2("data/VSV_gag_COUNTS.csv")
View(descript_gag)

rownames(descript_gag) <- paste0("X", descript_gag$sample)
rownames(count_gag) <- count_gag$Gene
count_gag$Gene <- NULL
identical(colnames(count_gag), rownames(descript_gag))

# PCA
tcount <- as.data.frame(t(count_gag))
tcount[] <- apply(tcount, 2, as.numeric)
View(tcount)
d <- which(apply(tcount, 2, var)==0)
tcount <- tcount[ , which(apply(tcount, 2, var) != 0)]

trans <- pca(tcount, ncomp = 10, center = TRUE, scale = TRUE)

plotIndiv(trans, comp = c(1, 2), ind.names = TRUE, 
          group = descript_gag$response_p24_visit9,
          legend = TRUE, title = 'VSV_gag batch')

##### Multi assay Experiment ####

#BiocManager::install("MultiAssayExperiment")
library(MultiAssayExperiment)
library(S4Vectors)

# Read tables
rna <- count_gag
coldata <- descript_gag
ac <- ("Supp/AntibodyTier.xml")

## Combine to a named list and call the ExperimentList constructor function
assayList <- list(Affy = exprdat, Methyl450k = methyldat, RNASeqGene = rnadat,
                  GISTIC = rangeSE)
## Use the ExperimentList constructor
ExpList <- ExperimentList(assayList)

VSVgag_list <- ExperimentList()





## 


