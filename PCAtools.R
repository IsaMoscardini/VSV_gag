rm(list = ls())
options(stringsAsFactors = F)


descript_gag <- read.csv2("data/VSV_gag_desc.csv")
count_gag <- read.csv2("data/VSV_gag_COUNTS.csv")

Anti <- read.csv("data/Antibod.csv")
Anam <- read.csv("data/Anam.csv")
View(Anti)

descript_gag$Cod <- Anam$Participant.ID[match(descript_gag$nID, Anam$ID)]
descript_gag$Age <- Anam$AGE..at.form.completion.[match(descript_gag$nID, Anam$ID)]
descript_gag$Gender <- Anam$Current.Gender[match(descript_gag$nID, Anam$ID)]
descript_gag$Pctpos_adj <- Anti$pctpos_adj[match(descript_gag$Cod, Anti$Participant.ID)]
descript_gag$Resp_Tcell <- Anti$RESPONSE[match(descript_gag$Cod, Anti$Participant.ID)]

ID_art <- c(122, 38, 52, 71, 203, 36, 126, 275)
descript_gag$Arthralgia <- ifelse(descript_gag$nID %in% ID_art, "Y", "N")

View(descript_gag)

