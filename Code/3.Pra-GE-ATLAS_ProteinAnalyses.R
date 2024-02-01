#####################################
# Pra-GE-ATLAS - Proteomic Analyses #
#####################################

#### 0.Load libraries and pkgs ####

library(pRocessomics)
library(sva)
library(limma)
library(basilisk)
reticulate::use_python("C:\\Users\\bboyl\\AppData\\Local\\Programs\\Python\\Python38")
library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(WGCNA)
library(GO.db)
library(myTAI)
library(RColorBrewer)
library(factoextra)
library(cdata)
library(scatterplot3d)
library(uwot)
library(UpSetR)
library(circlize)
library(ComplexHeatmap)
library(wesanderson)
library(myTAI)
library(preprocessCore)
library(rgexf)

#### 1.Data Prepare ####

#### 1.1 Import data ####

Pra_Tissues <- read.delim("clipboard", header = T)
Pra_FU_Total <- read.delim("clipboard", header = T)
Pra_UV_Total <- read.delim("clipboard", header = T)
Pra_UV_Nucleus <- read.delim("clipboard", header = T)
Pra_UV_Chloroplast <- read.delim("clipboard", header = T)
Pra_HS_Total <- read.delim("clipboard", header = T)
Pra_HS_Nucleus <- read.delim("clipboard", header = T)
Pra_HS_Chloroplast <- read.delim("clipboard", header = T)

#### 1.2 filter PD table and Colnames+TranscriptID ####

## filter PD table
Protein_list <- list(Tissues = Pra_Tissues,
                     FU_Total = Pra_FU_Total,
                     HS_Total = Pra_HS_Total,
                     UV_Total = Pra_UV_Total,
                     HS_Nucleus = Pra_HS_Nucleus,
                     UV_Nucleus = Pra_UV_Nucleus,
                     HS_Chloroplast = Pra_HS_Chloroplast,
                     UV_Chloroplast = Pra_UV_Chloroplast)
lapply(Protein_list,nrow)
Protein_list.filtered <- lapply(Protein_list, function(x){
  hits <- which(x$Coverage.... >= 5 & x$Exp..q.value..Combined <= 0.05 & x$X..Peptides >= 2 & x$X..Unique.Peptides >= 1 | x$Coverage.... >= 5 & x$Exp..q.value..Combined <= 0.05 & x$X..AAs < 50)
  x[hits,]
  })

## Subset interesting columns for downstream, rename it and ProtID to TranscriptID
lapply(Protein_list.filtered, colnames)
Protein_list.filtered$Tissues <- Protein_list.filtered$Tissues[,c(6,49:60)]
Protein_list.filtered$FU_Total <- Protein_list.filtered$FU_Total[,c(5,41:48)]
Protein_list.filtered$HS_Total <- Protein_list.filtered$HS_Total[,c(4,41:49)]
Protein_list.filtered$UV_Total <- Protein_list.filtered$UV_Total[,c(4,53:67)]
Protein_list.filtered$HS_Nucleus <- Protein_list.filtered$HS_Nucleus[,c(5,80:107)]
Protein_list.filtered$UV_Nucleus <- Protein_list.filtered$UV_Nucleus[,c(4,44:54)]
Protein_list.filtered$HS_Chloroplast <- Protein_list.filtered$HS_Chloroplast[,c(4,86:117)]
Protein_list.filtered$UV_Chloroplast <- Protein_list.filtered$UV_Chloroplast[,c(5,87:118)]
lapply(Protein_list.filtered, colnames)

colnames(Protein_list.filtered$Tissues) <- c("ProtID", "Tissues_AA_1", "Tissues_AA_2", "Tissues_AA_3",
                                             "Tissues_AJ_1", "Tissues_AJ_2", "Tissues_AJ_3",
                                             "Tissues_RJ_1", "Tissues_RJ_2", "Tissues_RJ_3",
                                             "Tissues_YF_1", "Tissues_YF_2", "Tissues_YF_3") 
colnames(Protein_list.filtered$FU_Total) <- c("ProtID", "FUTotal_C_1", "FUTotal_C_2", "FUTotal_C_3", "FUTotal_C_4",
                                              "FUTotal_FU_1", "FUTotal_FU_2", "FUTotal_FU_3", "FUTotal_FU_4") # import two replicates x treatment | 12 samples
colnames(Protein_list.filtered$HS_Total) <- c("ProtID", "HSTotal_C_1", "HSTotal_C_2", "HSTotal_C_3", 
                                              "HSTotal_T1_1", "HSTotal_T1_2",  "HSTotal_T1_3",
                                              "HSTotal_T3_1", "HSTotal_T3_2",  "HSTotal_T3_3") # import one replicate x treatment | 12
colnames(Protein_list.filtered$UV_Total) <- c("ProtID", "UVTotal_C_1", "UVTotal_C_2", "UVTotal_C_3", 
                                              "UVTotal_T1_1", "UVTotal_T1_2",  "UVTotal_T1_3",
                                              "UVTotal_T2_1", "UVTotal_T2_2",  "UVTotal_T2_3",
                                              "UVTotal_T3_1", "UVTotal_T3_2",  "UVTotal_T3_3",
                                              "UVTotal_R_1", "UVTotal_R_2", "UVTotal_R_3")
colnames(Protein_list.filtered$HS_Nucleus) <- c("ProtID", "HSNucleus_C_1", "HSNucleus_C_2", "HSNucleus_C_3",  "HSNucleus_C_4", 
                                              "HSNucleus_T1_1", "HSNucleus_T1_2",  "HSNucleus_T1_3", "HSNucleus_T1_4",
                                              "HSNucleus_T3_1", "HSNucleus_T3_2",  "HSNucleus_T3_3", "HSNucleus_T3_4",
                                              "HSNucleus_T5_1", "HSNucleus_T5_2",  "HSNucleus_T5_3", "HSNucleus_T5_4",
                                              "HSNucleus_T10_1", "HSNucleus_T10_2",  "HSNucleus_T10_3", "HSNucleus_T10_4",
                                              "HSNucleus_T5R_1", "HSNucleus_T5R_2",  "HSNucleus_T5R_3", "HSNucleus_T5R_4",
                                              "HSNucleus_CR_1", "HSNucleus_CR_2",  "HSNucleus_CR_3", "HSNucleus_CR_4")
colnames(Protein_list.filtered$UV_Nucleus) <- c("ProtID", "UVNucleus_T3_2", "UVNucleus_C_1", "UVNucleus_T1_2",  "UVNucleus_C_2", 
                                                "UVNucleus_C_3", "UVNucleus_T1_3",  "UVNucleus_T1_1", "UVNucleus_C_4",
                                                "UVNucleus_T1_4", "UVNucleus_T3_1",  "UVNucleus_T3_3") # imput one T3 | 12
# order columns by treatment
Protein_list.filtered$UV_Nucleus <- Protein_list.filtered$UV_Nucleus[,c(1,3,5,6,9,8,4,7,10,11,2,12)]
colnames(Protein_list.filtered$HS_Chloroplast) <- c("ProtID", 
                                                    "HSChloroplast_CE_1", "HSChloroplast_CE_2", "HSChloroplast_CE_3",  "HSChloroplast_CE_4", 
                                                    "HSChloroplast_CT_1", "HSChloroplast_CT_2",  "HSChloroplast_CT_3", "HSChloroplast_CT_4",
                                                    "HSChloroplast_TE1_1", "HSChloroplast_TE1_2",  "HSChloroplast_TE1_3", "HSChloroplast_TE1_4",
                                                    "HSChloroplast_TT1_1", "HSChloroplast_TT1_2",  "HSChloroplast_TT1_3", "HSChloroplast_TT1_4",
                                                    "HSChloroplast_TE3_1", "HSChloroplast_TE3_2",  "HSChloroplast_TE3_3", "HSChloroplast_TE3_4",
                                                    "HSChloroplast_TT3_1", "HSChloroplast_TT3_2",  "HSChloroplast_TT3_3", "HSChloroplast_TT3_4",
                                                    "HSChloroplast_TE5_1", "HSChloroplast_TE5_2",  "HSChloroplast_TE5_3", "HSChloroplast_TE5_4",
                                                    "HSChloroplast_TT5_1", "HSChloroplast_TT5_2",  "HSChloroplast_TT5_3", "HSChloroplast_TT5_4")
colnames(Protein_list.filtered$UV_Chloroplast) <- c("ProtID", 
                                                    "UVChloroplast_CE_1", "UVChloroplast_CE_2", "UVChloroplast_CE_3",  "UVChloroplast_CE_4", 
                                                    "UVChloroplast_CT_1", "UVChloroplast_CT_2",  "UVChloroplast_CT_3", "UVChloroplast_CT_4",
                                                    "UVChloroplast_TE1/2_1", "UVChloroplast_TE1/2_2",  "UVChloroplast_TE1/2_3", "UVChloroplast_TE1/2_4",
                                                    "UVChloroplast_TT1/2_1", "UVChloroplast_TT1/2_2",  "UVChloroplast_TT1/2_3", "UVChloroplast_TT1/2_4",
                                                    "UVChloroplast_TE1_1", "UVChloroplast_TE1_2",  "UVChloroplast_TE1_3", "UVChloroplast_TE1_4",
                                                    "UVChloroplast_TT1_1", "UVChloroplast_TT1_2",  "UVChloroplast_TT1_3", "UVChloroplast_TT1_4",
                                                    "UVChloroplast_TE2_1", "UVChloroplast_TE2_2",  "UVChloroplast_TE2_3", "UVChloroplast_TE2_4",
                                                    "UVChloroplast_TT2_1", "UVChloroplast_TT2_2",  "UVChloroplast_TT2_3", "UVChloroplast_TT2_4")
# New column Transcript ID; tranpose and add treatment column
DB <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/4.FinalAnnotation/0.final/Pra-GE-ATLAS_ReducedAnaltyses.txt", header = T, sep = "\t")
Protein_list.filtered.names <- lapply(Protein_list.filtered, function(x){
  x[,1] <- do.call(rbind, strsplit(x[,1], split = "_", fixed = T))[,1]
  x
})

s <- strsplit(DB$ProteinID, split = ",")
Transcript2Prot <- data.frame(TranscriptID = rep(DB$TranscriptID, sapply(s, length)), ProtID = unlist(s))
Transcript2Prot <- unique(Transcript2Prot)
# check all hits have Transcript-Correspondance
length(which((Protein_list.filtered.names$Tissues$ProtID %in% Transcript2Prot$ProtID) == TRUE))/nrow(Protein_list.filtered.names$Tissues)*100
length(which((Protein_list.filtered.names$FU_Total$ProtID %in% Transcript2Prot$ProtID) == TRUE))/nrow(Protein_list.filtered.names$FU_Total)*100
length(which((Protein_list.filtered.names$HS_Total$ProtID %in% Transcript2Prot$ProtID) == TRUE))/nrow(Protein_list.filtered.names$HS_Total)*100
length(which((Protein_list.filtered.names$UV_Total$ProtID %in% Transcript2Prot$ProtID) == TRUE))/nrow(Protein_list.filtered.names$UV_Total)*100
length(which((Protein_list.filtered.names$HS_Nucleus$ProtID %in% Transcript2Prot$ProtID) == TRUE))/nrow(Protein_list.filtered.names$HS_Nucleus)*100                                                                                                     
length(which((Protein_list.filtered.names$UV_Nucleus$ProtID %in% Transcript2Prot$ProtID) == TRUE))/nrow(Protein_list.filtered.names$UV_Nucleus)*100
length(which((Protein_list.filtered.names$HS_Chloroplast$ProtID %in% Transcript2Prot$ProtID) == TRUE))/nrow(Protein_list.filtered.names$HS_Chloroplast)*100
length(which((Protein_list.filtered.names$UV_Chloroplast$ProtID %in% Transcript2Prot$ProtID) == TRUE))/nrow(Protein_list.filtered.names$UV_Chloroplast)*100

Protein_list.filtered.names$Tissues$TranscriptID <- NA
Protein_list.filtered.names$FU_Total$TranscriptID <- NA
Protein_list.filtered.names$HS_Total$TranscriptID <- NA
Protein_list.filtered.names$UV_Total$TranscriptID <- NA
Protein_list.filtered.names$HS_Nucleus$TranscriptID <- NA
Protein_list.filtered.names$UV_Nucleus$TranscriptID <- NA
Protein_list.filtered.names$HS_Chloroplast$TranscriptID <- NA
Protein_list.filtered.names$UV_Chloroplast$TranscriptID <- NA

for(j in 1:length(Protein_list.filtered.names)){
  for(i in 1:nrow(Protein_list.filtered.names[[j]])){
    hit <- which(Transcript2Prot$ProtID == Protein_list.filtered.names[[j]]$ProtID[i])
    if(length(hit) == 1){
      Protein_list.filtered.names[[j]]$TranscriptID[i] <- Transcript2Prot$TranscriptID[hit]  
    }else if(length(hit) > 1){
      stop("WTF")
    }
  }
}

Protein_list.filtered.names.tranpose <- lapply(Protein_list.filtered.names, function(x){
  rownames(x) <- x[,1]
  x <- x[,-c(1,ncol(x))]
  fix <- t(x[])
  Treatment <- do.call(rbind,strsplit(rownames(fix), split = "_"))[,2]
  fixed.f <- cbind(Treatment, fix)
  fixed.f
})


#### 1.3 Pre-process data ####

write.table(Protein_list.filtered.names.tranpose$Tissues, "F:/ATLAS_Pra/Protein_module/Tissues2import.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(Protein_list.filtered.names.tranpose$FU_Total, "F:/ATLAS_Pra/Protein_module/FUTotal2import.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(Protein_list.filtered.names.tranpose$HS_Total, "F:/ATLAS_Pra/Protein_module/HSTotal2import.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(Protein_list.filtered.names.tranpose$UV_Total, "F:/ATLAS_Pra/Protein_module/UVTotal2import.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(Protein_list.filtered.names.tranpose$HS_Nucleus, "F:/ATLAS_Pra/Protein_module/HSNucleus2import.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(Protein_list.filtered.names.tranpose$UV_Nucleus, "F:/ATLAS_Pra/Protein_module/UVNucleus2import.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(Protein_list.filtered.names.tranpose$HS_Chloroplast, "F:/ATLAS_Pra/Protein_module/HSChloroplast2import.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(Protein_list.filtered.names.tranpose$UV_Chloroplast, "F:/ATLAS_Pra/Protein_module/UVChloroplast2import.txt", sep = "\t", quote = F, row.names = F, col.names = T)

Tissues <- read.delim("F:/ATLAS_Pra/Protein_module/Tissues2import.txt", header = T)
colnames(Tissues) <- gsub(pattern = ".", replacement = "-", x = colnames(Tissues), fixed = T)
FUTotal <- read.delim("F:/ATLAS_Pra/Protein_module/FUTotal2import.txt", header = T)
colnames(FUTotal) <- gsub(pattern = ".", replacement = "-", x = colnames(FUTotal), fixed = T)
HSTotal <- read.delim("F:/ATLAS_Pra/Protein_module/HSTotal2import.txt", header = T)
colnames(HSTotal) <- gsub(pattern = ".", replacement = "-", x = colnames(HSTotal), fixed = T)
UVTotal <- read.delim("F:/ATLAS_Pra/Protein_module/UVTotal2import.txt", header = T)
colnames(UVTotal) <- gsub(pattern = ".", replacement = "-", x = colnames(UVTotal), fixed = T)
HSNucleus <- read.delim("F:/ATLAS_Pra/Protein_module/HSNucleus2import.txt", header = T)
colnames(HSNucleus) <- gsub(pattern = ".", replacement = "-", x = colnames(HSNucleus), fixed = T)
UVNucleus <- read.delim("F:/ATLAS_Pra/Protein_module/UVNucleus2import.txt", header = T)
colnames(UVNucleus) <- gsub(pattern = ".", replacement = "-", x = colnames(UVNucleus), fixed = T)
HSChloroplast <- read.delim("F:/ATLAS_Pra/Protein_module/HSChloroplast2import.txt", header = T)
colnames(HSChloroplast) <- gsub(pattern = ".", replacement = "-", x = colnames(HSChloroplast), fixed = T)
UVChloroplast <- read.delim("F:/ATLAS_Pra/Protein_module/UVChloroplast2import.txt", header = T)
colnames(UVChloroplast) <- gsub(pattern = ".", replacement = "-", x = colnames(UVChloroplast), fixed = T)

ProtList.Tissues <- list(Tissues = Tissues)
class(ProtList.Tissues) <- "pRoDS"
ProtList.Tissues.preprocessed <- preprocess_omic_list(datalist = ProtList.Tissues, initialrow = 1, initialcolumn = 2, treatment1col = 1, treatment2col = NULL, treatment = 1, imputation = "RF", abdbal = "AvgIntensity", varsel = TRUE, parallel = TRUE, varselthld = 0.5,imputhld = 0.34)
write.table(ProtList.Tissues.preprocessed$Tissues, "F:/ATLAS_Pra/Protein_module/TissuesJIC.txt", sep = "\t", quote = F, row.names = F, col.names = T)
rownames(Protein_list.filtered.names$Tissues) <- Protein_list.filtered.names$Tissues$ProtID
Tissue.Preprocessed <- Protein_list.filtered.names$Tissues[colnames(ProtList.Tissues.preprocessed$Tissues)[-1],]
which((rownames(Tissues.Preprocessed) == rownames(t(ProtList.Tissues.preprocessed$Tissues[,-1]))) == FALSE) # check
Tissue.Preprocessed[,2:13] <- t(ProtList.Tissues.preprocessed$Tissues[,-1])[,1:12]
write.table(Tissue.Preprocessed, "F:/ATLAS_Pra/Protein_module/TissuesPreprocessed.txt", sep = "\t", quote = F, row.names = F, col.names = T)

ProtList.FUTotal <- list(FUTotal = FUTotal)
class(ProtList.FUTotal) <- "pRoDS"
ProtList.FUTotal.preprocessed <- preprocess_omic_list(datalist = ProtList.FUTotal, initialrow = 1, initialcolumn = 2, treatment1col = 1, treatment2col = NULL, treatment = 1, imputation = "RF", abdbal = "AvgIntensity", varsel = TRUE, parallel = TRUE, varselthld = 0.5,imputhld = 0.34)
write.table(ProtList.FUTotal.preprocessed$FUTotal, "F:/ATLAS_Pra/Protein_module/FUTotalJIC.txt", sep = "\t", quote = F, row.names = F, col.names = T)
rownames(Protein_list.filtered.names$FU_Total) <- Protein_list.filtered.names$FU_Total$ProtID
FUTotal.Preprocessed <- Protein_list.filtered.names$FU_Total[colnames(ProtList.FUTotal.preprocessed$FUTotal)[-1],]
which((rownames(FUTotal.Preprocessed) == rownames(t(ProtList.FUTotal.preprocessed$FUTotal[,-1]))) == FALSE) # check
FUTotal.Preprocessed[,2:9] <- t(ProtList.FUTotal.preprocessed$FUTotal[,-1])[,1:8]
write.table(FUTotal.Preprocessed, "F:/ATLAS_Pra/Protein_module/FUTotalPreprocessed.txt", sep = "\t", quote = F, row.names = F, col.names = T)

ProtList.HSTotal <- list(HSTotal = HSTotal)
class(ProtList.HSTotal) <- "pRoDS"
ProtList.HSTotal.preprocessed <- preprocess_omic_list(datalist = ProtList.HSTotal, initialrow = 1, initialcolumn = 2, treatment1col = 1, treatment2col = NULL, treatment = 1, imputation = "RF", abdbal = "AvgIntensity", varsel = TRUE, parallel = TRUE, varselthld = 0.5,imputhld = 0.34)
write.table(ProtList.HSTotal.preprocessed$HSTotal, "F:/ATLAS_Pra/Protein_module/HSTotalJIC.txt", sep = "\t", quote = F, row.names = F, col.names = T)
rownames(Protein_list.filtered.names$HS_Total) <- Protein_list.filtered.names$HS_Total$ProtID
HSTotal.Preprocessed <- Protein_list.filtered.names$HS_Total[colnames(ProtList.HSTotal.preprocessed$HSTotal)[-1],]
which((rownames(HSTotal.Preprocessed) == rownames(t(ProtList.HSTotal.preprocessed$HSTotal[,-1]))) == FALSE) # check
HSTotal.Preprocessed[,2:10] <- t(ProtList.HSTotal.preprocessed$HSTotal[,-1])[,1:9]
write.table(HSTotal.Preprocessed, "F:/ATLAS_Pra/Protein_module/HSTotalPreprocessed.txt", sep = "\t", quote = F, row.names = F, col.names = T)

ProtList.UVTotal <- list(UVTotal = UVTotal)
class(ProtList.UVTotal) <- "pRoDS"
ProtList.UVTotal.preprocessed <- preprocess_omic_list(datalist = ProtList.UVTotal, initialrow = 1, initialcolumn = 2, treatment1col = 1, treatment2col = NULL, treatment = 1, imputation = "RF", abdbal = "AvgIntensity", varsel = TRUE, parallel = TRUE, varselthld = 0.5,imputhld = 0.34)
write.table(ProtList.UVTotal.preprocessed$UVTotal, "F:/ATLAS_Pra/Protein_module/UVTotalJIC.txt", sep = "\t", quote = F, row.names = F, col.names = T)
rownames(Protein_list.filtered.names$UV_Total) <- Protein_list.filtered.names$UV_Total$ProtID
UVTotal.Preprocessed <- Protein_list.filtered.names$UV_Total[colnames(ProtList.UVTotal.preprocessed$UVTotal)[-1],]
which((rownames(UVTotal.Preprocessed) == rownames(t(ProtList.UVTotal.preprocessed$UVTotal[,-1]))) == FALSE) # check
UVTotal.Preprocessed[,2:16] <- t(ProtList.UVTotal.preprocessed$UVTotal[,-1])[,1:15]
write.table(UVTotal.Preprocessed, "F:/ATLAS_Pra/Protein_module/UVTotalPreprocessed.txt", sep = "\t", quote = F, row.names = F, col.names = T)

ProtList.HSNucleus <- list(HSNucleus = HSNucleus)
class(ProtList.HSNucleus) <- "pRoDS"
ProtList.HSNucleus.preprocessed <- preprocess_omic_list(datalist = ProtList.HSNucleus, initialrow = 1, initialcolumn = 2, treatment1col = 1, treatment2col = NULL, treatment = 1, imputation = "RF", abdbal = "AvgIntensity", varsel = TRUE, parallel = TRUE, varselthld = 0.5,imputhld = 0.34)
write.table(ProtList.HSNucleus.preprocessed$HSNucleus, "F:/ATLAS_Pra/Protein_module/HSNucleusJIC.txt", sep = "\t", quote = F, row.names = F, col.names = T)
rownames(Protein_list.filtered.names$HS_Nucleus) <- Protein_list.filtered.names$HS_Nucleus$ProtID
HSNucleus.Preprocessed <- Protein_list.filtered.names$HS_Nucleus[colnames(ProtList.HSNucleus.preprocessed$HSNucleus)[-1],]
which((rownames(HSNucleus.Preprocessed) == rownames(t(ProtList.HSNucleus.preprocessed$HSNucleus[,-1]))) == FALSE) # check
HSNucleus.Preprocessed[,2:29] <- t(ProtList.HSNucleus.preprocessed$HSNucleus[,-1])[,1:28]
write.table(HSNucleus.Preprocessed, "F:/ATLAS_Pra/Protein_module/HSNucleusPreprocessed.txt", sep = "\t", quote = F, row.names = F, col.names = T)

ProtList.UVNucleus <- list(UVNucleus = UVNucleus)
class(ProtList.UVNucleus) <- "pRoDS"
ProtList.UVNucleus.preprocessed <- preprocess_omic_list(datalist = ProtList.UVNucleus, initialrow = 1, initialcolumn = 2, treatment1col = 1, treatment2col = NULL, treatment = 1, imputation = "RF", abdbal = "AvgIntensity", varsel = TRUE, parallel = TRUE, varselthld = 0.5,imputhld = 0.34)
write.table(ProtList.UVNucleus.preprocessed$UVNucleus, "F:/ATLAS_Pra/Protein_module/UVNucleusJIC.txt", sep = "\t", quote = F, row.names = F, col.names = T)
rownames(Protein_list.filtered.names$UV_Nucleus) <- Protein_list.filtered.names$UV_Nucleus$ProtID
UVNucleus.Preprocessed <- Protein_list.filtered.names$UV_Nucleus[colnames(ProtList.UVNucleus.preprocessed$UVNucleus)[-1],]
which((rownames(UVNucleus.Preprocessed) == rownames(t(ProtList.UVNucleus.preprocessed$UVNucleus[,-1]))) == FALSE) # check
UVNucleus.Preprocessed[,2:12] <- t(ProtList.UVNucleus.preprocessed$UVNucleus[,-1])[,1:11]
write.table(UVNucleus.Preprocessed, "F:/ATLAS_Pra/Protein_module/UVNucleusPreprocessed.txt", sep = "\t", quote = F, row.names = F, col.names = T)

ProtList.HSChloroplast <- list(HSChloroplast = HSChloroplast)
class(ProtList.HSChloroplast) <- "pRoDS"
ProtList.HSChloroplast.preprocessed <- preprocess_omic_list(datalist = ProtList.HSChloroplast, initialrow = 1, initialcolumn = 2, treatment1col = 1, treatment2col = NULL, treatment = 1, imputation = "RF", abdbal = "AvgIntensity", varsel = TRUE, parallel = TRUE, varselthld = 0.5,imputhld = 0.34)
write.table(ProtList.HSChloroplast.preprocessed$HSChloroplast, "F:/ATLAS_Pra/Protein_module/HSChloroplastJIC.txt", sep = "\t", quote = F, row.names = F, col.names = T)
rownames(Protein_list.filtered.names$HS_Chloroplast) <- Protein_list.filtered.names$HS_Chloroplast$ProtID
HSChloroplast.Preprocessed <- Protein_list.filtered.names$HS_Chloroplast[colnames(ProtList.HSChloroplast.preprocessed$HSChloroplast)[-1],]
which((rownames(HSChloroplast.Preprocessed) == rownames(t(ProtList.HSChloroplast.preprocessed$HSChloroplast[,-1]))) == FALSE) # check
HSChloroplast.Preprocessed[,2:33] <- t(ProtList.HSChloroplast.preprocessed$HSChloroplast[,-1])[,1:32]
write.table(HSChloroplast.Preprocessed, "F:/ATLAS_Pra/Protein_module/HSChloroplastPreprocessed.txt", sep = "\t", quote = F, row.names = F, col.names = T)

ProtList.UVChloroplast <- list(UVChloroplast = UVChloroplast)
class(ProtList.UVChloroplast) <- "pRoDS"
ProtList.UVChloroplast.preprocessed <- preprocess_omic_list(datalist = ProtList.UVChloroplast, initialrow = 1, initialcolumn = 2, treatment1col = 1, treatment2col = NULL, treatment = 1, imputation = "RF", abdbal = "AvgIntensity", varsel = TRUE, parallel = TRUE, varselthld = 0.5,imputhld = 0.34)
write.table(ProtList.UVChloroplast.preprocessed$UVChloroplast, "F:/ATLAS_Pra/Protein_module/UVChloroplastJIC.txt", sep = "\t", quote = F, row.names = F, col.names = T)
rownames(Protein_list.filtered.names$UV_Chloroplast) <- Protein_list.filtered.names$UV_Chloroplast$ProtID
UVChloroplast.Preprocessed <- Protein_list.filtered.names$UV_Chloroplast[colnames(ProtList.UVChloroplast.preprocessed$UVChloroplast)[-1],]
which((rownames(UVChloroplast.Preprocessed) == rownames(t(ProtList.UVChloroplast.preprocessed$UVChloroplast[,-1]))) == FALSE) # check
UVChloroplast.Preprocessed[,2:33] <- t(ProtList.UVChloroplast.preprocessed$UVChloroplast[,-1])[,1:32]
write.table(UVChloroplast.Preprocessed, "F:/ATLAS_Pra/Protein_module/UVChloroplastPreprocessed.txt", sep = "\t", quote = F, row.names = F, col.names = T)

#### 1.4 Imputting replicates ####
## FU Total - impute two replicates x treatment == 12 total samples

View(ProtList.FUTotal.preprocessed$FUTotal)
ProtList.FUTotal.preprocessed.C <- ProtList.FUTotal.preprocessed
ProtList.FUTotal.preprocessed.C$FUTotal <- ProtList.FUTotal.preprocessed.C$FUTotal[1:4,]
ProtList.FUTotal.preprocessed.C$FUTotal[5,] <- rep(NA, ncol(ProtList.FUTotal.preprocessed.C$FUTotal))
ProtList.FUTotal.preprocessed.C$FUTotal[5,1] <- "C"
allmisscols <- which((sapply(ProtList.FUTotal.preprocessed.C$FUTotal, function(x) all(is.na(x) | x == 0 ))) == TRUE)
ProtList.FUTotal.preprocessed.C1 <- ProtList.FUTotal.preprocessed.C
ProtList.FUTotal.preprocessed.C2 <- ProtList.FUTotal.preprocessed.C
ProtList.FUTotal.preprocessed.C1.imputed <- preprocess_omic_list(datalist = ProtList.FUTotal.preprocessed.C1, initialrow = 1, initialcolumn = 2, treatment1col = 1, treatment2col = NULL, treatment = 1, imputation = "RF", abdbal = "none", varsel = FALSE, parallel = TRUE, imputhld = 0.34)
ProtList.FUTotal.preprocessed.C2.imputed <- preprocess_omic_list(datalist = ProtList.FUTotal.preprocessed.C2, initialrow = 1, initialcolumn = 2, treatment1col = 1, treatment2col = NULL, treatment = 1, imputation = "RF", abdbal = "none", varsel = FALSE, parallel = TRUE, imputhld = 0.34)
ProtList.FUTotal.preprocessed.C.imputed <- ProtList.FUTotal.preprocessed.C1.imputed
ProtList.FUTotal.preprocessed.C.imputed$FUTotal <- rbind(ProtList.FUTotal.preprocessed.C.imputed$FUTotal, ProtList.FUTotal.preprocessed.C2.imputed$FUTotal[5,])
allmisscols
MissCols <- data.frame("all-04-525948" = rep(0, 6),
                       "all-04-346386" = rep(0, 6),
                       "all-04-2104640" = rep(0, 6),
                       "all-04-1585198" = rep(0, 6))
colnames(MissCols) <-  gsub(".", "-", colnames(MissCols), fixed = T)
allmisscols
ProtList.FUTotal.preprocessed.C.imputed$FUTotal <- cbind(ProtList.FUTotal.preprocessed.C.imputed$FUTotal[,1:470], "all-04-525948" = MissCols[,1], 
      ProtList.FUTotal.preprocessed.C.imputed$FUTotal[,471:517],  "all-04-346386"= MissCols[,2],
      ProtList.FUTotal.preprocessed.C.imputed$FUTotal[,518:553], "all-04-2104640" = MissCols[,3],
      ProtList.FUTotal.preprocessed.C.imputed$FUTotal[,554:1092], "all-04-1585198" = MissCols[,4],
      ProtList.FUTotal.preprocessed.C.imputed$FUTotal[,1093:ncol(ProtList.FUTotal.preprocessed.C.imputed$FUTotal)])
which((colnames(ProtList.FUTotal.preprocessed$FUTotal) == colnames(ProtList.FUTotal.preprocessed.C.imputed$FUTotal)) == FALSE) # check
write.table(ProtList.FUTotal.preprocessed.C.imputed$FUTotal, "F:/ATLAS_Pra/Protein_module/0.Process/3.Imputed/FUTotal.Cimputed.txt", sep = "\t", quote = F, col.names = T)

ProtList.FUTotal.preprocessed.FU <- ProtList.FUTotal.preprocessed
ProtList.FUTotal.preprocessed.FU$FUTotal <- ProtList.FUTotal.preprocessed.FU$FUTotal[5:8,]
ProtList.FUTotal.preprocessed.FU$FUTotal[5,] <- rep(NA, ncol(ProtList.FUTotal.preprocessed.FU$FUTotal))
ProtList.FUTotal.preprocessed.FU$FUTotal[5,1] <- "FU"
allmisscols.FU <- which((sapply(ProtList.FUTotal.preprocessed.FU$FUTotal, function(x) all(is.na(x) | x == 0 ))) == TRUE)
ProtList.FUTotal.preprocessed.FU1 <- ProtList.FUTotal.preprocessed.FU
ProtList.FUTotal.preprocessed.FU2 <- ProtList.FUTotal.preprocessed.FU
ProtList.FUTotal.preprocessed.FU1.imputed <- preprocess_omic_list(datalist = ProtList.FUTotal.preprocessed.FU1, initialrow = 1, initialcolumn = 2, treatment1col = 1, treatment2col = NULL, treatment = 1, imputation = "RF", abdbal = "none", varsel = FALSE, parallel = TRUE, imputhld = 0.34)
ProtList.FUTotal.preprocessed.FU2.imputed <- preprocess_omic_list(datalist = ProtList.FUTotal.preprocessed.FU2, initialrow = 1, initialcolumn = 2, treatment1col = 1, treatment2col = NULL, treatment = 1, imputation = "RF", abdbal = "none", varsel = FALSE, parallel = TRUE, imputhld = 0.34)
ProtList.FUTotal.preprocessed.FU.imputed <- ProtList.FUTotal.preprocessed.FU1.imputed
#ProtList.FUTotal.preprocessed.FU.imputed$FUTotal <- rbind(ProtList.FUTotal.preprocessed.FU.imputed$FUTotal, ProtList.FUTotal.preprocessed.FU2.imputed$FUTotal[5,])
allmisscols.FU
MissCols <- data.frame("all-04-1356560" = rep(0, 5),
                       "all-04-852530" = rep(0, 5),
                       "all-04-328985" = rep(0, 5),
                       "all-04-1066968" = rep(0, 5),
                       "all-04-2103340" = rep(0, 5))
colnames(MissCols) <-  gsub(".", "-", colnames(MissCols), fixed = T)
allmisscols.FU
ProtList.FUTotal.preprocessed.FU.imputed$FUTotal <- cbind(ProtList.FUTotal.preprocessed.FU.imputed$FUTotal[,1:927], "all-04-1356560" = MissCols[,1], 
                                                         ProtList.FUTotal.preprocessed.FU.imputed$FUTotal[,928:966],  "all-04-852530"= MissCols[,2],
                                                         ProtList.FUTotal.preprocessed.FU.imputed$FUTotal[,967:1422], "all-04-328985" = MissCols[,3],
                                                         ProtList.FUTotal.preprocessed.FU.imputed$FUTotal[,1423:1438], "all-04-1066968" = MissCols[,4],
                                                         ProtList.FUTotal.preprocessed.FU.imputed$FUTotal[,1439:1473], "all-04-2103340" = MissCols[,5],
                                                         ProtList.FUTotal.preprocessed.FU.imputed$FUTotal[,1474:ncol(ProtList.FUTotal.preprocessed.FU.imputed$FUTotal)])
which((colnames(ProtList.FUTotal.preprocessed$FUTotal) == colnames(ProtList.FUTotal.preprocessed.FU.imputed$FUTotal)) == FALSE) # check
write.table(ProtList.FUTotal.preprocessed.FU.imputed$FUTotal, "F:/ATLAS_Pra/Protein_module/0.Process/3.Imputed/FUTotal.FUimputed.txt", sep = "\t", quote = F, col.names = T)

# final table fusion as in 1.3
FUTotal.Preprocessed.Imputed <- rbind(ProtList.FUTotal.preprocessed.C.imputed$FUTotal, ProtList.FUTotal.preprocessed.FU.imputed$FUTotal)
which((rownames(FUTotal.Preprocessed) == rownames(t(FUTotal.Preprocessed.Imputed[,-1]))) == FALSE) # check
FUTotal.Preprocessed.Imputed.F <- FUTotal.Preprocessed 
FUTotal.Preprocessed.Imputed.F$FUTotal_C_5 <- t(FUTotal.Preprocessed.Imputed[,-1])[,5]
FUTotal.Preprocessed.Imputed.F$FUTotal_C_6 <- t(FUTotal.Preprocessed.Imputed[,-1])[,6]
FUTotal.Preprocessed.Imputed.F$FUTotal_FU_5 <- t(FUTotal.Preprocessed.Imputed[,-1])[,11]
FUTotal.Preprocessed.Imputed.F <- FUTotal.Preprocessed.Imputed.F[,c(1:5,11:12,6:9,13,10)]
write.table(FUTotal.Preprocessed.Imputed.F, "F:/ATLAS_Pra/Protein_module/0.Process/3.Imputed/FUTotal.Preprocessed.Imputed.txt", sep = "\t", quote = F, row.names = F, col.names = T)

## HS Total - impute one replicate x treatment == 12 total samples

View(ProtList.HSTotal.preprocessed$HSTotal)
ProtList.HSTotal.preprocessed.C <- ProtList.HSTotal.preprocessed
ProtList.HSTotal.preprocessed.C$HSTotal <- ProtList.HSTotal.preprocessed.C$HSTotal[1:3,]
ProtList.HSTotal.preprocessed.C$HSTotal[4,] <- rep(NA, ncol(ProtList.HSTotal.preprocessed.C$HSTotal))
ProtList.HSTotal.preprocessed.C$HSTotal[4,1] <- "C"
allmisscols.HS.C <- which((sapply(ProtList.HSTotal.preprocessed.C$HSTotal, function(x) all(is.na(x) | x == 0 ))) == TRUE)
ProtList.HSTotal.preprocessed.C.imputed <- preprocess_omic_list(datalist = ProtList.HSTotal.preprocessed.C, initialrow = 1, initialcolumn = 2, treatment1col = 1, treatment2col = NULL, treatment = 1, imputation = "RF", abdbal = "none", varsel = FALSE, parallel = TRUE, imputhld = 0.34)
allmisscols.HS.C
MissCols <- data.frame(matrix(data = 0, nrow = 4, ncol = length(allmisscols.HS.C)))
colnames(MissCols) <- names(allmisscols.HS.C)
allmisscols.HS.C
ProtList.HSTotal.preprocessed.C.imputed$HSTotal <- cbind(ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,1:65], "all-04-1546141" = MissCols[,1],  "all-04-20870" = MissCols[,2], 
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,66:156],  "all-04-70791"= MissCols[,3],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,157:201], "all-04-511017" = MissCols[,4],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,202:216], "all-04-1361101" = MissCols[,5],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,217:253], "all-04-2368571" = MissCols[,6],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,254:320], "all-04-1249331" = MissCols[,7],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,321:343], "all-04-1487771" = MissCols[,8],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,344:376], "all-04-1531421" = MissCols[,9],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,377:388], "all-04-1363102" = MissCols[,10],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,389:(461-11)], "all-04-1029550" = MissCols[,11],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,451], "all-04-1012791" = MissCols[,12],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,(463-12+1):(470-13)], "all-04-1555365" = MissCols[,13],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,(470-13+1):(495-14)], "all-04-1363928" = MissCols[,14],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,(495-14+1):(517-15)], "all-04-344628" = MissCols[,15],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,503], "all-04-1380597" = MissCols[,16],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,(519-16+1):(541-17)], "all-04-1623427" = MissCols[,17],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,(541-17+1):(563-18)], "all-04-1945198" = MissCols[,18],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,(563-18+1):(603-19)], "all-04-1973711" = MissCols[,19],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,(603-19+1):(795-20)], "all-04-1023804" = MissCols[,20],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,(795-20+1):(843-21)], "all-04-1361235" = MissCols[,21], "all-04-492999" = MissCols[,22],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,(843-21+1):(849-23)], "all-04-1364366" = MissCols[,23],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,(849-23+1):(895-24)], "all-04-1039775" = MissCols[,24],
                                                         ProtList.HSTotal.preprocessed.C.imputed$HSTotal[,(895-24+1):ncol(ProtList.HSTotal.preprocessed.C.imputed$HSTotal)])

which((colnames(ProtList.HSTotal.preprocessed$HSTotal) == colnames(ProtList.HSTotal.preprocessed.C.imputed$HSTotal)) == FALSE) # check
colnames(ProtList.HSTotal.preprocessed.C.imputed$HSTotal)[462] <- colnames(ProtList.HSTotal.preprocessed$HSTotal)[462]
colnames(ProtList.HSTotal.preprocessed.C.imputed$HSTotal)[518] <- colnames(ProtList.HSTotal.preprocessed$HSTotal)[518]
which((colnames(ProtList.HSTotal.preprocessed$HSTotal) == colnames(ProtList.HSTotal.preprocessed.C.imputed$HSTotal)) == FALSE) # check
write.table(ProtList.HSTotal.preprocessed.C.imputed$HSTotal, "F:/ATLAS_Pra/Protein_module/0.Process/3.Imputed/HSTotal.Cimputed.txt", sep = "\t", quote = F, col.names = T)

View(ProtList.HSTotal.preprocessed$HSTotal)
ProtList.HSTotal.preprocessed.T1 <- ProtList.HSTotal.preprocessed
ProtList.HSTotal.preprocessed.T1$HSTotal <- ProtList.HSTotal.preprocessed.T1$HSTotal[4:6,]
ProtList.HSTotal.preprocessed.T1$HSTotal[4,] <- rep(NA, ncol(ProtList.HSTotal.preprocessed.T1$HSTotal))
ProtList.HSTotal.preprocessed.T1$HSTotal[4,1] <- "T1"
allmisscols.HS.T1 <- which((sapply(ProtList.HSTotal.preprocessed.T1$HSTotal, function(x) all(is.na(x) | x == 0 ))) == TRUE)
ProtList.HSTotal.preprocessed.T1.imputed <- preprocess_omic_list(datalist = ProtList.HSTotal.preprocessed.T1, initialrow = 1, initialcolumn = 2, treatment1col = 1, treatment2col = NULL, treatment = 1, imputation = "RF", abdbal = "none", varsel = FALSE, parallel = TRUE, imputhld = 0.34)
allmisscols.HS.T1
MissCols <- data.frame(matrix(data = 0, nrow = 4, ncol = length(allmisscols.HS.T1)))
colnames(MissCols) <- names(allmisscols.HS.T1)
allmisscols.HS.T1
ProtList.HSTotal.preprocessed.T1.imputed$HSTotal <- cbind(ProtList.HSTotal.preprocessed.T1.imputed$HSTotal[,1:66], "all-04-20870" = MissCols[,1], 
                                                         ProtList.HSTotal.preprocessed.T1.imputed$HSTotal[,(66+1):(405-2)],  "all-04-1549211"= MissCols[,2],
                                                         ProtList.HSTotal.preprocessed.T1.imputed$HSTotal[,(405-2+1):(470-3)], "all-04-1555365" = MissCols[,3],
                                                         ProtList.HSTotal.preprocessed.T1.imputed$HSTotal[,(470-3+1):(519-4)], "all-04-1380597" = MissCols[,4],
                                                         ProtList.HSTotal.preprocessed.T1.imputed$HSTotal[,(519-4+1):(706-5)], "all-04-1591626" = MissCols[,5],
                                                         ProtList.HSTotal.preprocessed.T1.imputed$HSTotal[,(706-5+1):ncol(ProtList.HSTotal.preprocessed.T1.imputed$HSTotal)])

which((colnames(ProtList.HSTotal.preprocessed$HSTotal) == colnames(ProtList.HSTotal.preprocessed.T1.imputed$HSTotal)) == FALSE) # check
write.table(ProtList.HSTotal.preprocessed.T1.imputed$HSTotal, "F:/ATLAS_Pra/Protein_module/0.Process/3.Imputed/HSTotal.T1imputed.txt", sep = "\t", quote = F, col.names = T)

View(ProtList.HSTotal.preprocessed$HSTotal)
ProtList.HSTotal.preprocessed.T3 <- ProtList.HSTotal.preprocessed
ProtList.HSTotal.preprocessed.T3$HSTotal <- ProtList.HSTotal.preprocessed.T3$HSTotal[7:9,]
ProtList.HSTotal.preprocessed.T3$HSTotal[4,] <- rep(NA, ncol(ProtList.HSTotal.preprocessed.T3$HSTotal))
ProtList.HSTotal.preprocessed.T3$HSTotal[4,1] <- "T3"
allmisscols.HS.T3 <- which((sapply(ProtList.HSTotal.preprocessed.T3$HSTotal, function(x) all(is.na(x) | x == 0 ))) == TRUE)
ProtList.HSTotal.preprocessed.T3.imputed <- preprocess_omic_list(datalist = ProtList.HSTotal.preprocessed.T3, initialrow = 1, initialcolumn = 2, treatment1col = 1, treatment2col = NULL, treatment = 1, imputation = "RF", abdbal = "none", varsel = FALSE, parallel = TRUE, imputhld = 0.34)
allmisscols.HS.T3
MissCols <- data.frame(matrix(data = 0, nrow = 4, ncol = length(allmisscols.HS.T3)))
colnames(MissCols) <- names(allmisscols.HS.T3)
allmisscols.HS.T3
ProtList.HSTotal.preprocessed.T3.imputed$HSTotal <- cbind(ProtList.HSTotal.preprocessed.T3.imputed$HSTotal[,1:193], "all-04-1038509" = MissCols[,1], 
                                                          ProtList.HSTotal.preprocessed.T3.imputed$HSTotal[,(193+1):(222-2)],  "all-04-1599001"= MissCols[,2],
                                                          ProtList.HSTotal.preprocessed.T3.imputed$HSTotal[,(222-2+1):(322-3)], "all-04-1955465" = MissCols[,3],
                                                          ProtList.HSTotal.preprocessed.T3.imputed$HSTotal[,(322-3+1):(348-4)], "all-04-876315" = MissCols[,4],
                                                          ProtList.HSTotal.preprocessed.T3.imputed$HSTotal[,(348-4+1):(463-5)], "all-04-1012791" = MissCols[,5],
                                                          ProtList.HSTotal.preprocessed.T3.imputed$HSTotal[,(463-5+1):(477-6)],  "all-04-1023028"= MissCols[,6],
                                                          ProtList.HSTotal.preprocessed.T3.imputed$HSTotal[,(477-6+1):(526-7)], "all-04-2124392" = MissCols[,7],
                                                          ProtList.HSTotal.preprocessed.T3.imputed$HSTotal[,(526-7+1):(575-8)], "all-04-2333530" = MissCols[,8],
                                                          ProtList.HSTotal.preprocessed.T3.imputed$HSTotal[,(575-8+1):(662-9)], "all-04-1363635" = MissCols[,9],
                                                          ProtList.HSTotal.preprocessed.T3.imputed$HSTotal[,(662-9+1):(818-10)],  "all-04-1207341"= MissCols[,10],
                                                          ProtList.HSTotal.preprocessed.T3.imputed$HSTotal[,(818-10+1):(843-11)], "all-04-1361235" = MissCols[,11],
                                                          ProtList.HSTotal.preprocessed.T3.imputed$HSTotal[,(843-11+1):(864-12)], "all-04-1933409" = MissCols[,12],
                                                          ProtList.HSTotal.preprocessed.T3.imputed$HSTotal[,(864-12+1):ncol(ProtList.HSTotal.preprocessed.T3.imputed$HSTotal)])

which((colnames(ProtList.HSTotal.preprocessed$HSTotal) == colnames(ProtList.HSTotal.preprocessed.T3.imputed$HSTotal)) == FALSE) # check
write.table(ProtList.HSTotal.preprocessed.T3.imputed$HSTotal, "F:/ATLAS_Pra/Protein_module/0.Process/3.Imputed/HSTotal.T3imputed.txt", sep = "\t", quote = F, col.names = T)


# final table fusion as in 1.3
HSTotal.Preprocessed.Imputed <- rbind(ProtList.HSTotal.preprocessed.C.imputed$HSTotal, ProtList.HSTotal.preprocessed.T1.imputed$HSTotal, ProtList.HSTotal.preprocessed.T3.imputed$HSTotal)
which((rownames(HSTotal.Preprocessed) == rownames(t(HSTotal.Preprocessed.Imputed[,-1]))) == FALSE) # check
HSTotal.Preprocessed.Imputed.F <- HSTotal.Preprocessed 
HSTotal.Preprocessed.Imputed.F$HSTotal_C_4 <- t(HSTotal.Preprocessed.Imputed[,-1])[,4]
HSTotal.Preprocessed.Imputed.F$HSTotal_T1_4 <- t(HSTotal.Preprocessed.Imputed[,-1])[,8]
HSTotal.Preprocessed.Imputed.F$HSTotal_T3_4 <- t(HSTotal.Preprocessed.Imputed[,-1])[,12]
HSTotal.Preprocessed.Imputed.F <- HSTotal.Preprocessed.Imputed.F[,c(1:4,12,5:7,13,8:10,14,11)]
write.table(HSTotal.Preprocessed.Imputed.F, "F:/ATLAS_Pra/Protein_module/0.Process/3.Imputed/HSTotal.Preprocessed.Imputed.txt", sep = "\t", quote = F, row.names = F, col.names = T)

## UV Nucleus - impute one replicate in T3 treatment == 12 total samples

View(ProtList.UVNucleus.preprocessed$UVNucleus)
ProtList.UVNucleus.preprocessed.T3 <- ProtList.UVNucleus.preprocessed
ProtList.UVNucleus.preprocessed.T3$UVNucleus <- ProtList.UVNucleus.preprocessed.T3$UVNucleus[9:11,]
ProtList.UVNucleus.preprocessed.T3$UVNucleus[4,] <- rep(NA, ncol(ProtList.UVNucleus.preprocessed.T3$UVNucleus))
ProtList.UVNucleus.preprocessed.T3$UVNucleus[4,1] <- "T3"
allmisscols.UVNucleus.T3 <- which((sapply(ProtList.UVNucleus.preprocessed.T3$UVNucleus, function(x) all(is.na(x) | x == 0 ))) == TRUE)
ProtList.UVNucleus.preprocessed.T3.imputed <- preprocess_omic_list(datalist = ProtList.UVNucleus.preprocessed.T3, initialrow = 1, initialcolumn = 2, treatment1col = 1, treatment2col = NULL, treatment = 1, imputation = "RF", abdbal = "none", varsel = FALSE, parallel = TRUE, imputhld = 0.34)
allmisscols.UVNucleus.T3
MissCols <- data.frame(matrix(data = 0, nrow = 4, ncol = length(allmisscols.UVNucleus.T3)))
colnames(MissCols) <- names(allmisscols.UVNucleus.T3)
allmisscols.UVNucleus.T3
ProtList.UVNucleus.preprocessed.T3.imputed$UVNucleus <- cbind(ProtList.UVNucleus.preprocessed.T3.imputed$UVNucleus[,1:170], "all-04-904959" = MissCols[,1], 
                                                         ProtList.UVNucleus.preprocessed.T3.imputed$UVNucleus[,171:(345-2)],  "all-04-1552990"= MissCols[,2],
                                                         ProtList.UVNucleus.preprocessed.T3.imputed$UVNucleus[,(345-2+1):ncol(ProtList.UVNucleus.preprocessed.T3.imputed$UVNucleus)])

which((colnames(ProtList.UVNucleus.preprocessed$UVNucleus) == colnames(ProtList.UVNucleus.preprocessed.T3.imputed$UVNucleus)) == FALSE) # check
write.table(ProtList.UVNucleus.preprocessed.T3.imputed$UVNucleus, "F:/ATLAS_Pra/Protein_module/0.Process/3.Imputed/UVNucleus.T3imputed.txt", sep = "\t", quote = F, col.names = T)

# final table fusion as in 1.3
UVNucleus.Preprocessed.Imputed <- rbind(ProtList.UVNucleus.preprocessed$UVNucleus[1:8,], ProtList.UVNucleus.preprocessed.T3.imputed$UVNucleus)
which((rownames(UVNucleus.Preprocessed) == rownames(t(UVNucleus.Preprocessed.Imputed[,-1]))) == FALSE) # check
UVNucleus.Preprocessed.Imputed.F <- UVNucleus.Preprocessed
UVNucleus.Preprocessed.Imputed.F$UVNucleus_T3_4 <- t(UVNucleus.Preprocessed.Imputed[,-1])[,12]
UVNucleus.Preprocessed.Imputed.F <- UVNucleus.Preprocessed.Imputed.F[,c(1:12,14,13)]
write.table(UVNucleus.Preprocessed.Imputed.F, "F:/ATLAS_Pra/Protein_module/0.Process/3.Imputed/UVNucleus.Preprocessed.Imputed.txt", sep = "\t", quote = F, row.names = F, col.names = T)

#### 1.5 Merge all proteomes ####

ProtID <- c(Tissue.Preprocessed$ProtID, FUTotal.Preprocessed.Imputed.F$ProtID, HSTotal.Preprocessed.Imputed.F$ProtID,
  UVTotal.Preprocessed$ProtID, HSNucleus.Preprocessed$ProtID, UVNucleus.Preprocessed.Imputed.F$ProtID,
  HSChloroplast.Preprocessed$ProtID, UVChloroplast.Preprocessed$ProtID)

TranscriptID <- c(Tissue.Preprocessed$TranscriptID, FUTotal.Preprocessed.Imputed.F$TranscriptID, HSTotal.Preprocessed.Imputed.F$TranscriptID,
            UVTotal.Preprocessed$TranscriptID, HSNucleus.Preprocessed$TranscriptID, UVNucleus.Preprocessed.Imputed.F$TranscriptID,
            HSChloroplast.Preprocessed$TranscriptID, UVChloroplast.Preprocessed$TranscriptID)

ProtAtlas <- data.frame(ProtID = ProtID, TranscriptID = TranscriptID)
ProtAtlas <- unique(ProtAtlas)
ProtAtlas <- ProtAtlas[order(ProtAtlas$ProtID),]

Tissue.Preprocessed <- Tissue.Preprocessed[order(Tissue.Preprocessed$ProtID),]
FUTotal.Preprocessed.Imputed.F <- FUTotal.Preprocessed.Imputed.F[order(FUTotal.Preprocessed.Imputed.F$ProtID),]
HSTotal.Preprocessed.Imputed.F <- HSTotal.Preprocessed.Imputed.F[order(HSTotal.Preprocessed.Imputed.F$ProtID),]
UVTotal.Preprocessed <- UVTotal.Preprocessed[order(UVTotal.Preprocessed$ProtID),]
HSNucleus.Preprocessed <- HSNucleus.Preprocessed[order(HSNucleus.Preprocessed$ProtID),]
UVNucleus.Preprocessed.Imputed.F <- UVNucleus.Preprocessed.Imputed.F[order(UVNucleus.Preprocessed.Imputed.F$ProtID),]
HSChloroplast.Preprocessed <- HSChloroplast.Preprocessed[order(HSChloroplast.Preprocessed$ProtID),]
UVChloroplast.Preprocessed <- UVChloroplast.Preprocessed[order(UVChloroplast.Preprocessed$ProtID),]

# all columns to 0s
ProtAtlas$Tissues_AA_1 <- 0
ProtAtlas$Tissues_AA_2 <- 0
ProtAtlas$Tissues_AA_3 <- 0
ProtAtlas$Tissues_AJ_1 <- 0
ProtAtlas$Tissues_AJ_2 <- 0
ProtAtlas$Tissues_AJ_3 <- 0
ProtAtlas$Tissues_RJ_1 <- 0
ProtAtlas$Tissues_RJ_2 <- 0
ProtAtlas$Tissues_RJ_3 <- 0
ProtAtlas$Tissues_YF_1 <- 0
ProtAtlas$Tissues_YF_2 <- 0
ProtAtlas$Tissues_YF_3 <- 0

ProtAtlas$FUTotal_C_1 <- 0
ProtAtlas$FUTotal_C_2 <- 0
ProtAtlas$FUTotal_C_3 <- 0
ProtAtlas$FUTotal_C_4 <- 0
ProtAtlas$FUTotal_C_5 <- 0
ProtAtlas$FUTotal_FU_1 <- 0
ProtAtlas$FUTotal_FU_2 <- 0
ProtAtlas$FUTotal_FU_3 <- 0
ProtAtlas$FUTotal_FU_4 <- 0
ProtAtlas$FUTotal_FU_5 <- 0

ProtAtlas$HSTotal_C_1 <- 0
ProtAtlas$HSTotal_C_2 <- 0
ProtAtlas$HSTotal_C_3 <- 0
ProtAtlas$HSTotal_C_4 <- 0
ProtAtlas$HSTotal_T1_1 <- 0
ProtAtlas$HSTotal_T1_2 <- 0
ProtAtlas$HSTotal_T1_3 <- 0
ProtAtlas$HSTotal_T1_4 <- 0
ProtAtlas$HSTotal_T3_1 <- 0
ProtAtlas$HSTotal_T3_2 <- 0
ProtAtlas$HSTotal_T3_3 <- 0
ProtAtlas$HSTotal_T3_4 <- 0

ProtAtlas$UVTotal_C_1 <- 0
ProtAtlas$UVTotal_C_2 <- 0
ProtAtlas$UVTotal_C_3 <- 0
ProtAtlas$UVTotal_T1_1 <- 0
ProtAtlas$UVTotal_T1_2 <- 0
ProtAtlas$UVTotal_T1_3 <- 0
ProtAtlas$UVTotal_T2_1 <- 0
ProtAtlas$UVTotal_T2_2 <- 0
ProtAtlas$UVTotal_T2_3 <- 0
ProtAtlas$UVTotal_T3_1 <- 0
ProtAtlas$UVTotal_T3_2 <- 0
ProtAtlas$UVTotal_T3_3 <- 0
ProtAtlas$UVTotal_R_1 <- 0
ProtAtlas$UVTotal_R_2 <- 0
ProtAtlas$UVTotal_R_3 <- 0

ProtAtlas$HSNucleus_C_1 <- 0
ProtAtlas$HSNucleus_C_2 <- 0
ProtAtlas$HSNucleus_C_3 <- 0
ProtAtlas$HSNucleus_C_4 <- 0
ProtAtlas$HSNucleus_T1_1 <- 0
ProtAtlas$HSNucleus_T1_2 <- 0
ProtAtlas$HSNucleus_T1_3 <- 0
ProtAtlas$HSNucleus_T1_4 <- 0
ProtAtlas$HSNucleus_T3_1 <- 0
ProtAtlas$HSNucleus_T3_2 <- 0
ProtAtlas$HSNucleus_T3_3 <- 0
ProtAtlas$HSNucleus_T3_4 <- 0
ProtAtlas$HSNucleus_T5_1 <- 0
ProtAtlas$HSNucleus_T5_2 <- 0
ProtAtlas$HSNucleus_T5_3 <- 0
ProtAtlas$HSNucleus_T5_4 <- 0
ProtAtlas$HSNucleus_T10_1 <- 0
ProtAtlas$HSNucleus_T10_2 <- 0
ProtAtlas$HSNucleus_T10_3 <- 0
ProtAtlas$HSNucleus_T10_4 <- 0
ProtAtlas$HSNucleus_T5R_1 <- 0
ProtAtlas$HSNucleus_T5R_2 <- 0
ProtAtlas$HSNucleus_T5R_3 <- 0
ProtAtlas$HSNucleus_T5R_4 <- 0
ProtAtlas$HSNucleus_CR_1 <- 0
ProtAtlas$HSNucleus_CR_2 <- 0
ProtAtlas$HSNucleus_CR_3 <- 0
ProtAtlas$HSNucleus_CR_4 <- 0

ProtAtlas$UVNucleus_C_1 <- 0
ProtAtlas$UVNucleus_C_2 <- 0
ProtAtlas$UVNucleus_C_3 <- 0
ProtAtlas$UVNucleus_C_4 <- 0
ProtAtlas$UVNucleus_T1_1 <- 0
ProtAtlas$UVNucleus_T1_2 <- 0
ProtAtlas$UVNucleus_T1_3 <- 0
ProtAtlas$UVNucleus_T1_4 <- 0
ProtAtlas$UVNucleus_T3_1 <- 0
ProtAtlas$UVNucleus_T3_2 <- 0
ProtAtlas$UVNucleus_T3_3 <- 0
ProtAtlas$UVNucleus_T3_4 <- 0

ProtAtlas$HSChloroplast_CE_1 <- 0
ProtAtlas$HSChloroplast_CE_2 <- 0
ProtAtlas$HSChloroplast_CE_3 <- 0
ProtAtlas$HSChloroplast_CE_4 <- 0
ProtAtlas$HSChloroplast_CT_1 <- 0
ProtAtlas$HSChloroplast_CT_2 <- 0
ProtAtlas$HSChloroplast_CT_3 <- 0
ProtAtlas$HSChloroplast_CT_4 <- 0
ProtAtlas$HSChloroplast_TE1_1 <- 0
ProtAtlas$HSChloroplast_TE1_2 <- 0
ProtAtlas$HSChloroplast_TE1_3 <- 0
ProtAtlas$HSChloroplast_TE1_4 <- 0
ProtAtlas$HSChloroplast_TT1_1 <- 0
ProtAtlas$HSChloroplast_TT1_2 <- 0
ProtAtlas$HSChloroplast_TT1_3 <- 0
ProtAtlas$HSChloroplast_TT1_4 <- 0
ProtAtlas$HSChloroplast_TE3_1 <- 0
ProtAtlas$HSChloroplast_TE3_2 <- 0
ProtAtlas$HSChloroplast_TE3_3 <- 0
ProtAtlas$HSChloroplast_TE3_4 <- 0
ProtAtlas$HSChloroplast_TT3_1 <- 0
ProtAtlas$HSChloroplast_TT3_2 <- 0
ProtAtlas$HSChloroplast_TT3_3 <- 0
ProtAtlas$HSChloroplast_TT3_4 <- 0
ProtAtlas$HSChloroplast_TE5_1 <- 0
ProtAtlas$HSChloroplast_TE5_2 <- 0
ProtAtlas$HSChloroplast_TE5_3 <- 0
ProtAtlas$HSChloroplast_TE5_4 <- 0
ProtAtlas$HSChloroplast_TT5_1 <- 0
ProtAtlas$HSChloroplast_TT5_2 <- 0
ProtAtlas$HSChloroplast_TT5_3 <- 0
ProtAtlas$HSChloroplast_TT5_4 <- 0

ProtAtlas$UVChloroplast_CE_1 <- 0
ProtAtlas$UVChloroplast_CE_2 <- 0
ProtAtlas$UVChloroplast_CE_3 <- 0
ProtAtlas$UVChloroplast_CE_4 <- 0
ProtAtlas$UVChloroplast_CT_1 <- 0
ProtAtlas$UVChloroplast_CT_2 <- 0
ProtAtlas$UVChloroplast_CT_3 <- 0
ProtAtlas$UVChloroplast_CT_4 <- 0
ProtAtlas$`UVChloroplast_TE1/2_1` <- 0
ProtAtlas$`UVChloroplast_TE1/2_2` <- 0
ProtAtlas$`UVChloroplast_TE1/2_3` <- 0
ProtAtlas$`UVChloroplast_TE1/2_4` <- 0
ProtAtlas$`UVChloroplast_TT1/2_1` <- 0
ProtAtlas$`UVChloroplast_TT1/2_2` <- 0
ProtAtlas$`UVChloroplast_TT1/2_3` <- 0
ProtAtlas$`UVChloroplast_TT1/2_4` <- 0
ProtAtlas$UVChloroplast_TE1_1 <- 0
ProtAtlas$UVChloroplast_TE1_2 <- 0
ProtAtlas$UVChloroplast_TE1_3 <- 0
ProtAtlas$UVChloroplast_TE1_4 <- 0
ProtAtlas$UVChloroplast_TT1_1 <- 0
ProtAtlas$UVChloroplast_TT1_2 <- 0
ProtAtlas$UVChloroplast_TT1_3 <- 0
ProtAtlas$UVChloroplast_TT1_4 <- 0
ProtAtlas$UVChloroplast_TE2_1 <- 0
ProtAtlas$UVChloroplast_TE2_2 <- 0
ProtAtlas$UVChloroplast_TE2_3 <- 0
ProtAtlas$UVChloroplast_TE2_4 <- 0
ProtAtlas$UVChloroplast_TT2_1 <- 0
ProtAtlas$UVChloroplast_TT2_2 <- 0
ProtAtlas$UVChloroplast_TT2_3 <- 0
ProtAtlas$UVChloroplast_TT2_4 <- 0

# fill all columns and export to data.frame
which((Tissue.Preprocessed$ProtID %in% ProtAtlas$ProtID) == FALSE)
ProtAtlas[ProtAtlas$ProtID %in% Tissue.Preprocessed$ProtID,c(3:14)] <- Tissue.Preprocessed[,c(2:13)]

which((FUTotal.Preprocessed.Imputed.F$ProtID %in% ProtAtlas$ProtID) == FALSE)
ProtAtlas[ProtAtlas$ProtID %in% FUTotal.Preprocessed.Imputed.F$ProtID,15:24] <- FUTotal.Preprocessed.Imputed.F[,c(2:6,8:12)]

which((HSTotal.Preprocessed.Imputed.F$ProtID %in% ProtAtlas$ProtID) == FALSE)
ProtAtlas[ProtAtlas$ProtID %in% HSTotal.Preprocessed.Imputed.F$ProtID,25:36] <- HSTotal.Preprocessed.Imputed.F[,c(2:13)]

which((UVTotal.Preprocessed$ProtID %in% ProtAtlas$ProtID) == FALSE)
ProtAtlas[ProtAtlas$ProtID %in% UVTotal.Preprocessed$ProtID,37:51] <- UVTotal.Preprocessed[,c(2:16)]

which((HSNucleus.Preprocessed$ProtID %in% ProtAtlas$ProtID) == FALSE)
ProtAtlas[ProtAtlas$ProtID %in% HSNucleus.Preprocessed$ProtID,52:79] <- HSNucleus.Preprocessed[,c(2:29)]

which((UVNucleus.Preprocessed.Imputed.F$ProtID %in% ProtAtlas$ProtID) == FALSE)
ProtAtlas[ProtAtlas$ProtID %in% UVNucleus.Preprocessed.Imputed.F$ProtID,80:91] <- UVNucleus.Preprocessed.Imputed.F[,c(2:13)]

which((HSChloroplast.Preprocessed$ProtID %in% ProtAtlas$ProtID) == FALSE)
ProtAtlas[ProtAtlas$ProtID %in% HSChloroplast.Preprocessed$ProtID,92:123] <- HSChloroplast.Preprocessed[,c(2:33)]

which((UVChloroplast.Preprocessed$ProtID %in% ProtAtlas$ProtID) == FALSE)
ProtAtlas[ProtAtlas$ProtID %in% UVChloroplast.Preprocessed$ProtID,124:155] <- UVChloroplast.Preprocessed[,c(2:33)]

ProtAtlas.jic <- ProtAtlas
write.table(ProtAtlas, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/Protein_module/ProtAtlas.Analyses.txt", sep = "\t", quote = F, col.names = T)
write.table(ProtAtlas, file = "F:/ATLAS_Pra/Protein_module/0.Process/ProtAtlas.Analyses.txt", sep = "\t", quote = F, col.names = T)



#### 2. Global characterization ####

PROTAS <- read.delim(file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/Protein_module/ProtAtlas.Analyses.txt", header = T, sep = "\t")
PROTAS <- t(PROTAS)
colnames(PROTAS) <- PROTAS[1,]
PROTAS <- PROTAS[-c(1:2),]
Samples <- paste0("F",rownames(PROTAS))

Tissues <- c(rep("Needle", 6), rep("Root", 3), rep("Buds", 13), rep("Needle", 131))
Stress <- c(rep("Control", 17), rep("FU", 5), rep("Control", 4), rep("Heat", 8), rep("Control", 3), rep("UV", 9), rep("Recovery", 3), 
            rep("Control", 4), rep("Heat", 16), rep("Recovery", 8), 
            rep("Control", 4),  rep("UV", 8), 
            rep("Control", 8), rep("Heat", 24),
            rep("Control", 8), rep("UV", 24))
Population <- c(rep("None", 17), rep("None", 5), rep("None", 4), rep("None", 8), rep("None", 3), rep("None", 9), rep("None", 3), 
            rep("None", 4), rep("None", 16), rep("None", 8), 
            rep("None", 4),  rep("None", 8), 
            rep("E", 4), rep("T", 4), rep("E", 4), rep("T", 4), rep("E", 4), rep("T", 4), rep("E", 4), rep("T", 4),
            rep("E", 4), rep("T", 4), rep("E", 4), rep("T", 4), rep("E", 4), rep("T", 4), rep("E", 4), rep("T", 4)
            )
Technique <- c(rep("Total", 17), rep("Total", 5), rep("Total", 4), rep("Total", 8), rep("Total", 3), rep("Total", 9), rep("Total", 3), 
               rep("Nucleus", 4), rep("Nucleus", 16), rep("Nucleus", 8), 
               rep("Nucleus", 4),  rep("Nucleus", 8), 
               rep("Chloroplast", 8), rep("Chloroplast", 24),
               rep("Chloroplast", 8), rep("Chloroplast", 24))
Intensity <- c(rep("Control", 17), rep("T1", 5), rep("Control", 4), rep("T1", 4), rep("T3", 4), rep("Control", 3), rep("T1", 3), rep("T3", 3), rep("T5", 3), rep("Recovery", 3), 
               rep("Control", 4), rep("T1", 4), rep("T3", 4), rep("T5", 4), rep("T10", 4), rep("Recovery", 8), 
               rep("Control", 4),  rep("T1", 4), rep("T3", 4), 
               rep("Control", 8), rep("T1", 8), rep("T3", 8), rep("T5", 8),
               rep("Control", 8), rep("T1", 8), rep("T3", 8), rep("T5", 8))

Stress.Intensity <- paste(Stress, Intensity, sep = "_")
Population.Stress <- paste(Population, Stress, sep = "_")
Technique.Stress <- paste(Technique, Stress, sep = "_")
Population.Intensity <- paste(Population, Intensity, sep = "_")
Technique.Intensity <- paste(Technique, Intensity, sep = "_")
Stress.Intensity <- gsub("Control_Control", "Control", Stress.Intensity)
Stress.Intensity <- gsub("Recovery_Recovery", "Recovery", Stress.Intensity)
labelitass <- list(Stress = Stress, Intensity = Intensity, Population = Population, Technique = Technique, Tissues = Tissues, Stress.Intensity = Stress.Intensity, Population.Stress = Population.Stress, Population.Intensity = Population.Intensity,
                   Technique.Stress = Technique.Stress, Technique.Intensity = Technique.Intensity)
ncolors <- lapply(labelitass, function(x) length(levels(as.factor(x))))
fcolors <- c(RColorBrewer::brewer.pal(12, "Paired") ,"black",  "peachpuff3", "lightgoldenrod")
color.labelitass <- list(Stress = c("lightsteelblue3", "darkgreen", "salmon3", "lightsteelblue1", "purple"),
                         Intensity = c("lightsteelblue3", "lightsteelblue1", "salmon1", "salmon4", "salmon2", "salmon3"),
                         Population = c("purple", "black", "gold"),
                         Technique = c("darkgreen", "purple", "gold"),
                         Tissues = c("darkgreen", "purple", "gold"),
                         Stress.Intensity = c("lightsteelblue3", "darkgreen", "salmon1", "salmon4", "salmon2", "salmon3", "lightsteelblue1", "mediumpurple1", "mediumpurple2", "mediumpurple3"),
                         Population.Stress = fcolors[1:11],
                         Population.Intensity =  fcolors[1:14],
                         Technique.Stress = fcolors[1:12],
                         Technique.Intensity = fcolors[1:15])


#### Descriptives ####

PROTAS <- apply(PROTAS, MARGIN = 2, FUN = as.numeric)

PROTAS.zscore <- apply(PROTAS, MARGIN = 2, FUN = function(x){
  x <- (x - mean(x))/sd(x)
  x
})

PROTAS.log10 <- apply(PROTAS, MARGIN = 2, FUN = function(x){
  x <- log10(x + 1.1)
  x
})

PROTAS.log10.zscore <- apply(PROTAS, MARGIN = 2, FUN = function(x){
  x <- log10(x + 1)
  x <- (x - mean(x))/sd(x)
  x
})

PROTAS.quantiles <- apply(PROTAS, MARGIN = 2, FUN = function(x){
  x <- normalizeQuantiles(x, ties = TRUE)
  x
})

PROTAS.raw <- PROTAS

rawinputs <- list(raw = PROTAS.raw, zscore = PROTAS.zscore, log10 = PROTAS.log10, log10zscore = PROTAS.log10.zscore, quantile = PROTAS.quantiles)

#### PCA ####

inputs <- lapply(rawinputs, function(x) prcomp(x, center = F, scale. = F, rank. = 10))
fviz_eig(inputs$raw)
fviz_eig(inputs$zscore)
fviz_eig(inputs$log10)
fviz_eig(inputs$log10zscore)
fviz_eig(inputs$quantile)

for(j in 1:length(inputs)){
  for(i in 1:length(labelitass)){
    data <- as.data.frame(inputs[[j]]$x)
    meas_vars <- colnames(data)
    control.Table <- data.frame(expand.grid(meas_vars, meas_vars, stringsAsFactors = FALSE))
    colnames(control.Table) <- c("x", "y")
    control.Table <- cbind(data.frame(pair_key = paste(control.Table[[1]],
                                                       control.Table[[2]], sep = "-"),
                                      stringsAsFactors = FALSE),
                           control.Table)
    data <- cbind(label = labelitass[[i]], data)
    data_aug <- rowrecs_to_blocks(data, control.Table, columnsToCopy = "label")
    splt <- strsplit(data_aug$pair_key, split = "-", fixed = T)
    data_aug$xv <- vapply(splt, function(si) si[[1]], character(1))
    data_aug$yv <- vapply(splt, function(si) si[[2]], character(1))
    data_aug$xv <- factor(as.character(data_aug$xv), meas_vars)
    data_aug$yv <- factor(as.character(data_aug$yv), meas_vars)
    if(length(levels(as.factor(labelitass[[i]]))) >= 6){
      print(ggplot(data_aug, aes(x = x, y = y)) +
              geom_point(aes(color=label, shape ="circle")) + scale_color_manual(values = color.labelitass[[i]]) + 
              facet_grid(yv~xv, labeller = label_both, scales = "free")  +
              ggtitle("PCA Biplot Matrix") + 
              ylab(NULL) + xlab(NULL) + theme_light() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()))
    }else if(length(levels(as.factor(labelitass[[i]]))) < 6){
      print(ggplot(data_aug, aes(x = x, y = y)) +
              geom_point(aes(color=label, shape =label)) + scale_color_manual(values = color.labelitass[[i]]) +
              facet_grid(yv~xv, labeller = label_both, scales = "free") +
              ggtitle("PCA Biplot Matrix") + 
              ylab(NULL) + xlab(NULL) + theme_light() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()))
    }
    
    a <- readline(prompt = "Enter PCA number: ")
    b <- readline(prompt = "Enter PCA number: ")
    c <- readline(prompt = "Enter PCA number: ")
    LNs <- as.numeric(c(a,b,c))
    par("ask" = TRUE)
    colors <- color.labelitass[[i]][as.numeric(as.factor(data$label))]
    #print(plot_ly(as.data.frame(data), x = data[,LNs[1]+1], y = data[,LNs[2]+1], z = data[,LNs[3]+1],
    #              color = data$label, colors = colors) %>%
    #        add_markers() %>%
    #        layout(scene = list(xaxis = list(title = paste0("PCA", a)),
    #                            yaxis = list(title = paste0("PCA", b)),
    #                            zaxis = list(title = paste0("PCA",c)))))
    
    print(scatterplot3d(x = data[,LNs[1]+1], y = data[,LNs[2]+1], z = data[,LNs[3]+1], color = colors, pch = 16, xlab = paste0("PCA", LNs[1]), ylab = paste0("PCA", LNs[2]), zlab = paste0("PCA", LNs[3]), type = "h"))
    legend("right", legend = data$label[!duplicated(data$label)],
           col =  colors[!duplicated(colors)], pch = 16)
    
  } # i conditions
} # j inpu

#### UMAP ####

set.seed(404)
inputs <- lapply(rawinputs, function(x) umap(x, n_components = 5, n_threads = 7))
inputs.pca <- lapply(rawinputs, function(x) umap(x, n_components = 5, n_threads = 7, pca = 50))
inputs <- lapply(inputs, function(x){ 
  colnames(x) <- paste0("UMAP", 1:ncol(x))
  x
})

inputs.pca <- lapply(inputs.pca, function(x){ 
  colnames(x) <- paste0("UMAP", 1:ncol(x))
  x
})

INPUTS.F <- list(WI = inputs.pca, WO = inputs)

write.table(INPUTS.F$WO$log10, file = "F:/ATLAS_Pra/Protein_module/0.Process/UMAPtable_WOlog10.txt", sep = "\t", col.names = T, row.names = T)

for(j in 1:length(INPUTS.F)){
  for(z in 1:length(INPUTS.F[[j]])){
    for(i in 1:length(labelitass)){
      data <- as.data.frame(INPUTS.F[[j]][[z]])
      meas_vars <- colnames(data)
      control.Table <- data.frame(expand.grid(meas_vars, meas_vars, stringsAsFactors = FALSE))
      colnames(control.Table) <- c("x", "y")
      control.Table <- cbind(data.frame(pair_key = paste(control.Table[[1]],
                                                         control.Table[[2]], sep = "-"),
                                        stringsAsFactors = FALSE),
                             control.Table)
      data <- cbind(label = labelitass[[i]], data)
      data_aug <- rowrecs_to_blocks(data, control.Table, columnsToCopy = "label")
      splt <- strsplit(data_aug$pair_key, split = "-", fixed = T)
      data_aug$xv <- vapply(splt, function(si) si[[1]], character(1))
      data_aug$yv <- vapply(splt, function(si) si[[2]], character(1))
      data_aug$xv <- factor(as.character(data_aug$xv), meas_vars)
      data_aug$yv <- factor(as.character(data_aug$yv), meas_vars)
      if(length(levels(as.factor(labelitass[[i]]))) >= 6){
        print(ggplot(data_aug, aes(x = x, y = y)) +
                geom_point(aes(color=label, shape ="circle")) + scale_color_manual(values = color.labelitass[[i]]) + 
                facet_grid(yv~xv, labeller = label_both, scales = "free") +
                ggtitle("UMAP Biplot Matrix") + 
                ylab(NULL) + xlab(NULL) + theme_light() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()))
      }else if(length(levels(as.factor(labelitass[[i]]))) < 6){
        print(ggplot(data_aug, aes(x = x, y = y)) +
                geom_point(aes(color=label, shape =label)) + scale_color_manual(values = color.labelitass[[i]]) + 
                facet_grid(yv~xv, labeller = label_both, scales = "free") +
                ggtitle("UMAP Biplot Matrix") + 
                ylab(NULL) + xlab(NULL) + theme_light() + theme(axis.text.x = element_blank(), axis.text.y = element_blank()))
      }
      
      a <- readline(prompt = "Enter UMAP number: ")
      b <- readline(prompt = "Enter UMAP number: ")
      c <- readline(prompt = "Enter UMAP number: ")
      LNs <- as.numeric(c(a,b,c))
      par("ask" = TRUE)
      colors <- color.labelitass[[i]][as.numeric(as.factor(data$label))]
      #print(plot_ly(as.data.frame(data), x = data[,LNs[1]+1], y = data[,LNs[2]+1], z = data[,LNs[3]+1],
      #              color = data$label, colors = colors) %>%
      #        add_markers() %>%
      #        layout(scene = list(xaxis = list(title = paste0("UMAP", a)),
      #                            yaxis = list(title = paste0("UMAP", b)),
      #                            zaxis = list(title = paste0("UMAP",c)))))
      
      print(scatterplot3d(x = data[,LNs[1]+1], y = data[,LNs[2]+1], z = data[,LNs[3]+1], color = colors, pch = 16, xlab = paste0("UMAP", LNs[1]), ylab = paste0("UMAP", LNs[2]), zlab = paste0("UMAP", LNs[3]), type = "h"))
      legend("right", legend = data$label[!duplicated(data$label)],
             col =  colors[!duplicated(colors)], pch = 16)
      
    } # labelitass
  } # intra-WI-WO
} # INPUTS.F
# trying to global biplot picture

#### Differential proteins ####
## Too many contrast; we only select the most informative ones 
## Needle vs Roots vs Floral buds (New tissues set) | 1:12 Tissues
## FU vs Control (total) | 13:22 FU
## UV vs Control (Total) | 35:49 UV Total
## UV vs Control (Nucleus) | 78:89 UV Nucleus
## UV vs Control (Chloroplast) | 122:153 UVChloroplast
## Heat vs Control (Total) | 23:34 Heat total
## Heat vs Control (Nucleus) | 50:77 Heat Nucleus
## Heat vs Control (Chloroplast) | 90:121 Heat Chloroplast
## Heat vs UV (Total) | c(27:34, 38:46) Heat vs UV Total
## Heat vs UV (Nucleus) | c(55:69,82:88) Heat vs UV Nucleus
## Heat vs UV (Chloroplast) | c(98:121, 130:153) Heat vs UV Chloroplast
## UV vs Recovery (UV) | PREV DONE
## Heat vs Recovery (Heat) | PREV DONE 


inputs <- list(log10 = PROTAS.log10[c(98:121, 130:153),], log10.zscore = PROTAS.log10.zscore[c(98:121, 130:153),])
Differential <- list()
sva <- TRUE
conditions <- list(Tissues = labelitass$Tissues, Stress = labelitass$Stress, Technique = labelitass$Technique)
conditions <- list(Heat_UV_Chloroplast = labelitass$Stress[c(98:121, 130:153)])
output.path <- "F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/"


for(l in 2:length(inputs)){
  for(v in 1:length(conditions)){
    con <- file(paste0(output.path, names(conditions[v]), "_sva_", names(inputs[l]),".log"))
    sink(con, append=TRUE)
    sink(con, append=TRUE, type="message")
    input <- t(inputs[[l]])
    colnames(input) <- Samples[c(98:121, 130:153)]
    groups <- as.numeric(as.factor(conditions[[v]]))
    comb.contrast <- combn(groups[!duplicated(groups)], 2)
    cat(paste0("\n" ,names(inputs[l]), sep = ": ", length(levels(as.factor(groups))), " treatments and ", ncol(comb.contrast), " two-by-two contrasts \n"))
    groups.list <- list()
    for(z in groups[!duplicated(groups)]){
      groups.list[[z]] <- which(groups == z)
      }
    for(i in 1:ncol(comb.contrast)){
      cat("\n Performing differential expression ... \n")
      # Build null-model and model
      samples <- c(groups.list[[comb.contrast[1,i]]], groups.list[[comb.contrast[2,i]]])
      cat(paste0("\n Dataset: ", names(inputs[l]), " norm, ", names(conditions[v]), " grouping, ", paste0(levels(as.factor(conditions[[v]][samples])), collapse = "-"), " contrasts."))
      PhenoData <- data.frame(samples = colnames(input[,samples]), Treatment = conditions[[v]][samples])
      mod = model.matrix(~0+as.factor(Treatment), data=PhenoData)
      colnames(mod) <- c(levels(as.factor(PhenoData$Treatment)))
      mod0 = model.matrix(~1,data=PhenoData)
      colnames(mod0) <- "samples"
      # Side_note01: In the future remove infinte and zero values if svas are giving problems
      if(!sva){
        # Fit model
        fit <-lmFit(input[,samples],mod)
        # Select contrast
        contrastss <- c(paste0(colnames(mod)[2], "-", colnames(mod)[1]))
        contrast.matrix <- makeContrasts(contrastss, levels=mod)
        fit2 <- contrasts.fit(fit,contrast.matrix)
        fit2 <- eBayes(fit2)
        # Filter and export
        DE_all <- limma::topTable(fit2, number = Inf, coef = 1, sort.by = "logFC")
        DE_filtered <- limma::topTable(fit2, number = Inf, coef = 1, p.value = 0.05, sort.by = "logFC")
        cat(paste0("\n \t ", contrastss, sep = " | ", nrow(input), " total proteins , ", nrow(DE_filtered), " differential proteins (0.05) \n"))
        file.name <- paste0(names(inputs[l]), "_",contrastss)
        write.table(DE_filtered, file = paste0(output.path, names(conditions[v]), "_", file.name, "_filtered_", names(conditions[v]),".txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = F)
        write.table(DE_all, file = paste0(output.path, names(conditions[v]), "_", file.name, "_all_", names(conditions[v]),".txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = F)
        Differential[[paste0(file.name, "_filtered_", names(conditions[v]), "_nosva")]] <- DE_filtered
        Differential[[paste0(file.name,"_all_", names(conditions[v]), "_nosva")]] <- DE_all
      }else if(sva){
        # sva
        n.sv = num.sv(as.matrix(input[,samples]),mod,method="leek")
        cat(paste0("\n \t sva: ", n.sv, " unknown batch effects founded \n"))
        if(n.sv == 0){
          # Fit model
          fit <-lmFit(input[,samples],mod)
          # Select contrast
          contrastss <- c(paste0(colnames(mod)[2], "-", colnames(mod)[1]))
          contrast.matrix <- makeContrasts(contrastss, levels=mod)
          fit2 <- contrasts.fit(fit,contrast.matrix)
          fit2 <- eBayes(fit2)
          # Filter and export
          DE_all <- limma::topTable(fit2, number = Inf, coef = 1, sort.by = "logFC")
          DE_filtered <- limma::topTable(fit2, number = Inf, coef = 1, p.value = 0.05, sort.by = "logFC")
          cat(paste0("\n \t ", contrastss, sep = " | ", nrow(input), " total proteins , ", nrow(DE_filtered), " differential proteins (0.05) \n"))
          file.name <- paste0(names(inputs[l]), "_",contrastss)
          write.table(DE_filtered, file = paste0(output.path, names(conditions[v]), "_", file.name, "_filtered_", names(conditions[v]),".txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = F)
          write.table(DE_all, file = paste0(output.path, names(conditions[v]), "_", file.name, "_all_", names(conditions[v]),".txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = F)
          Differential[[paste0(file.name, "_filtered_", names(conditions[v]))]] <- DE_filtered
          Differential[[paste0(file.name,"_all_", names(conditions[v]))]] <- DE_all
        }else if(n.sv > 0){
          # Model correction by sva
          svobj = sva(as.matrix(input[,samples]),mod,mod0,n.sv=n.sv)
          colnames(svobj$sv) <- paste(rep("col",ncol(svobj$sv)),c(1:ncol(svobj$sv)),sep="")
          modSv = cbind(mod,svobj$sv)
          # Fit model
          fit <-lmFit(input[,samples],modSv)
          # Select contrast
          contrastss <- c(paste0(colnames(modSv)[2], "-", colnames(modSv)[1]))
          contrast.matrix <- makeContrasts(contrastss, levels=modSv)
          fit2 <- contrasts.fit(fit,contrast.matrix)
          fit2 <- eBayes(fit2)
          # Filter and export
          DE_all <- limma::topTable(fit2, number = Inf, coef = 1, sort.by = "logFC")
          DE_filtered <- limma::topTable(fit2, number = Inf, coef = 1, p.value = 0.05, sort.by = "logFC")
          cat(paste0("\n \t ", contrastss, sep = " | ", nrow(input), " total proteins , ", nrow(DE_filtered), " differential proteins (0.05) \n"))
          file.name <- paste0(names(inputs[l]), "_",contrastss)
          write.table(DE_filtered, file = paste0(output.path, names(conditions[v]), "_", file.name, "_filtered_", names(conditions[v]),".txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = F)
          write.table(DE_all, file = paste0(output.path, names(conditions[v]), "_", file.name, "_all_", names(conditions[v]),".txt"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = F)
          Differential[[paste0(file.name, "_filtered_", names(conditions[v]))]] <- DE_filtered
          Differential[[paste0(file.name, "_all_", names(conditions[v]))]] <- DE_all
        }
      }
    }
    
    sink() 
    sink(type="message")
  }
}


#### Enrich heatmaps ####

## Create annotation table per ProteinID using transcriptID database

PROTAS.annotation <- data.frame(ProtID = colnames(PROTAS), Mercator.Bin = NA, GO = NA, Symbol = NA)
Annotation <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/4.FinalAnnotation/0.final/Pra-GE-ATLAS_dammitMercatorIP5EggnogProtID.txt")

## Symbol tunning for later use in volcano plots

Annotation$Symbol <- NA
table(Annotation$EGGNOGMAPPER_DESC)[1:10]
Annotation$EGGNOGMAPPER_DESC <- gsub("-,", "", fixed = T, x = Annotation$EGGNOGMAPPER_DESC)
Annotation$EGGNOGMAPPER_DESC <- gsub(",-", "", fixed = T, x = Annotation$EGGNOGMAPPER_DESC)
Annotation$EGGNOGMAPPER_DESC <- gsub("-", "", fixed = T, x = Annotation$EGGNOGMAPPER_DESC)
Annotation$EGGNOGMAPPER_DESC[Annotation$EGGNOGMAPPER_DESC == ""] <- NA
table(Annotation$EGGNOGMAPPER_DESC)[1:10]
table(Annotation$IP5_Description)[1:75]
Annotation$IP5_Description <- gsub("-,", "", fixed = T, x = Annotation$IP5_Description)
Annotation$IP5_Description <- gsub(",-", "", fixed = T, x = Annotation$IP5_Description)
Annotation$IP5_Description <- gsub("-", "", fixed = T, x = Annotation$IP5_Description)
Annotation$IP5_Description[Annotation$IP5_Description == ""] <- NA
table(Annotation$IP5_Description)[1:10]
Symbols <- gsub("-,", "", x = Annotation$EGGNOGMAPPER_SYMBOL[is.na(Annotation$Symbol)], fixed = T)
Symbols <- gsub(",-", "", x = Symbols, fixed = T)
Symbols <- gsub("-", "", x = Symbols, fixed = T)
Symbols[Symbols == ""] <- NA
Annotation$Symbol <- Symbols
Annotation$Symbol[is.na(Annotation$Symbol)] <- Annotation$EGGNOGMAPPER_DESC[is.na(Annotation$Symbol)]
Annotation$Symbol[is.na(Annotation$Symbol)] <- Annotation$IP5_Description[is.na(Annotation$Symbol)]
Annotation$Symbol[is.na(Annotation$Symbol)] <- Annotation$IP5_SignatureAccesion[is.na(Annotation$Symbol)]
Annotation$Symbol[is.na(Annotation$Symbol)] <- Annotation$TranscriptID[is.na(Annotation$Symbol)]

## Quality check again

ALL <- unlist(strsplit(x = Annotation$ProteinID, split = ",", fixed = T))
ALL <- ALL[!is.na(ALL)]
ALL <- unique(ALL)
length(which((PROTAS.annotation$ProtID %in% ALL) == FALSE)) # all hits are contained in database annotation

## ProteinID annotation table generation

check <- separate_longer_delim(Annotation, ProteinID, delim = ",")

for(i in 1:nrow(PROTAS.annotation)){
  hit <- which(check$ProteinID == PROTAS.annotation$ProtID[i])
  if(length(hit) == 0){stop("Something went worng")}
  if(length(hit) > 1){
    Mercator.hit <- unique(check$MercatorBin[hit])
    if(length(Mercator.hit) == 1){
      PROTAS.annotation$Mercator.Bin[i] <- Mercator.hit
    }else if(length(Mercator.hit) > 1){
      PROTAS.annotation$Mercator.Bin[i] <- paste(Mercator.hit, collapse = ",")
    }
    GO.hit <- unique(check$EGGNOGMAPPER_GO[hit])
    if(length(GO.hit) == 1){
      PROTAS.annotation$GO[i] <- GO.hit
    }else if(length(GO.hit) > 1){
      PROTAS.annotation$GO[i] <- paste(GO.hit, collapse = ",")
    }
    Symbol.hit <- unique(check$Symbol[hit])
    if(length(Symbol.hit) == 1){
      PROTAS.annotation$Symbol[i] <- Symbol.hit
    }else if(length(Symbol.hit) > 1){
      PROTAS.annotation$Symbol[i] <- paste(Symbol.hit, collapse = ",")
    }
  }else if(length(hit) == 1){
    PROTAS.annotation$Mercator.Bin[i] <- check$MercatorBin[hit]
    PROTAS.annotation$GO[i] <- check$EGGNOGMAPPER_GO[hit]
    PROTAS.annotation$Symbol[i] <- check$Symbol[hit]
  }
  
}

write.table(x = PROTAS.annotation, file = "F:/ATLAS_Pra/Protein_module/0.Process/Tables/PROTAS.Annotation.txt", quote = F, sep = "\t")
rownames(PROTAS.annotation) <- PROTAS.annotation$ProtID

## Define sets

Tissues <- list(Needle_Buds = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Tissues/Tissues_log10.zscore_Needle-Buds_all_Tissues.txt", header = T, sep = "\t", row.names = 1),
                Root_Needle = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Tissues/Tissues_log10.zscore_Root-Needle_all_Tissues.txt", header = T, sep = "\t", row.names = 1),
                Root_Buds = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Tissues/Tissues_log10.zscore_Root-Buds_all_Tissues.txt", header = T, sep = "\t", row.names = 1))
Tissues.res <- list()

Stress <- list(Total_FU = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Fusarium/Fusarium_log10.zscore_FU-Control_all_Fusarium.txt", header = T, sep = "\t", row.names = 1),
               Total_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/UVTotal/UVTotal_log10.zscore_UV-Control_all_UVTotal.txt", header = T, sep = "\t", row.names = 1),
               Total_Heat = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/HeatTotal/HeatTotal_log10.zscore_Heat-Control_all_HeatTotal.txt", header = T, sep = "\t", row.names = 1),
               Nucleus_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/UVNucleus/UVNucleus_log10.zscore_UV-Control_all_UVNucleus.txt", header = T, sep = "\t", row.names = 1),
               Nucleus_Heat = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/HeatNucleus/HeatNucleus_log10.zscore_Heat-Control_all_HeatNucleus.txt", header = T, sep = "\t", row.names = 1),
               Chloroplast_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/UVChloroplast/UVChloroplast_log10.zscore_UV-Control_all_UVChloroplast.txt", header = T, sep = "\t", row.names = 1),
               Chloroplast_Heat = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/HeatChloroplast/HeatChloroplast_log10.zscore_Heat-Control_all_HeatChloroplast.txt", header = T, sep = "\t", row.names = 1),
               Recovery_Heat = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/HeatNucleus/HeatNucleus_log10.zscore_Recovery-Heat_all_HeatNucleus.txt", header = T, sep = "\t", row.names = 1),
               Recovery_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/UVTotal/UVTotal_log10.zscore_UV-Recovery_all_UVTotal.txt", header = T, sep = "\t", row.names = 1),
               Total_Heat_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Heat_UV_Total/Heat_UV_Total_log10.zscore_UV-Heat_all_Heat_UV_Total.txt", header = T, sep = "\t", row.names = 1),
               Nucleus_Heat_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Heat_UV_Nucleus/Heat_UV_Nucleus_log10.zscore_UV-Heat_all_Heat_UV_Nucleus.txt", header = T, sep = "\t", row.names = 1),
               Chloroplast_Heat_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Heat_UV_Chloroplast/Heat_UV_Chloroplast_log10.zscore_UV-Heat_all_Heat_UV_Chloroplast.txt", header = T, sep = "\t", row.names = 1))
Stress.res <- list()

## Compute enrichments

outputs <- "F:/ATLAS_Pra/Protein_module/0.Process/Figures/Descriptives/Proteins/Enrich/"
output.path <- "F:/ATLAS_Pra/Protein_module/0.Process/Figures/Descriptives/Proteins/Enrich/"
PROTAS.annotation.bin <- separate_longer_delim(PROTAS.annotation, Mercator.Bin, delim = ",")
Bin_db <- list()
bin.names <- levels(as.factor(PROTAS.annotation.bin$Mercator.Bin))
for(i in bin.names){
  Bin_db[[i]] <- PROTAS.annotation.bin$ProtID[PROTAS.annotation.bin$Mercator.Bin == i]
}
names(Bin_db) <- c("PS", "Redox hom", "Phytohormone act", "Chromatin org", "Cell division", "DNA dmg response", "RNA biosynthesis", "RNA processing", 
                   "Protein biosynthesis", "Protein modification", "Protein hom",
                   "Cellular respiration", "Cytoskeleton org", "Cell wall org", "Vessicle trafficking", "Protein translocation", "Solute transport",
                   "Nutrient uptake", "External stimuli resp", "Multi-process reg", "Plant reproduction", 
                   "Carbohydrate met", "Clade-specific met", "Unknown", 
                   "Aminoacid met", "Lipid met", "Enzyme classification", "Nucleotide met", "Coenzyme met", "Polyamine met", "Secondary met")

for(i in names(Tissues)){
  data <- Tissues[[i]]
  data$ID <- rownames(data)
  data <- data[,c("ID", "logFC")]
  data <- deframe(data)
  fGseaRES <- fgsea(pathways = Bin_db, stats = data, minSize = 5, maxSize = 1000)
  fgseaResTidy <- fGseaRES %>%
  as_tibble() %>%
  arrange(desc(NES))
  write.table(fgseaResTidy[,-ncol(fgseaResTidy)], file = paste0(output.path, i, "_fgsea.txt"), col.names = T, sep = "\t", quote = F)
  Tissues.res[[i]] <- fgseaResTidy 
}

# Store Tissues.res for heatmap plotting together with Stress.res fgsea and WGCNA modules ORA

for(i in names(Stress)){
  data <- Stress[[i]]
  data$ID <- rownames(data)
  data <- data[,c("ID", "logFC")]
  data <- deframe(data)
  fGseaRES <- fgsea(pathways = Bin_db, stats = data, minSize = 5, maxSize = 1000)
  fgseaResTidy <- fGseaRES %>%
    as_tibble() %>%
    arrange(desc(NES))
  write.table(fgseaResTidy[,-ncol(fgseaResTidy)], file = paste0(output.path, i, "_fgsea.txt"), col.names = T, sep = "\t", quote = F)
  Stress.res[[i]] <- fgseaResTidy 
}

# Store Stress.res for heatmap plotting together with Tissues.res fgsea and WGCNA modules ORA

#### Plot Heatmaps ####

col_significance <- colorRamp2(c(0,1), c("white", "dodgerblue4"))
col_quantity <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

## Tissues

Tissues.res[["Needle_Buds"]] <- read.delim("clipboard")
Tissues.res[["Root_Needle"]] <- read.delim("clipboard")
Tissues.res[["Root_Buds"]] <- read.delim("clipboard")

Heatmap.tissue.Q <- data.frame(Factor = Tissues.res$Needle_Buds$pathway, Needle_Buds = NA, Root_Needle = NA, Root_Buds = NA)
Heatmap.tissue.S <- data.frame(Factor = Tissues.res$Needle_Buds$pathway, Needle_Buds = NA, Root_Needle = NA, Root_Buds = NA)

for(i in 1:nrow(Heatmap.tissue.Q)){
  for(j in colnames(Heatmap.tissue.Q)[-1]){
    Heatmap.tissue.Q[i,j] <- Tissues.res[[j]]$NES[which(Tissues.res[[j]]$pathway == Heatmap.tissue.Q$Factor[i])]
  }
}
for(i in 1:nrow(Heatmap.tissue.S)){
  for(j in colnames(Heatmap.tissue.S)[-1]){
    Heatmap.tissue.S[i,j] <- Tissues.res[[j]]$padj[which(Tissues.res[[j]]$pathway == Heatmap.tissue.Q$Factor[i])]
    if(Heatmap.tissue.S[i,j] <= 0.1){
      Heatmap.tissue.S[i,j] <- 0
    }else if(Heatmap.tissue.S[i,j] > 0.1){
      Heatmap.tissue.S[i,j] <- 1
    }
  }
}

rownames(Heatmap.tissue.Q) <- Heatmap.tissue.Q$Factor
rownames(Heatmap.tissue.S) <- Heatmap.tissue.S$Factor

ht.Q.tissue <- Heatmap(Heatmap.tissue.Q[,2:ncol(Heatmap.tissue.Q)], name = "Tissues.NES", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = F, cluster_columns = F, cluster_rows = F,
                             bottom_annotation = HeatmapAnnotation(Location = c("Total", "Total", "Total"), Tissues1 = c("Needle", "Root", "Root"), Tissues2 = c("Bud", "Needle", "Bud"), col = list(Tissues1 = c("Needle"="#1B839E", "Root"="#D95F02"), Tissues2 = c("Bud"="#7570B3", "Needle"="#1B839E"), Location = c("Total"="darkgoldenrod"))),
                             top_annotation = HeatmapAnnotation(Type = c(rep("Quantity", 3)), col=list(Type = c("Quantity" = "#374E55FF"))))
ht.S.tissue <- Heatmap(Heatmap.tissue.S[,2:ncol(Heatmap.tissue.S)], name = "Tissues.S", col = col_significance, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = F, cluster_columns = F,  cluster_rows = F,
                       bottom_annotation = HeatmapAnnotation(Location = c("Total", "Total", "Total"), Tissues1 = c("Needle", "Root", "Root"), Tissues2 = c("Bud", "Needle", "Bud"), Location = c("Total", "Total", "Total"), col = list(Tissues1 = c("Needle"="#1B839E", "Root"="#D95F02"), Tissues2 = c("Bud"="#7570B3", "Needle"="#1B839E"), Location = c("Total"="darkgoldenrod"))),
                       top_annotation = HeatmapAnnotation(Type = c(rep("Significance", 3)), col=list(Type = c("Significance" = "#DF8F44FF"))))

ht.S.tissue <- Heatmap(Heatmap.tissue.S[,2:ncol(Heatmap.tissue.S)], name = "Tissues.S", col = col_significance, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = F, cluster_columns = F,  cluster_rows = F,
                       bottom_annotation = HeatmapAnnotation(Tissues1 = c("Needle", "Root", "Root"), Tissues2 = c("Bud", "Needle", "Buds"), col = list(Tissues1 = c("Needle"="#1B839E", "Root"="#D95F02"), Treatment = c("Bud"="#7570B3", "Needle"="#1B839E"))),
                       top_annotation = HeatmapAnnotation(Type = c(rep("Significance", 3)), col=list(Type = c("Significance" = "#DF8F44FF"))))

ht.tissue <- draw(ht.Q.tissue + ht.S.tissue, ht_gap = unit(1, "cm"))

## Stress
ORDER.jic <- Tissues.res$Needle_Buds$pathway
Stress.res[["Total_FU"]] <- read.delim("clipboard")
Stress.res[["Total_UV"]] <- read.delim("clipboard")
Stress.res[["Total_Heat"]] <- read.delim("clipboard")
Stress.res[["Nucleus_UV"]] <- read.delim("clipboard")
Stress.res[["Nucleus_Heat"]] <- read.delim("clipboard")
Stress.res[["Chloroplast_UV"]] <- read.delim("clipboard")
Stress.res[["Chloroplast_Heat"]] <- read.delim("clipboard")
Stress.res[["Recovery_Heat"]] <- read.delim("clipboard")
Stress.res[["Recovery_UV"]] <- read.delim("clipboard")
Stress.res[["Total_Heat_UV"]] <- read.delim("clipboard")
Stress.res[["Nucleus_Heat_UV"]] <- read.delim("clipboard")
Stress.res[["Chloroplast_Heat_UV"]] <- read.delim("clipboard")

Heatmap.Stress.Q <- data.frame(Factor = Tissues.res$Needle_Buds$pathway, Total_FU = NA, Total_UV = NA, Total_Heat = NA,
                               Nucleus_UV = NA, Nucleus_Heat = NA, Chloroplast_UV = NA, Chloroplast_Heat = NA,
                               Recovery_Heat = NA, Recovery_UV = NA, Total_Heat_UV = NA, Nucleus_Heat_UV = NA, Chloroplast_Heat_UV = NA)
Heatmap.Stress.S <- data.frame(Factor = Tissues.res$Needle_Buds$pathway, Total_FU = NA, Total_UV = NA, Total_Heat = NA,
                               Nucleus_UV = NA, Nucleus_Heat = NA, Chloroplast_UV = NA, Chloroplast_Heat = NA,
                               Recovery_Heat = NA, Recovery_UV = NA, Total_Heat_UV = NA, Nucleus_Heat_UV = NA, Chloroplast_Heat_UV = NA)

for(i in 1:nrow(Heatmap.Stress.Q)){
  for(j in colnames(Heatmap.Stress.Q)[-1]){
    Heatmap.Stress.Q[i,j] <- Stress.res[[j]]$NES[which(Stress.res[[j]]$pathway == Heatmap.Stress.Q$Factor[i])]
  }
}
for(i in 1:nrow(Heatmap.Stress.S)){
  for(j in colnames(Heatmap.Stress.S)[-1]){
    Heatmap.Stress.S[i,j] <- Stress.res[[j]]$padj[which(Stress.res[[j]]$pathway == Heatmap.Stress.Q$Factor[i])]
    if(Heatmap.Stress.S[i,j] <= 0.1 & !is.na(Heatmap.Stress.S[i,j])){
      Heatmap.Stress.S[i,j] <- 0
    }else if(Heatmap.Stress.S[i,j] > 0.1  & !is.na(Heatmap.Stress.S[i,j])){
      Heatmap.Stress.S[i,j] <- 1
    }else if(is.na(Heatmap.Stress.S[i,j])){
      Heatmap.Stress.S[i,j] <- 1
      }
  }
}

rownames(Heatmap.Stress.Q) <- Heatmap.Stress.Q$Factor
rownames(Heatmap.Stress.S) <- Heatmap.Stress.S$Factor

ht.Q.stress <- Heatmap(Heatmap.Stress.Q[,2:ncol(Heatmap.Stress.Q)], name = "Stress.NES", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = F, cluster_columns = F, cluster_rows = F,
                       bottom_annotation = HeatmapAnnotation(Location = c("Total", "Total", "Total", "Nucleus", "Nucleus", "Chloroplast", "Chloroplast", "Nucleus", "Total", "Total", "Nucleus", "Chloroplast"), Stress1 = c("FU", "UV", "Heat", "UV", "Heat", "UV", "Heat", "Recovery", "UV", "UV", "UV", "UV"), Stress2 = c("Control", "Control", "Control", "Control", "Control", "Control", "Control", "Heat", "Recovery", "Heat", "Heat", "Heat"), col = list(Location = c("Total"="darkgoldenrod", "Nucleus"="darkmagenta", "Chloroplast"="darkgreen" ), Stress1 = c("FU"="#7E6148FF", "UV"="#3C5488FF", "Heat"="#E64B35FF", "Recovery"="#4DBBD5FF"), Stress2 = c("Control"="#00A087FF", "Heat"="#E64B35FF", "Recovery"="#4DBBD5FF"))),
                       top_annotation = HeatmapAnnotation(Type = c(rep("Quantity", 12)), col=list(Type = c("Quantity" = "#374E55FF"))))
ht.S.stress <- Heatmap(Heatmap.Stress.S[,2:ncol(Heatmap.Stress.S)], name = "Stress.S", col = col_significance, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = F, cluster_columns = F, cluster_rows = F,
                       bottom_annotation = HeatmapAnnotation(Location = c("Total", "Total", "Total", "Nucleus", "Nucleus", "Chloroplast", "Chloroplast", "Nucleus", "Total", "Total", "Nucleus", "Chloroplast"), Stress1 = c("FU", "UV", "Heat", "UV", "Heat", "UV", "Heat", "Recovery", "UV", "UV", "UV", "UV"), Stress2 = c("Control", "Control", "Control", "Control", "Control", "Control", "Control", "Heat", "Recovery", "Heat", "Heat", "Heat"), col = list(Location = c("Total"="darkgoldenrod", "Nucleus"="darkmagenta", "Chloroplast"="darkgreen" ), Stress1 = c("FU"="#7E6148FF", "UV"="#3C5488FF", "Heat"="#E64B35FF", "Recovery"="#4DBBD5FF"), Stress2 = c("Control"="#00A087FF", "Heat"="#E64B35FF", "Recovery"="#4DBBD5FF"))),
                       top_annotation = HeatmapAnnotation(Type = c(rep("Significance", 12)), col=list(Type = c("Significance" = "#DF8F44FF"))))


ht.S.stress <- Heatmap(Heatmap.tissue.S[,2:ncol(Heatmap.tissue.S)], name = "Tissues.S", col = col_significance, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = F, cluster_columns = F,  cluster_rows = F,
                       bottom_annotation = HeatmapAnnotation(Tissues1 = c("Needle", "Root", "Root"), Tissues2 = c("Bud", "Needle", "Bud"), col = list(Tissues1 = c("Needle"="#1B839E", "Root"="#D95F02"), Tissues2 = c("Bud"="#7570B3", "Needle"="#1B839E"))),
                       top_annotation = HeatmapAnnotation(Type = c(rep("Significance", 3)), col=list(Type = c("Significance" = "#DF8F44FF"))))

ht.S.tissue <- Heatmap(Heatmap.tissue.S[,2:ncol(Heatmap.tissue.S)], name = "Tissues.S", col = col_significance, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = F, cluster_columns = F,  cluster_rows = F,
                       bottom_annotation = HeatmapAnnotation(Tissues1 = c("Needle", "Root", "Root"), Tissues2 = c("Bud", "Needle", "Buds"), col = list(Tissues1 = c("Needle"="#1B839E", "Root"="#D95F02"), Treatment = c("Bud"="#7570B3", "Needle"="#1B839E"))),
                       top_annotation = HeatmapAnnotation(Type = c(rep("Significance", 3)), col=list(Type = c("Significance" = "#DF8F44FF"))))

ht.stress <- draw(ht.Q.stress + ht.S.stress, ht_gap = unit(1, "cm"))

draw(ht.Q.tissue + ht.S.tissue + ht.Q.stress + ht.S.stress, ht_gap = unit(1, "cm"))

#### Volcanos ####

Tissues <- list(Needle_Buds = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Tissues/Tissues_log10.zscore_Needle-Buds_all_Tissues.txt", header = T, sep = "\t", row.names = 1),
                Root_Needle = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Tissues/Tissues_log10.zscore_Root-Needle_all_Tissues.txt", header = T, sep = "\t", row.names = 1),
                Root_Buds = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Tissues/Tissues_log10.zscore_Root-Buds_all_Tissues.txt", header = T, sep = "\t", row.names = 1))

Stress <- list(Total_FU = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Fusarium/Fusarium_log10.zscore_FU-Control_all_Fusarium.txt", header = T, sep = "\t", row.names = 1),
               Total_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/UVTotal/UVTotal_log10.zscore_UV-Control_all_UVTotal.txt", header = T, sep = "\t", row.names = 1),
               Total_Heat = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/HeatTotal/HeatTotal_log10.zscore_Heat-Control_all_HeatTotal.txt", header = T, sep = "\t", row.names = 1),
               Nucleus_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/UVNucleus/UVNucleus_log10.zscore_UV-Control_all_UVNucleus.txt", header = T, sep = "\t", row.names = 1),
               Nucleus_Heat = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/HeatNucleus/HeatNucleus_log10.zscore_Heat-Control_all_HeatNucleus.txt", header = T, sep = "\t", row.names = 1),
               Chloroplast_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/UVChloroplast/UVChloroplast_log10.zscore_UV-Control_all_UVChloroplast.txt", header = T, sep = "\t", row.names = 1),
               Chloroplast_Heat = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/HeatChloroplast/HeatChloroplast_log10.zscore_Heat-Control_all_HeatChloroplast.txt", header = T, sep = "\t", row.names = 1),
               Recovery_Heat = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/HeatNucleus/HeatNucleus_log10.zscore_Recovery-Heat_all_HeatNucleus.txt", header = T, sep = "\t", row.names = 1),
               Recovery_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/UVTotal/UVTotal_log10.zscore_UV-Recovery_all_UVTotal.txt", header = T, sep = "\t", row.names = 1),
               Total_Heat_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Heat_UV_Total/Heat_UV_Total_log10.zscore_UV-Heat_all_Heat_UV_Total.txt", header = T, sep = "\t", row.names = 1),
               Nucleus_Heat_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Heat_UV_Nucleus/Heat_UV_Nucleus_log10.zscore_UV-Heat_all_Heat_UV_Nucleus.txt", header = T, sep = "\t", row.names = 1),
               Chloroplast_Heat_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Heat_UV_Chloroplast/Heat_UV_Chloroplast_log10.zscore_UV-Heat_all_Heat_UV_Chloroplast.txt", header = T, sep = "\t", row.names = 1))

outputs <- "F:/ATLAS_Pra/Protein_module/0.Process/Figures/Descriptives/Proteins/Volcano/"
pal <- c( "darkgray", "lightgoldenrod3", "lightsteelblue", "lightskyblue4")
output.path <- "F:/ATLAS_Pra/Protein_module/0.Process/Figures/Descriptives/Proteins/Volcano/"

for(i in 1:length(Tissues)){
  Tissues[[i]]$Symbol <- PROTAS.annotation[rownames(Tissues[[i]]),4]
  Tissues[[i]]$logpvalue <- -log10(Tissues[[i]]$adj.P.Val)
  Tissues[[i]]$group <- "NotSignificant"
  Tissues[[i]]$group[which(Tissues[[i]]$adj.P.Val <= 0.05 & abs(Tissues[[i]]$logFC) < 1.5 )] <- "Significant"
  Tissues[[i]]$group[which(Tissues[[i]]$adj.P.Val > 0.05 & abs(Tissues[[i]]$logFC) >= 1.5 )] <- "FoldChange"
  Tissues[[i]]$group[which(Tissues[[i]]$adj.P.Val <= 0.05 & abs(Tissues[[i]]$logFC) >= 1.5 )] <- "Significant&FoldChange"
  top_peaks_down <- Tissues[[i]][with(Tissues[[i]], order(logFC, adj.P.Val, decreasing = FALSE)),][1:5,]
  top_peaks_top <- Tissues[[i]][with(Tissues[[i]], order(logFC, adj.P.Val, decreasing = TRUE)),][1:5,]
  top_peaks <- rbind(top_peaks_top, top_peaks_down)
  a <- list()
  for(V in seq_len(nrow(top_peaks))){
    m <- top_peaks[V, ]
    a[[V]] <- list(
      x = m[["logFC"]],
      y = m[["logpvalue"]],
      text = m[["Symbol"]],
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 0.3,
      ax = -10,
      ay = 50,
      size = 0.5 
    )
  }
  p1 <- plot_ly(data = Tissues[[i]], x = Tissues[[i]]$logFC, y = Tissues[[i]]$logpvalue, text = Tissues[[i]]$Symbol, mode = "markers", color = Tissues[[i]]$group, colors = pal, xaxis = "FoldChange[log2]",  yaxis = "p-value[-log10]") %>% 
    layout(title = names(Tissues[i])) %>% 
    layout(annotations = a, yaxis = list(range = c(0, 80)), xaxis = list(range = c(-10, 10)))
  write.table(top_peaks, file = paste0(output.path, names(Tissues[i]), "_TopPeaks_2volcano.txt"), quote = F, sep = "\t")
  save_image(p = p1, file = paste0(output.path, names(Tissues[i]), "_2volcano.pdf"))
}

for(i in 1:length(Stress)){
  Stress[[i]]$Symbol <- PROTAS.annotation[rownames(Stress[[i]]),4]
  Stress[[i]]$logpvalue <- -log10(Stress[[i]]$adj.P.Val)
  Stress[[i]]$group <- "NotSignificant"
  Stress[[i]]$group[which(Stress[[i]]$adj.P.Val <= 0.05 & abs(Stress[[i]]$logFC) < 1.5 )] <- "Significant"
  Stress[[i]]$group[which(Stress[[i]]$adj.P.Val > 0.05 & abs(Stress[[i]]$logFC) >= 1.5 )] <- "FoldChange"
  Stress[[i]]$group[which(Stress[[i]]$adj.P.Val <= 0.05 & abs(Stress[[i]]$logFC) >= 1.5 )] <- "Significant&FoldChange"
  top_peaks_down <- Stress[[i]][with(Stress[[i]], order(logFC, adj.P.Val, decreasing = FALSE)),][1:5,]
  top_peaks_top <- Stress[[i]][with(Stress[[i]], order(logFC, adj.P.Val, decreasing = TRUE)),][1:5,]
  top_peaks <- rbind(top_peaks_top, top_peaks_down)
  a <- list()
  for(V in seq_len(nrow(top_peaks))){
    m <- top_peaks[V, ]
    a[[V]] <- list(
      x = m[["logFC"]],
      y = m[["logpvalue"]],
      text = m[["Symbol"]],
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 0.3,
      ax = -10,
      ay = 50,
      size = 0.5 
    )
  }
  p1 <- plot_ly(data = Stress[[i]], x = Stress[[i]]$logFC, y = Stress[[i]]$logpvalue, text = Stress[[i]]$Symbol, mode = "markers", color = Stress[[i]]$group, colors = pal, xaxis = "FoldChange[log2]",  yaxis = "p-value[-log10]") %>% 
    layout(title = names(Stress[i])) %>% 
    layout(annotations = a, yaxis = list(range = c(0, 80)), xaxis = list(range = c(-10, 10)))
  write.table(top_peaks, file = paste0(output.path, names(Stress[i]), "_TopPeaks_2volcano.txt"), quote = F, sep = "\t")
  save_image(p = p1, file = paste0(output.path, names(Stress[i]), "_2volcano.pdf"))
}

#### Upset of differential prots ####

## Tissues

Tissues <- list(Needle_Buds = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Tissues/Tissues_log10.zscore_Needle-Buds_filtered_Tissues.txt", header = T, sep = "\t", row.names = 1),
                Root_Needle = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Tissues/Tissues_log10.zscore_Root-Needle_filtered_Tissues.txt", header = T, sep = "\t", row.names = 1),
                Root_Buds = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Tissues/Tissues_log10.zscore_Root-Buds_filtered_Tissues.txt", header = T, sep = "\t", row.names = 1))

Tissues <- lapply(X = Tissues, FUN = function(x){
  rownames(x)
})

upset(fromList(Tissues), sets = names(Tissues), order.by = "freq", 
      keep.order = TRUE, empty.intersections = T, mb.ratio = c(0.6,0.4), sets.bar.color = "navajowhite3")

## Stress

Stress <- list(Total_FU = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Fusarium/Fusarium_log10.zscore_FU-Control_filtered_Fusarium.txt", header = T, sep = "\t", row.names = 1),
              Total_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/UVTotal/UVTotal_log10.zscore_UV-Control_filtered_UVTotal.txt", header = T, sep = "\t", row.names = 1),
              Total_Heat = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/HeatTotal/HeatTotal_log10.zscore_Heat-Control_filtered_HeatTotal.txt", header = T, sep = "\t", row.names = 1),
              Nucleus_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/UVNucleus/UVNucleus_log10.zscore_UV-Control_filtered_UVNucleus.txt", header = T, sep = "\t", row.names = 1),
              Nucleus_Heat = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/HeatNucleus/HeatNucleus_log10.zscore_Heat-Control_filtered_HeatNucleus.txt", header = T, sep = "\t", row.names = 1),
              Chloroplast_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/UVChloroplast/UVChloroplast_log10.zscore_UV-Control_filtered_UVChloroplast.txt", header = T, sep = "\t", row.names = 1),
              Chloroplast_Heat = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/HeatChloroplast/HeatChloroplast_log10.zscore_Heat-Control_filtered_HeatChloroplast.txt", header = T, sep = "\t", row.names = 1),
              Recovery_Heat = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/HeatNucleus/HeatNucleus_log10.zscore_Recovery-Heat_filtered_HeatNucleus.txt", header = T, sep = "\t", row.names = 1),
              Recovery_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/UVTotal/UVTotal_log10.zscore_UV-Recovery_filtered_UVTotal.txt", header = T, sep = "\t", row.names = 1),
              Total_Heat_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Heat_UV_Total/Heat_UV_Total_log10.zscore_UV-Heat_filtered_Heat_UV_Total.txt", header = T, sep = "\t", row.names = 1),
              Nucleus_Heat_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Heat_UV_Nucleus/Heat_UV_Nucleus_log10.zscore_UV-Heat_filtered_Heat_UV_Nucleus.txt", header = T, sep = "\t", row.names = 1),
              Chloroplast_Heat_UV = read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Tables/sva+limma/Heat_UV_Chloroplast/Heat_UV_Chloroplast_log10.zscore_UV-Heat_filtered_Heat_UV_Chloroplast.txt", header = T, sep = "\t", row.names = 1))

Stress <- lapply(X = Stress, FUN = function(x){
  rownames(x)
})

upset(fromList(Stress), sets = names(Stress), order.by = "freq", 
      keep.order = TRUE, empty.intersections = T, mb.ratio = c(0.6,0.4), sets.bar.color = "navajowhite3", )

#### Network Analyses + Module corr + Module Enrich Heatmaps ####

options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 15)
output.path <- "F:/ATLAS_Pra/Protein_module/0.Process/Figures/Descriptives/Proteins/WGCNA/"
setwd(output.path)

## Data Prepare and checks

data0 <- read.delim("F:/ATLAS_Pra/Protein_module/0.Process/ProtAtlas.Analyses.txt", header = T, sep = "\t")
data0 <- t(data0)
colnames(data0) <- data0[1,]
data0 <- data0[-c(1:2),]
ROWN <- rownames(data0)
data0 <- apply(data0, MARGIN = 2, FUN = as.numeric)
rownames(data0) <- ROWN

data0.zscore <- apply(data0, MARGIN = 2, FUN = function(x){
  x <- (x - mean(x))/sd(x)
  x
})

data0.log10 <- apply(data0, MARGIN = 2, FUN = function(x){
  x <- log10(x + 1.1)
  x
})

data0.log10.zscore <- apply(data0, MARGIN = 2, FUN = function(x){
  x <- log10(x + 1)
  x <- (x - mean(x))/sd(x)
  x
})

data0.quantiles <- apply(data0, MARGIN = 2, FUN = function(x){
  x <- normalizeQuantiles(x, ties = TRUE)
  x
})

data0.raw <- data0

data0.list <- list(raw = data0.raw, log10 = data0.log10, log10.zscore = data0.log10.zscore, quantiles = data0.quantiles)
sample0.metadata <- data.frame(sample_ID = rownames(data0), Stress = labelitass$Stress, Intensity = labelitass$Intensity,
                               Population = labelitass$Population, Technique = labelitass$Technique, Tissues = labelitass$Tissues,
                               Stress.Intensity = labelitass$Stress.Intensity, Population.Stress = labelitass$Population.Stress,
                               Population.Intensity = labelitass$Population.Intensity, Technique.Stress = labelitass$Technique.Stress,
                               Technique.Intensity = labelitass$Technique.Intensity
                               ) # labelitass from UMAP/PCA section

match(sample0.metadata$sample_ID, rownames(data0))

## check outliers: NO OUTLIERS; log10 best transformation

goodSamplesGenes(datExpr = data0.list$log10) # ALL TRUE/OK

sampleTree = hclust(dist(data0.list$log10), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

sampleTree2 = hclust(dist(data0.list$log10), method = "average")
sample0.metadata.factor <- sample0.metadata 
sample0.metadata.factor[,2:ncol(sample0.metadata.factor)] <- apply(sample0.metadata.factor[,2:ncol(sample0.metadata.factor)], MARGIN = 2, FUN = function(x){
  as.numeric(as.factor(x))
})
rownames(sample0.metadata.factor) <- sample0.metadata.factor$sample_ID
traitColors = numbers2colors(sample0.metadata.factor[,-1], signed = FALSE)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(sample0.metadata.factor[,-1]), 
                    main = "Sample dendrogram and trait heatmap")

## Choose soft threshold parameter

powers = c(c(1:100), seq(from = 22, to=100, by=2))
sft = pickSoftThreshold(data0.list$log10, powerVector = powers, verbose = 5, dataIsExpr = TRUE, networkType = "signed", blockSize = 30) 

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

# REASONS: proteomic data does not reach scale-free threshold:
# - scale free networks in biology are rare https://www.nature.com/articles/s41467-019-08746-5: not the best thing for thresholding
# - data could exhibit a strong driver that makes a subset of the samples globally different from the rest
# causing high correlation among large groups of genes which invalidates the assumtion of the scale-free topology approximation
# In this case, the lack of scale-free topology fit turns out to be caused by interesting biological factors (as showed in UMAP and OUTLIERS in SampleOutliersTraits)
# so appropiate sof-thresholding power can be chosen based on the number of samples specified by a table update in December 2017.
# Less than 20;20-30;30-40;>40 Unsigned and hybrid: 9;8;7;6
# Less than 20;20-30;30-40;>40 Signed: 18;16;14;12
# CONCLUSION: Signed hybrid power 7
# CONCLUSION: Signed power 12

softPower.hybrid <- 7
softPower.signed <- 12

## Turning data into topological overlap matrix

adjacency.hybrid = adjacency(datExpr = data0.list$log10, power = softPower.hybrid, type = "signed hybrid", corFnc = "bicor")
adjacency.signed = adjacency(datExpr = data0.list$log10, power = softPower.signed, type = "signed", corFnc = "bicor")

TOM.hybrid = TOMsimilarity(adjacency.hybrid)
TOM.signed = TOMsimilarity(adjacency.signed)

dissTOM.hybrid = 1 - TOM.hybrid
dissTOM.signed = 1 - TOM.signed

## Grouping genes in modules

geneTree.hybrid = hclust(as.dist(dissTOM.hybrid), method = "average")
geneTree.signed = hclust(as.dist(dissTOM.signed), method = "average");

dynamicMods.hybrid = cutreeDynamic(dendro = geneTree.hybrid, distM = dissTOM.hybrid, deepSplit = 2, 
                            pamRespectsDendro = FALSE, minClusterSize = 30)
dynamicMods.signed = cutreeDynamic(dendro = geneTree.signed, distM = dissTOM.signed, deepSplit = 2, 
                                   pamRespectsDendro = FALSE, minClusterSize = 30)

table(dynamicMods.hybrid)
length(table(dynamicMods.hybrid)) 

table(dynamicMods.signed)
length(table(dynamicMods.signed)) 

dynamicColors.hybrid = labels2colors(dynamicMods.hybrid)
table(dynamicColors.hybrid)

dynamicColors.signed = labels2colors(dynamicMods.signed)
table(dynamicColors.signed)

plotDendroAndColors(geneTree.hybrid, dynamicColors.hybrid, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
plotDendroAndColors(geneTree.signed, dynamicColors.signed, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")

## Merge modules

MEList.hybrid = moduleEigengenes(data0.list$log10, colors = dynamicColors.hybrid)
MEs.hybrid = MEList.hybrid$eigengenes
MEDiss.hybrid = 1-cor(MEs.hybrid)
METree.hybrid = hclust(as.dist(MEDiss.hybrid), method = "average")
plot(METree.hybrid, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres.hybrid = 0.3
abline(h=MEDissThres.hybrid, col = "red")
merge.hybrid = mergeCloseModules(data0.list$log10, dynamicColors.hybrid, cutHeight = MEDissThres.hybrid, verbose = 3) 
mergedColors.hybrid = merge.hybrid$colors  
mergedMEs.hybrid = merge.hybrid$newMEs  
plotDendroAndColors(geneTree.hybrid, cbind(dynamicColors.hybrid, mergedColors.hybrid), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05) 
write.table(merge.hybrid$oldMEs,file="Hybrid.oldMEs.txt", quote = F, sep = "\t")
write.table(merge.hybrid$newMEs,file="Hybrid.newMEs.txt", quote = F, sep = "\t")

MEList.signed = moduleEigengenes(data0.list$log10, colors = dynamicColors.signed)
MEs.signed = MEList.signed$eigengenes
MEDiss.signed = 1-cor(MEs.signed)
METree.signed = hclust(as.dist(MEDiss.signed), method = "average")
plot(METree.signed, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres.signed = 0.3
abline(h=MEDissThres.signed, col = "red")
merge.signed = mergeCloseModules(data0.list$log10, dynamicColors.signed, cutHeight = MEDissThres.signed, verbose = 3) 
mergedColors.signed = merge.signed$colors  
mergedMEs.signed = merge.signed$newMEs  
plotDendroAndColors(geneTree.signed, cbind(dynamicColors.signed, mergedColors.signed), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05) 
write.table(merge.signed$oldMEs,file="Signed.oldMEs.txt", quote = F, sep = "\t")
write.table(merge.signed$newMEs,file="Signed.newMEs.txt", quote = F, sep = "\t")

moduleColors.hybrid = mergedColors.hybrid
moduleColors.signed = mergedColors.signed

colorOrder = c("grey", standardColors(50))

moduleLabels.hybrid = match(moduleColors.hybrid, colorOrder)-1
MEs.hybrid = mergedMEs.hybrid

moduleLabels.signed = match(moduleColors.signed, colorOrder)-1
MEs.signed = mergedMEs.signed

## Associating modules and phenotypes

## hybrid (only this one because signed is not the main purpose)

sample0.metadata$Technique.Tissue.Stress.Intensity <- paste0(sample0.metadata$Technique, ".", sample0.metadata$Tissues, ".", sample0.metadata$Stress, ".", sample0.metadata$Intensity) 
datTrait <- binarizeCategoricalColumns(sample0.metadata[,ncol(sample0.metadata)], includePairwise = F, includeLevelVsAll = T, dropFirstLevelVsAll = F)
datTrait <- binarizeCategoricalColumns(sample0.metadata, includePairwise = F, includeLevelVsAll = T, dropFirstLevelVsAll = F)

nGenes = ncol(data0.list$log10)
nSamples = nrow(data0.list$log10)

MEs0.hybrid = moduleEigengenes(data0.list$log10, moduleColors.hybrid)$eigengenes
MEs.hybrid = orderMEs(MEs0.hybrid)
moduleTraitCor.hybrid = cor(MEs.hybrid, datTrait, use = "p")
moduleTraitPvalue.hybrid = corPvalueStudent(moduleTraitCor.hybrid, nSamples)
BH=0.05/25

moduleTraitCor.hybrid.jc <- moduleTraitCor.hybrid
moduleTraitPvalue.hybrid.jc <- moduleTraitPvalue.hybrid

for(i in 1:nrow(moduleTraitCor.hybrid)){
  for(j in 1:ncol(moduleTraitCor.hybrid)){
    if(moduleTraitPvalue.hybrid[i,j] > BH){moduleTraitCor.hybrid[i,j] <- 0}
  }
}
moduleTraitCor.hybrid <- moduleTraitCor.hybrid[,!(colSums(moduleTraitCor.hybrid) == 0)]
moduleTraitPvalue.hybrid <- moduleTraitPvalue.hybrid[,colnames(moduleTraitCor.hybrid)]

textMatrix.hybrid =  paste(signif(moduleTraitCor.hybrid, 2), "\n(",
                    signif(moduleTraitPvalue.hybrid, 1), ")", sep = "")
dim(textMatrix.hybrid) = dim(moduleTraitCor.hybrid)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor.hybrid,
               xLabels = colnames(moduleTraitCor.hybrid),
               yLabels = names(MEs.hybrid),
               ySymbols = names(MEs.hybrid),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix.hybrid,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


sort(table(mergedColors.hybrid), decreasing = T)
rownames(moduleTraitCor.hybrid.jc) == colnames(MEs.hybrid)
rownames(moduleTraitCor.hybrid) <- c("M09(69)", "M07(107)", "M04(417)", "M05(153)", "M03(640)", "M02(1127)", "M11(34)", "M01(4870)", "M06(133)", "M10(54)", "M08(78)",
                                     "M12(15)")
colnames(MEs.hybrid) <- rownames(moduleTraitCor.hybrid)

write.table(moduleTraitCor.hybrid, file = "F:/ATLAS_Pra/Protein_module/0.Process/Figures/Descriptives/Proteins/WGCNA/moduleTraitCor.hybrid.filtered.renamed.txt", quote = F, sep = "\t", row.names = T, col.names = T)
write.table(MEs.hybrid, file = "F:/ATLAS_Pra/Protein_module/0.Process/Figures/Descriptives/Proteins/WGCNA/MEs.hybrid.filtered.renamed.txt", quote = F, sep = "\t", row.names = T, col.names = T)

## plotting final TraitCor heatmap

col_trait <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

ht.Trait.hor <- Heatmap(t(moduleTraitCor.hybrid), name = "Trait.cor", col = col_fun, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F,
                      right_annotation = rowAnnotation(Location = c(rep("Chloroplast", 7), rep("Nucleus", 8), rep("Total", 10)),  Tissue = c(rep("Needle", 15), rep("Bud", 2), rep("Needle", 7), rep("Root",1)), 
                                                           Stress = c("Control", rep("Heat", 3), rep("UV", 3), "Control", rep("Heat", 4), "Recovery", rep("UV", 2), "Control", "FU", "Control", "Heat", "Heat", "Recovery", rep("UV", 3), "Control"), 
                                                           Intensity = c("T0", "T1", "T2", "T3", "T1", "T2", "T3", "T0", "T1", "T4", "T2", "T3", "T0", "T1", "T2", "T0", "T1", "T0", "T1", "T2", "T0", "T1", "T2", "T3", "T0" ),
                                                           col = list(Location = c("Total"="darkgoldenrod", "Nucleus"="darkmagenta", "Chloroplast"="darkgreen" ),
                                                                      Tissue = c("Root"="#D95F02", "Bud"="#7570B3", "Needle"="#1B839E"),
                                                                      Stress = c("Control"="#00A087FF", "Recovery"="#4DBBD5FF", "FU"="#7E6148FF", "UV"="#3C5488FF", "Heat"="#E64B35FF"),
                                                                      Intensity = c("T0"="gray88", "T1"="gray68", "T2"="gray48", "T3"="gray28", "T4"="gray18", "T5"="gray8")) 
                                                           ),
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      if(t(moduleTraitCor.hybrid)[i, j] > 0 | t(moduleTraitCor.hybrid)[i, j] < 0)
                        grid.text(sprintf("%.1f", t(moduleTraitCor.hybrid)[i, j]), x, y, gp = gpar(fontsize = 10))
                    })
ht.Trait.ver <- Heatmap(moduleTraitCor.hybrid, name = "Trait.cor", col = col_fun, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F,
                    bottom_annotation = HeatmapAnnotation(Location = c(rep("Chloroplast", 7), rep("Nucleus", 8), rep("Total", 10)),  Tissue = c(rep("Needle", 15), rep("Bud", 2), rep("Needle", 7), rep("Root",1)), 
                                                     Stress = c("Control", rep("Heat", 3), rep("UV", 3), "Control", rep("Heat", 4), "Recovery", rep("UV", 2), "Control", "FU", "Control", "Heat", "Heat", "Recovery", rep("UV", 3), "Control"), 
                                                     Intensity = c("T0", "T1", "T2", "T3", "T1", "T2", "T3", "T0", "T1", "T4", "T2", "T3", "T0", "T1", "T2", "T0", "T1", "T0", "T1", "T2", "T0", "T1", "T2", "T3", "T0" ),
                                                     col = list(Location = c("Total"="darkgoldenrod", "Nucleus"="darkmagenta", "Chloroplast"="darkgreen" ),
                                                                Tissue = c("Root"="#D95F02", "Bud"="#7570B3", "Needle"="#1B839E"),
                                                                Stress = c("Control"="#00A087FF", "Recovery"="#4DBBD5FF", "FU"="#7E6148FF", "UV"="#3C5488FF", "Heat"="#E64B35FF"),
                                                                Intensity = c("T0"="gray88", "T1"="gray68", "T2"="gray48", "T3"="gray28", "T4"="gray18", "T5"="gray8")) 
                    ),
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      if(moduleTraitCor.hybrid[i, j] > 0 | moduleTraitCor.hybrid[i, j] < 0)
                        grid.text(sprintf("%.1f", moduleTraitCor.hybrid[i, j]), x, y, gp = gpar(fontsize = 10))
                    })

## Enrichment based on module membership

MEs.hybrid <- read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Figures/Descriptives/Proteins/WGCNA/MEs.hybrid.filtered.renamed.txt")
geneModuleMembership = as.data.frame(cor(data0.list$log10, MEs.hybrid, use = "p"))

Membership.res <- list()
for(i in colnames(geneModuleMembership)){
  data <- data.frame(ID = rownames(geneModuleMembership), logFC = geneModuleMembership[,i])
  data <- deframe(data)
  fGseaRES <- fgsea(pathways = Bin_db, stats = data, minSize = 5, maxSize = 1000, nPermSimple = 100000)
  fgseaResTidy <- fGseaRES %>%
    as_tibble() %>%
    arrange(desc(NES))
  write.table(fgseaResTidy[,-ncol(fgseaResTidy)], file = paste0(output.path, i, "_fgsea.txt"), col.names = T, sep = "\t", quote = F)
  Membership.res[[i]] <- fgseaResTidy 
}


Heatmap.Modules.Q <- data.frame(Factor = Tissues.res$Needle_Buds$pathway, 
                               "M09(69)" = NA, "M07(107)" = NA, "M04(417)" = NA,
                               "M05(153)" = NA, "M03(640)" = NA, "M02(1127)" = NA, "M11(34)" = NA,
                               "M01(4870)" = NA, "M06(133)" = NA, "M10(54)" = NA, "M08(78)" = NA, "M12(15)" = NA)
Heatmap.Modules.S <- data.frame(Factor = Tissues.res$Needle_Buds$pathway, 
                                "M09(69)" = NA, "M07(107)" = NA, "M04(417)" = NA,
                                "M05(153)" = NA, "M03(640)" = NA, "M02(1127)" = NA, "M11(34)" = NA,
                                "M01(4870)" = NA, "M06(133)" = NA, "M10(54)" = NA, "M08(78)" = NA, "M12(15)" = NA)

for(i in 1:nrow(Heatmap.Modules.Q)){
  for(j in colnames(Heatmap.Modules.Q)[-1]){
    Heatmap.Modules.Q[i,j] <- Membership.res[[j]]$NES[which(Membership.res[[j]]$pathway == Heatmap.Modules.Q$Factor[i])]
  }
}
for(i in 1:nrow(Heatmap.Modules.S)){
  for(j in colnames(Heatmap.Modules.S)[-1]){
    Heatmap.Modules.S[i,j] <- Membership.res[[j]]$padj[which(Membership.res[[j]]$pathway == Heatmap.Modules.Q$Factor[i])]
    if(Heatmap.Modules.S[i,j] <= 0.1 & !is.na(Heatmap.Modules.S[i,j])){
      Heatmap.Modules.S[i,j] <- 0
    }else if(Heatmap.Modules.S[i,j] > 0.1  & !is.na(Heatmap.Modules.S[i,j])){
      Heatmap.Modules.S[i,j] <- 1
    }else if(is.na(Heatmap.Modules.S[i,j])){
      Heatmap.Modules.S[i,j] <- 1
    }
  }
}

rownames(Heatmap.Modules.Q) <- Heatmap.Modules.Q$Factor
rownames(Heatmap.Modules.S) <- Heatmap.Modules.S$Factor

ht.Q.Modules <- Heatmap(Heatmap.Modules.Q[,2:ncol(Heatmap.Modules.Q)], name = "Modules.NES", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F,
                       top_annotation = HeatmapAnnotation(Type = c(rep("Quantity", 12)), col=list(Type = c("Quantity" = "#374E55FF"))))
ht.S.Modules <- Heatmap(Heatmap.Modules.S[,2:ncol(Heatmap.Modules.S)], name = "Modules.S", col = col_significance, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F,
                       top_annotation = HeatmapAnnotation(Type = c(rep("Significance", 12)), col=list(Type = c("Significance" = "#DF8F44FF"))))


ht.S.stress <- Heatmap(Heatmap.tissue.S[,2:ncol(Heatmap.tissue.S)], name = "Tissues.S", col = col_significance, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = F, cluster_columns = F,  cluster_rows = F,
                       bottom_annotation = HeatmapAnnotation(Tissues1 = c("Needle", "Root", "Root"), Tissues2 = c("Bud", "Needle", "Bud"), col = list(Tissues1 = c("Needle"="#1B839E", "Root"="#D95F02"), Tissues2 = c("Bud"="#7570B3", "Needle"="#1B839E"))),
                       top_annotation = HeatmapAnnotation(Type = c(rep("Significance", 3)), col=list(Type = c("Significance" = "#DF8F44FF"))))

ht.S.tissue <- Heatmap(Heatmap.tissue.S[,2:ncol(Heatmap.tissue.S)], name = "Tissues.S", col = col_significance, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = F, cluster_columns = F,  cluster_rows = F,
                       bottom_annotation = HeatmapAnnotation(Tissues1 = c("Needle", "Root", "Root"), Tissues2 = c("Bud", "Needle", "Buds"), col = list(Tissues1 = c("Needle"="#1B839E", "Root"="#D95F02"), Treatment = c("Bud"="#7570B3", "Needle"="#1B839E"))),
                       top_annotation = HeatmapAnnotation(Type = c(rep("Significance", 3)), col=list(Type = c("Significance" = "#DF8F44FF"))))


draw(ht.Q.tissue + ht.S.tissue + ht.Q.stress + ht.S.stress + ht.Q.Modules + ht.S.Modules, ht_gap = unit(1, "cm"))

## Export to plot network gephi

modGenes = colnames(data0.list$log10)
modTOM = TOM.hybrid
dimnames(modTOM) = list(modGenes, modGenes)
modTOMSignificantes = which(modTOM>0.1)
geneInfo0 = data.frame(ESTs = modGenes,
                       moduleColor = moduleColors.hybrid
                       )
geneOrder = order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo.csv")
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = "CytoscapeEdgeFile07.txt",
                               nodeFile = "CytoscapeNodeFile07.txt",
                               weighted = TRUE,
                               threshold = 0.7,
                               nodeNames = modGenes,
                               nodeAttr = moduleColors.hybrid)

nodes <- read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Figures/Descriptives/Proteins/WGCNA/CytoscapeNodeFile01.txt", header = T)
edges <- read.delim("F:/ATLAS_Pra/Protein_module/0.Process/Figures/Descriptives/Proteins/WGCNA/CytoscapeEdgeFile01.txt", header = T)

nodes$ID <- rownames(nodes)
nodes <- nodes[,c(4,1,3)]
colnames(nodes) <- c("id", "label", "Module")
nodes$Module <- gsub("turquoise", "M01", x = nodes$Module, fixed = T)
nodes$Module <- gsub("blue", "M02", x = nodes$Module, fixed = T)
nodes$Module <- gsub("grey60", "M03", x = nodes$Module, fixed = T)
nodes$Module <- gsub("salmon", "M04", x = nodes$Module, fixed = T)
nodes$Module <- gsub("black", "M05", x = nodes$Module, fixed = T)
nodes$Module <- gsub("pink", "M06", x = nodes$Module, fixed = T)
nodes$Module <- gsub("magenta", "M07", x = nodes$Module, fixed = T)
nodes$Module <- gsub("purple", "M08", x = nodes$Module, fixed = T)
nodes$Module <- gsub("greenyellow", "M09", x = nodes$Module, fixed = T)
nodes$Module <- gsub("lightcyan", "M10", x = nodes$Module, fixed = T)
nodes$Module <- gsub("lightgreen", "M11", x = nodes$Module, fixed = T)
nodes$Module <- gsub("grey", "M12", x = nodes$Module, fixed = T)
sort(table(nodes$Module), decreasing = T)

nodes.Mod <- nodes[nodes$Module %in% c("M02", "M03", "M05", "M07"),]

edges <- edges[which(edges$fromNode %in% nodes.Mod$label & edges$toNode %in% nodes.Mod$label),]

nodes$id[which(nodes$label == "all-04-1000034")]
edges$source <- match(edges$fromNode, nodes$label)
edges$target <- match(edges$toNode, nodes$label)

edges <- edges[,c(7,8,1,2,3)]
colnames(edges) <- c("source", "target", "source_label", "target_label", "weight")

write.table(nodes.Mod, file = "F:/ATLAS_Pra/Protein_module/0.Process/Figures/Descriptives/Proteins/WGCNA/PROTAS01_M02357_nodes_gephi.csv", col.names = T, row.names = F, sep = ",")
write.table(edges, file = "F:/ATLAS_Pra/Protein_module/0.Process/Figures/Descriptives/Proteins/WGCNA/PROTAS01_M02357_edges_gephi.csv", col.names = T, row.names = F, sep = ",")

write.gexf(nodes = nodes[,1:2], edges = edges[,1:2], nodesAtt = nodes, edgesAtt = edges, 
           output = "F:/ATLAS_Pra/Protein_module/0.Process/Figures/Descriptives/Proteins/WGCNA/PROTAS04_gephi.gexf",
           defaultedgetype = "undirected")

nodes.Mod$Symbol <- NA
for(i in 1:nrow(nodes.Mod)){
  nodes.Mod$Symbol[i] <- PROTAS.annotation$Symbol[which(PROTAS.annotation$ProtID == nodes.Mod$label[i])]
}

write.table(nodes.Mod, file = "F:/ATLAS_Pra/Protein_module/0.Process/Figures/Descriptives/Proteins/WGCNA/PROTAS01_M02357_nodes_gephi_ANNOT.csv", col.names = T, row.names = F, sep = ",")

#### Evolutionary Proteomics ####

output.path <- "F:/ATLAS_Pra/Protein_module/0.Process/Figures/Descriptives/Proteins/EvoProteomics/"

## Annotate abundance tables with gene age

PROTAS <- read.delim(file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/Protein_module/ProtAtlas.Analyses.txt", header = T, sep = "\t")
Ages <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/4.GenEra/Output/Pinus_radiata_proteome/3347_gene_ages.tsv", header = T, sep = "\t")
Ages$X.gene <- do.call(rbind, strsplit(x = Ages$X.gene, split = ".", fixed = T))[,1] 
PROTAS$Phylostratum <- NA
PROTAS <- PROTAS[PROTAS$ProtID %in% Ages$X.gene,]

for(i in 1:nrow(PROTAS)){
  age <- unique(Ages$rank[which(Ages$X.gene == PROTAS$ProtID[i])])
  if(length(age) == 1){
    PROTAS$Phylostratum[i] <- age
  }else if(length(age) > 1){
    PROTAS$Phylostratum[i] <- min(age)
  }else if(length(age) == 0){
    stop("Something went wrong boy ... \n")
  }
}

PROTAS <- PROTAS[,c(ncol(PROTAS),1:ncol(PROTAS)-1)]
colnames(PROTAS)[2] <- "GeneID"
PROTAS <- PROTAS[,-3]
PROTAS <- PROTAS[!is.na(PROTAS$Phylostratum),]

## Prepare both transforms log10 and log10.zscore

EVOPROTAS.log10 <- PROTAS
EVOPROTAS.log10.zscore <- PROTAS

EVOPROTAS.log10[,3:ncol(EVOPROTAS.log10)] <- apply(EVOPROTAS.log10[,3:ncol(EVOPROTAS.log10)], MARGIN = 2, FUN = function(x){
  x <- log10(x + 1.1)
  x
})

EVOPROTAS.log10.zscore[,3:ncol(EVOPROTAS.log10.zscore)] <- apply(EVOPROTAS.log10.zscore[,3:ncol(EVOPROTAS.log10.zscore)], MARGIN = 2, FUN = function(x){
  x <- log10(x + 1)
  x <- (x - mean(x))/sd(x)
  x
})

## Tissues

EVO.Tissues.raw <- data.frame(Phylostratum = PROTAS$Phylostratum, GeneID = PROTAS$GeneID,
                                Adult_Needle = rowMeans(PROTAS[,3:5]), 
                                Juvenile_Needle = rowMeans(PROTAS[,6:8]),
                                Roots = rowMeans(PROTAS[,9:11]),
                                Bud = rowMeans(PROTAS[,12:14])
)

EVO.Tissues.Means.log10 <- EVO.Tissues.raw
EVO.Tissues.Means.log10[,3:ncol(EVO.Tissues.Means.log10)] <- apply(EVO.Tissues.Means.log10[,3:ncol(EVO.Tissues.Means.log10)], MARGIN = 2, FUN = function(x){
  x <- log10(x + 1.1)
  x
})


EVO.Tissues.log10 <- data.frame(Phylostratum = EVOPROTAS.log10$Phylostratum, GeneID = EVOPROTAS.log10$GeneID,
                                Adult_Needle = rowMeans(EVOPROTAS.log10[,3:5]), 
                                Juvenile_Needle = rowMeans(EVOPROTAS.log10[,6:8]),
                                Roots = rowMeans(EVOPROTAS.log10[,9:11]),
                                Bud = rowMeans(EVOPROTAS.log10[,12:14])
                                )

EVO.Tissues.log10.zscore <- data.frame(Phylostratum = EVOPROTAS.log10.zscore$Phylostratum, GeneID = EVOPROTAS.log10.zscore$GeneID,
                                Adult_Needle = rowMeans(EVOPROTAS.log10.zscore[,3:5]), 
                                Juvenile_Needle = rowMeans(EVOPROTAS.log10.zscore[,6:8]),
                                Roots = rowMeans(EVOPROTAS.log10.zscore[,9:11]),
                                Bud = rowMeans(EVOPROTAS.log10.zscore[,12:14])
)


is.ExpressionSet(EVO.Tissues.log10) 
is.ExpressionSet(EVO.Tissues.log10.zscore) 

PlotSignature(EVO.Tissues.raw, measure = "TAI", TestStatistic = "FlatLineTest", xlab = "Tissues", ylab = "PAI", permutations = 10000)
PlotSignature(EVO.Tissues.log10, measure = "TAI", TestStatistic = "FlatLineTest", xlab = "Tissues", ylab = "PAI", permutations = 10000)
PlotSignature(EVO.Tissues.log10.zscore, measure = "TAI", TestStatistic = "FlatLineTest", xlab = "Tissues", ylab = "PAI", permutations = 10000)

PlotSignature(EVO.Tissues.Means.log10, measure = "TAI", TestStatistic = "FlatLineTest", xlab = "Tissues", ylab = "PAI", permutations = 10000)

## Stress Total

EVO.StressTotal.raw <- data.frame(Phylostratum = PROTAS$Phylostratum, GeneID = PROTAS$GeneID,
                                    Control_FU = rowMeans(PROTAS[,15:19]), 
                                    FU = rowMeans(PROTAS[,20:24]),
                                    Control_Heat = rowMeans(PROTAS[,25:28]),
                                    Heat_T1 = rowMeans(PROTAS[,29:32]),
                                    Heat_T3 = rowMeans(PROTAS[,33:36]),
                                    Control_UV = rowMeans(PROTAS[,37:39]),
                                    UV_T1 = rowMeans(PROTAS[,40:42]),
                                    UV_T2 = rowMeans(PROTAS[,43:45]),
                                    UV_T3 = rowMeans(PROTAS[,46:48]),
                                    UV_R = rowMeans(PROTAS[,49:51])
)

EVO.Stress.Means.log10 <- EVO.StressTotal.raw
EVO.Stress.Means.log10[,3:ncol(EVO.Stress.Means.log10)] <- apply(EVO.Stress.Means.log10[,3:ncol(EVO.Stress.Means.log10)], MARGIN = 2, FUN = function(x){
  x <- log10(x + 1.1)
  x
})

EVO.StressTotal.log10 <- data.frame(Phylostratum = EVOPROTAS.log10$Phylostratum, GeneID = EVOPROTAS.log10$GeneID,
                                Control_FU = rowMeans(EVOPROTAS.log10[,15:19]), 
                                FU = rowMeans(EVOPROTAS.log10[,20:24]),
                                Control_Heat = rowMeans(EVOPROTAS.log10[,25:28]),
                                Heat_T1 = rowMeans(EVOPROTAS.log10[,29:32]),
                                Heat_T3 = rowMeans(EVOPROTAS.log10[,33:36]),
                                Control_UV = rowMeans(EVOPROTAS.log10[,37:39]),
                                UV_T1 = rowMeans(EVOPROTAS.log10[,40:42]),
                                UV_T2 = rowMeans(EVOPROTAS.log10[,43:45]),
                                UV_T3 = rowMeans(EVOPROTAS.log10[,46:48]),
                                UV_R = rowMeans(EVOPROTAS.log10[,49:51])
)

EVO.StressTotal.log10.zscore <- data.frame(Phylostratum = EVOPROTAS.log10.zscore$Phylostratum, GeneID = EVOPROTAS.log10.zscore$GeneID,
                                    Control_FU = rowMeans(EVOPROTAS.log10.zscore[,15:19]), 
                                    FU = rowMeans(EVOPROTAS.log10.zscore[,20:24]),
                                    Control_Heat = rowMeans(EVOPROTAS.log10.zscore[,25:28]),
                                    Heat_T1 = rowMeans(EVOPROTAS.log10.zscore[,29:32]),
                                    Heat_T3 = rowMeans(EVOPROTAS.log10.zscore[,33:36]),
                                    Control_UV = rowMeans(EVOPROTAS.log10.zscore[,37:39]),
                                    UV_T1 = rowMeans(EVOPROTAS.log10.zscore[,40:42]),
                                    UV_T2 = rowMeans(EVOPROTAS.log10.zscore[,43:45]),
                                    UV_T3 = rowMeans(EVOPROTAS.log10.zscore[,46:48]),
                                    UV_R = rowMeans(EVOPROTAS.log10.zscore[,49:51])
)


is.ExpressionSet(EVO.StressTotal.raw)
is.ExpressionSet(EVO.StressTotal.log10) 
is.ExpressionSet(EVO.StressTotal.log10.zscore) 

PlotSignature(EVO.StressTotal.raw, measure = "TAI", TestStatistic = "FlatLineTest", xlab = "Tissues", ylab = "PAI",  permutations = 10000)
PlotSignature(EVO.StressTotal.log10, measure = "TAI", TestStatistic = "FlatLineTest", xlab = "Tissues", ylab = "PAI", permutations = 10000)
PlotSignature(EVO.StressTotal.log10.zscore, measure = "TAI", TestStatistic = "FlatLineTest", xlab = "Tissues", ylab = "PAI",  permutations = 10000)

PlotSignature(EVO.Stress.Means.log10, measure = "TAI", TestStatistic = "FlatLineTest", xlab = "Tissues", ylab = "PAI", permutations = 10000)

## UV Heat Total, Nucleus, Chloroplast

EVO.StressComp.raw <- data.frame(Phylostratum = PROTAS$Phylostratum, GeneID = PROTAS$GeneID,
                                  Total_Control_Heat = rowMeans(PROTAS[,25:28]),
                                  Total_Heat_T1 = rowMeans(PROTAS[,29:32]),
                                  Total_Heat_T3 = rowMeans(PROTAS[,33:36]),
                                  Total_Control_UV = rowMeans(PROTAS[,37:39]),
                                  Total_UV_T1 = rowMeans(PROTAS[,40:42]),
                                  Total_UV_T2 = rowMeans(PROTAS[,43:45]),
                                  Total_UV_T3 = rowMeans(PROTAS[,46:48]),
                                  Total_UV_R = rowMeans(PROTAS[,49:51]),
                                  Nucleus_Control_Heat = rowMeans(PROTAS[,52:55]),
                                  Nucleus_Heat_T1 = rowMeans(PROTAS[,56:59]),
                                  Nucleus_Heat_T3 = rowMeans(PROTAS[,60:63]),
                                  Nucleus_Heat_T5 = rowMeans(PROTAS[,64:67]),
                                  Nucleus_Heat_T10 = rowMeans(PROTAS[,68:71]),
                                  Nucleus_Heat_T5R = rowMeans(PROTAS[,72:75]),
                                  Nucleus_Heat_CR = rowMeans(PROTAS[,76:79]),
                                  Nucleus_Control_UV = rowMeans(PROTAS[,80:83]),
                                  Nucleus_UV_T1 = rowMeans(PROTAS[,84:87]),
                                  Nucleus_UV_T3 = rowMeans(PROTAS[,88:91]),
                                  Chloroplast_Control_Heat = rowMeans(PROTAS[,92:99]),
                                  Chloroplast_Heat_T1 = rowMeans(PROTAS[,100:107]),
                                  Chloroplast_Heat_T3 = rowMeans(PROTAS[,108:115]),
                                  Chloroplast_Heat_T5 = rowMeans(PROTAS[,116:123]),
                                  Chloroplast_Control_UV = rowMeans(PROTAS[,124:131]),
                                  Chloroplast_UV_T1 = rowMeans(PROTAS[,132:139]),
                                  Chloroplast_UV_T2 = rowMeans(PROTAS[,140:147]),
                                  Chloroplast_UV_T3 = rowMeans(PROTAS[,148:155])
)

EVO.StressComp.Means.log10 <- EVO.StressComp.raw
EVO.StressComp.Means.log10[,3:ncol(EVO.StressComp.Means.log10)] <- apply(EVO.StressComp.Means.log10[,3:ncol(EVO.StressComp.Means.log10)], MARGIN = 2, FUN = function(x){
  x <- log10(x + 1.1)
  x
})

EVO.StressComp.log10 <- data.frame(Phylostratum = EVOPROTAS.log10$Phylostratum, GeneID = EVOPROTAS.log10$GeneID,
                                    Total_Control_Heat = rowMeans(EVOPROTAS.log10[,25:28]),
                                    Total_Heat_T1 = rowMeans(EVOPROTAS.log10[,29:32]),
                                    Total_Heat_T3 = rowMeans(EVOPROTAS.log10[,33:36]),
                                    Total_Control_UV = rowMeans(EVOPROTAS.log10[,37:39]),
                                    Total_UV_T1 = rowMeans(EVOPROTAS.log10[,40:42]),
                                    Total_UV_T2 = rowMeans(EVOPROTAS.log10[,43:45]),
                                    Total_UV_T3 = rowMeans(EVOPROTAS.log10[,46:48]),
                                    Total_UV_R = rowMeans(EVOPROTAS.log10[,49:51]),
                                    Nucleus_Control_Heat = rowMeans(EVOPROTAS.log10[,52:55]),
                                    Nucleus_Heat_T1 = rowMeans(EVOPROTAS.log10[,56:59]),
                                    Nucleus_Heat_T3 = rowMeans(EVOPROTAS.log10[,60:63]),
                                    Nucleus_Heat_T5 = rowMeans(EVOPROTAS.log10[,64:67]),
                                    Nucleus_Heat_T10 = rowMeans(EVOPROTAS.log10[,68:71]),
                                    Nucleus_Heat_T5R = rowMeans(EVOPROTAS.log10[,72:75]),
                                    Nucleus_Heat_CR = rowMeans(EVOPROTAS.log10[,76:79]),
                                    Nucleus_Control_UV = rowMeans(EVOPROTAS.log10[,80:83]),
                                    Nucleus_UV_T1 = rowMeans(EVOPROTAS.log10[,84:87]),
                                    Nucleus_UV_T3 = rowMeans(EVOPROTAS.log10[,88:91]),
                                    Chloroplast_Control_Heat = rowMeans(EVOPROTAS.log10[,92:99]),
                                    Chloroplast_Heat_T1 = rowMeans(EVOPROTAS.log10[,100:107]),
                                    Chloroplast_Heat_T3 = rowMeans(EVOPROTAS.log10[,108:115]),
                                    Chloroplast_Heat_T5 = rowMeans(EVOPROTAS.log10[,116:123]),
                                    Chloroplast_Control_UV = rowMeans(EVOPROTAS.log10[,124:131]),
                                    Chloroplast_UV_T1 = rowMeans(EVOPROTAS.log10[,132:139]),
                                    Chloroplast_UV_T2 = rowMeans(EVOPROTAS.log10[,140:147]),
                                    Chloroplast_UV_T3 = rowMeans(EVOPROTAS.log10[,148:155])
)
EVO.StressComp.log10.zscore <- data.frame(Phylostratum = EVOPROTAS.log10.zscore$Phylostratum, GeneID = EVOPROTAS.log10.zscore$GeneID,
                                           Total_Control_Heat = rowMeans(EVOPROTAS.log10.zscore[,25:28]),
                                           Total_Heat_T1 = rowMeans(EVOPROTAS.log10.zscore[,29:32]),
                                           Total_Heat_T3 = rowMeans(EVOPROTAS.log10.zscore[,33:36]),
                                           Total_Control_UV = rowMeans(EVOPROTAS.log10.zscore[,37:39]),
                                           Total_UV_T1 = rowMeans(EVOPROTAS.log10.zscore[,40:42]),
                                           Total_UV_T2 = rowMeans(EVOPROTAS.log10.zscore[,43:45]),
                                           Total_UV_T3 = rowMeans(EVOPROTAS.log10.zscore[,46:48]),
                                           Total_UV_R = rowMeans(EVOPROTAS.log10.zscore[,49:51]),
                                           Nucleus_Control_Heat = rowMeans(EVOPROTAS.log10.zscore[,52:55]),
                                           Nucleus_Heat_T1 = rowMeans(EVOPROTAS.log10.zscore[,56:59]),
                                           Nucleus_Heat_T3 = rowMeans(EVOPROTAS.log10.zscore[,60:63]),
                                           Nucleus_Heat_T5 = rowMeans(EVOPROTAS.log10.zscore[,64:67]),
                                           Nucleus_Heat_T10 = rowMeans(EVOPROTAS.log10.zscore[,68:71]),
                                           Nucleus_Heat_T5R = rowMeans(EVOPROTAS.log10.zscore[,72:75]),
                                           Nucleus_Heat_CR = rowMeans(EVOPROTAS.log10.zscore[,76:79]),
                                           Nucleus_Control_UV = rowMeans(EVOPROTAS.log10.zscore[,80:83]),
                                           Nucleus_UV_T1 = rowMeans(EVOPROTAS.log10.zscore[,84:87]),
                                           Nucleus_UV_T3 = rowMeans(EVOPROTAS.log10.zscore[,88:91]),
                                           Chloroplast_Control_Heat = rowMeans(EVOPROTAS.log10.zscore[,92:99]),
                                           Chloroplast_Heat_T1 = rowMeans(EVOPROTAS.log10.zscore[,100:107]),
                                           Chloroplast_Heat_T3 = rowMeans(EVOPROTAS.log10.zscore[,108:115]),
                                           Chloroplast_Heat_T5 = rowMeans(EVOPROTAS.log10.zscore[,116:123]),
                                           Chloroplast_Control_UV = rowMeans(EVOPROTAS.log10.zscore[,124:131]),
                                           Chloroplast_UV_T1 = rowMeans(EVOPROTAS.log10.zscore[,132:139]),
                                           Chloroplast_UV_T2 = rowMeans(EVOPROTAS.log10.zscore[,140:147]),
                                           Chloroplast_UV_T3 = rowMeans(EVOPROTAS.log10.zscore[,148:155])
)


is.ExpressionSet(EVO.StressComp.raw)
is.ExpressionSet(EVO.StressComp.log10) 
is.ExpressionSet(EVO.StressComp.log10.zscore) 

PlotSignature(EVO.StressComp.raw, measure = "TAI", TestStatistic = "FlatLineTest", xlab = "Tissues", ylab = "PAI",  permutations = 10000)
PlotSignature(EVO.StressComp.log10, measure = "TAI", TestStatistic = "FlatLineTest", xlab = "Tissues", ylab = "PAI", permutations = 10000)
PlotSignature(EVO.StressComp.log10.zscore, measure = "TAI", TestStatistic = "FlatLineTest", xlab = "Tissues", ylab = "PAI",  permutations = 10000)

PlotSignature(EVO.StressComp.Means.log10, measure = "TAI", TestStatistic = "FlatLineTest", xlab = "Tissues", ylab = "PAI", permutations = 10000)

## Population E Population T

EVO.StressPops.raw <- data.frame(Phylostratum = PROTAS$Phylostratum, GeneID = PROTAS$GeneID,
                                 Chloroplast_Control_Heat_E = rowMeans(PROTAS[,92:95]),
                                 Chloroplast_Heat_T1_E = rowMeans(PROTAS[,100:103]),
                                 Chloroplast_Heat_T3_E = rowMeans(PROTAS[,108:111]),
                                 Chloroplast_Heat_T5_E = rowMeans(PROTAS[,116:119]),
                                 Chloroplast_Control_Heat_T = rowMeans(PROTAS[,96:99]),
                                 Chloroplast_Heat_T1_T = rowMeans(PROTAS[,104:107]),
                                 Chloroplast_Heat_T3_T = rowMeans(PROTAS[,112:115]),
                                 Chloroplast_Heat_T5_T = rowMeans(PROTAS[,120:123]),
                                 Chloroplast_Control_UV_E = rowMeans(PROTAS[,124:127]),
                                 Chloroplast_UV_T1_E = rowMeans(PROTAS[,132:135]),
                                 Chloroplast_UV_T2_E = rowMeans(PROTAS[,140:143]),
                                 Chloroplast_UV_T3_E = rowMeans(PROTAS[,148:151]),
                                 Chloroplast_Control_UV_T = rowMeans(PROTAS[,128:131]),
                                 Chloroplast_UV_T1_T = rowMeans(PROTAS[,136:139]),
                                 Chloroplast_UV_T2_T = rowMeans(PROTAS[,144:147]),
                                 Chloroplast_UV_T3_T = rowMeans(PROTAS[,152:155])
)

EVO.StressPops.Means.log10 <- EVO.StressPops.raw
EVO.StressPops.Means.log10[,3:ncol(EVO.StressPops.Means.log10)] <- apply(EVO.StressPops.Means.log10[,3:ncol(EVO.StressPops.Means.log10)], MARGIN = 2, FUN = function(x){
  x <- log10(x + 1.1)
  x
})

EVO.StressPops.log10 <- data.frame(Phylostratum = EVOPROTAS.log10$Phylostratum, GeneID = EVOPROTAS.log10$GeneID,
                                   Chloroplast_Control_Heat_E = rowMeans(EVOPROTAS.log10[,92:95]),
                                   Chloroplast_Heat_T1_E = rowMeans(EVOPROTAS.log10[,100:103]),
                                   Chloroplast_Heat_T3_E = rowMeans(EVOPROTAS.log10[,108:111]),
                                   Chloroplast_Heat_T5_E = rowMeans(EVOPROTAS.log10[,116:119]),
                                   Chloroplast_Control_Heat_T = rowMeans(EVOPROTAS.log10[,96:99]),
                                   Chloroplast_Heat_T1_T = rowMeans(EVOPROTAS.log10[,104:107]),
                                   Chloroplast_Heat_T3_T = rowMeans(EVOPROTAS.log10[,112:115]),
                                   Chloroplast_Heat_T5_T = rowMeans(EVOPROTAS.log10[,120:123]),
                                   Chloroplast_Control_UV_E = rowMeans(EVOPROTAS.log10[,124:127]),
                                   Chloroplast_UV_T1_E = rowMeans(EVOPROTAS.log10[,132:135]),
                                   Chloroplast_UV_T2_E = rowMeans(EVOPROTAS.log10[,140:143]),
                                   Chloroplast_UV_T3_E = rowMeans(EVOPROTAS.log10[,148:151]),
                                   Chloroplast_Control_UV_T = rowMeans(EVOPROTAS.log10[,128:131]),
                                   Chloroplast_UV_T1_T = rowMeans(EVOPROTAS.log10[,136:139]),
                                   Chloroplast_UV_T2_T = rowMeans(EVOPROTAS.log10[,144:147]),
                                   Chloroplast_UV_T3_T = rowMeans(EVOPROTAS.log10[,152:155])
)
EVO.StressPops.log10.zscore <- data.frame(Phylostratum = EVOPROTAS.log10.zscore$Phylostratum, GeneID = EVOPROTAS.log10.zscore$GeneID,
                                           Chloroplast_Control_Heat_E = rowMeans(EVOPROTAS.log10.zscore[,92:95]),
                                           Chloroplast_Heat_T1_E = rowMeans(EVOPROTAS.log10.zscore[,100:103]),
                                           Chloroplast_Heat_T3_E = rowMeans(EVOPROTAS.log10.zscore[,108:111]),
                                           Chloroplast_Heat_T5_E = rowMeans(EVOPROTAS.log10.zscore[,116:119]),
                                           Chloroplast_Control_Heat_T = rowMeans(EVOPROTAS.log10.zscore[,96:99]),
                                           Chloroplast_Heat_T1_T = rowMeans(EVOPROTAS.log10.zscore[,104:107]),
                                           Chloroplast_Heat_T3_T = rowMeans(EVOPROTAS.log10.zscore[,112:115]),
                                           Chloroplast_Heat_T5_T = rowMeans(EVOPROTAS.log10.zscore[,120:123]),
                                           Chloroplast_Control_UV_E = rowMeans(EVOPROTAS.log10.zscore[,124:127]),
                                           Chloroplast_UV_T1_E = rowMeans(EVOPROTAS.log10.zscore[,132:135]),
                                           Chloroplast_UV_T2_E = rowMeans(EVOPROTAS.log10.zscore[,140:143]),
                                           Chloroplast_UV_T3_E = rowMeans(EVOPROTAS.log10.zscore[,148:151]),
                                           Chloroplast_Control_UV_T = rowMeans(EVOPROTAS.log10.zscore[,128:131]),
                                           Chloroplast_UV_T1_T = rowMeans(EVOPROTAS.log10.zscore[,136:139]),
                                           Chloroplast_UV_T2_T = rowMeans(EVOPROTAS.log10.zscore[,144:147]),
                                           Chloroplast_UV_T3_T = rowMeans(EVOPROTAS.log10.zscore[,152:155])
)


is.ExpressionSet(EVO.StressPops.raw)
is.ExpressionSet(EVO.StressPops.log10) 
is.ExpressionSet(EVO.StressPops.log10.zscore) 

PlotSignature(EVO.StressPops.raw, measure = "TAI", TestStatistic = "FlatLineTest", xlab = "Tissues", ylab = "PAI",  permutations = 10000)
PlotSignature(EVO.StressPops.log10, measure = "TAI", TestStatistic = "FlatLineTest", xlab = "Tissues", ylab = "PAI", permutations = 1000)
PlotSignature(EVO.StressPops.log10.zscore, measure = "TAI", TestStatistic = "FlatLineTest", xlab = "Tissues", ylab = "PAI",  permutations = 10000)

PlotSignature(EVO.StressPops.Means.log10, measure = "TAI", TestStatistic = "FlatLineTest", xlab = "Tissues", ylab = "PAI", permutations = 10000)

#### 3. Integration  ####

reticulate::use_python("C:\\Users\\bboyl\\AppData\\Local\\Programs\\Python\\Python38\\")
setwd("F:/ATLAS_Pra/Protein_module/0.Process/Figures/MOFA/")

#### 3.1 MOFA multigroup Stresses Total | FU; UV; Heat ####

ProtAtlas <- ProtAtlas.jic

### 3.0 Per Location UV HS alone

### 3.0 Per population

## 3.1.1 Format data and import to MOFA

MOFAStressTotal <- log1p(rawinputs$raw)
rownames(MOFAStressTotal) <- Samples
MOFAStressTotal <- t(MOFAStressTotal)
MOFAStressTotal <- MOFAStressTotal[,c(22:ncol(MOFAStressTotal))]
ProteinsList <- list(Proteins = MOFAStressTotal)
# careful HVF by each group and then take the union to take into account the group effect before selecting HVFs (MOFA log recommendation)
ProteinsList.vars.FU <- apply(ProteinsList$Proteins[,1:10], 1, var)
FU <- names(which(ProteinsList.vars.FU > 0))
FU <- names(sort(apply(ProteinsList$Proteins[,1:10], 1, var),decreasing = T)[1:500])
ProteinsList.vars.HS <- apply(ProteinsList$Proteins[,c(11:22, 38:65, 78:109)], 1, var)
HS <- names(which(ProteinsList.vars.HS > 0))
HS <- names(sort(apply(ProteinsList$Proteins[,c(11:22, 38:65, 78:109)], 1, var),decreasing = T)[1:500])
ProteinsList.vars.UV <- apply(ProteinsList$Proteins[,c(23:37, 66:77, 110:141)], 1, var)
UV <- names(which(ProteinsList.vars.UV > 0))
UV <- names(sort(apply(ProteinsList$Proteins[,c(23:37, 66:77, 110:141)], 1, var),decreasing = T)[1:500])
features <- unique(c(FU,HS,UV))
ProteinsList$Proteins <- ProteinsList$Proteins[rownames(ProteinsList$Proteins) %in% features,] 
ProteinsList$Proteins <- as.matrix(ProteinsList$Proteins)

groups <- c(rep("FU",10), rep("HS",12), rep("UV",15), rep("HS", length(38:65)), rep("UV", length(66:77)), 
            rep("HS", length(78:109)), rep("UV", length(110:141)))

MOFAgrouped <- create_mofa(data = ProteinsList, groups = groups, save_metadata = TRUE)

## 3.1.2 Prepare MOFA and run

plot_data_overview(MOFAgrouped)
sample_metadata <- data.frame(
  sample = colnames(ProteinsList$Proteins),
  stress = c(rep("C",5), rep("S",5), rep("C",4), rep("S",8), rep("C",3), rep("S",9), rep("R",3)),
  treatment =c(rep("C",5), rep("FU",5), rep("C",4), rep("HS",8), rep("C",3), rep("UV",9), rep("UVR",3)),
  intensity = c(labelitass$Intensity[13:49]),
  group = groups,
  sample_number = c(1:ncol(ProteinsList$Proteins))
)
sample_metadata <- data.frame(
  sample = colnames(ProteinsList$Proteins),
  stress = c(labelitass$Stress[13:153]),
  treatment =c(labelitass$Stress[13:153]),
  intensity = c(labelitass$Intensity[13:153]),
  location = c(labelitass$Technique[13:153]),
  population = c(labelitass$Intensity[13:153]),
  group = groups,
  sample_number = c(1:ncol(ProteinsList$Proteins))
)
sample_metadata$intensity <- gsub(pattern = "T3", replacement = "T2", sample_metadata$intensity, fixed = T)
sample_metadata$intensity <- gsub(pattern = "T5", replacement = "T3", sample_metadata$intensity, fixed = T)
sample_metadata$intensity <- gsub(pattern = "T10", replacement = "T4", sample_metadata$intensity, fixed = T)
sample_metadata$stress <- gsub(pattern = "FU", replacement = "Stress", sample_metadata$stress, fixed = T)
sample_metadata$stress <- gsub(pattern = "Heat", replacement = "Stress", sample_metadata$stress, fixed = T)
sample_metadata$stress <- gsub(pattern = "UV", replacement = "Stress", sample_metadata$stress, fixed = T)

data_opts <- get_default_data_options(MOFAgrouped)
model_opts <- get_default_model_options(MOFAgrouped)
model_opts$num_factors <- 10
train_opts <- get_default_training_options(MOFAgrouped)
train_opts$maxiter <- "1000000"
train_opts$convergence_mode <- "slow"
train_opts$drop_factor_threshold <- 0.05
train_opts$seed <- 404
outfileGrouped <- "F:/ATLAS_Pra/Protein_module/0.Process/Figures/MOFA/StressTotal/MOFAStressTotal_log101_V0.hdf5"

MOFAgrouped <- prepare_mofa(MOFAgrouped, data_options = data_opts, model_options = model_opts, training_options = train_opts)
MOFAgrouped.trained <- run_mofa(MOFAgrouped, outfile = outfileGrouped, use_basilisk = F)

## 3.1.3 Plot General

plot_factor_cor(MOFAgrouped.trained)
plot_variance_explained(MOFAgrouped.trained, x = "group", y = "factor")
plot_variance_explained(MOFAgrouped.trained, y = "factor")
plot_variance_explained(MOFAgrouped.trained, y = "factor", plot_total = T)

MOFAgrouped.trained@cache$variance_explained
MOFAgrouped.trained@cache$variance_explained$r2_per_factor$female
MOFAgrouped.trained@cache$variance_explained$r2_per_factor$male

a <- plot_variance_explained(MOFAgrouped.trained, x="view", y="factor")
b <- plot_variance_explained(MOFAgrouped.trained, x="group", y="factor", plot_total = T)[[2]]
gridExtra::grid.arrange(b, a, ncol = 1, nrow = 2)

## 3.1.4 LFs Inspection

# Single

MOFAgrouped.trained@samples_metadata <- sample_metadata

plot_factor(MOFAgrouped.trained, 
            factor = 1,
            color_by = "location",
            shape_by = "stress", dot_size = 3
)

# Multiple

plot_factors(MOFAgrouped.trained, 
             factors = c(1,2,3,4),
             color_by = "group", dot_size = 8, shape_by = "group"
)

plot_factors(MOFAgrouped.trained, 
             factors = c(1,2,3,4),
             color_by = "stress", dot_size = 3, shape_by = "location"
)

# 3.1.1 Format data and import to MOFA

MOFAStressTotal <- log1p(rawinputs$raw)
rownames(MOFAStressTotal) <- Samples
MOFAStressTotal <- t(MOFAStressTotal)
MOFAStressTotal <- MOFAStressTotal[,c(22:ncol(MOFAStressTotal))]
ProteinsList <- list(Proteins = MOFAStressTotal)
# careful HVF by each group and then take the union to take into account the group effect before selecting HVFs (MOFA log recommendation)
ProteinsList.vars.FU <- apply(ProteinsList$Proteins[,1:10], 1, var)
FU <- names(which(ProteinsList.vars.FU > 0))
FU <- names(sort(apply(ProteinsList$Proteins[,1:10], 1, var),decreasing = T)[1:500])
ProteinsList.vars.HS <- apply(ProteinsList$Proteins[,c(11:22, 38:65, 78:109)], 1, var)
HS <- names(which(ProteinsList.vars.HS > 0))
HS <- names(sort(apply(ProteinsList$Proteins[,c(11:22, 38:65, 78:109)], 1, var),decreasing = T)[1:500])
ProteinsList.vars.UV <- apply(ProteinsList$Proteins[,c(23:37, 66:77, 110:141)], 1, var)
UV <- names(which(ProteinsList.vars.UV > 0))
UV <- names(sort(apply(ProteinsList$Proteins[,c(23:37, 66:77, 110:141)], 1, var),decreasing = T)[1:500])
features <- unique(c(FU,HS,UV))
ProteinsList$Proteins <- ProteinsList$Proteins[rownames(ProteinsList$Proteins) %in% features,] 
ProteinsList$Proteins <- as.matrix(ProteinsList$Proteins)

groups <- c(rep("FU",10), rep("HS",12), rep("UV",15), rep("HS", length(38:65)), rep("UV", length(66:77)), 
            rep("HS", length(78:109)), rep("UV", length(110:141)))

MOFAgrouped <- create_mofa(data = ProteinsList, groups = groups, save_metadata = TRUE)

## 3.1.2 Prepare MOFA and run

plot_data_overview(MOFAgrouped)
sample_metadata <- data.frame(
  sample = colnames(ProteinsList$Proteins),
  stress = c(rep("C",5), rep("S",5), rep("C",4), rep("S",8), rep("C",3), rep("S",9), rep("R",3)),
  treatment =c(rep("C",5), rep("FU",5), rep("C",4), rep("HS",8), rep("C",3), rep("UV",9), rep("UVR",3)),
  intensity = c(labelitass$Intensity[13:49]),
  group = groups,
  sample_number = c(1:ncol(ProteinsList$Proteins))
)
sample_metadata <- data.frame(
  sample = colnames(ProteinsList$Proteins),
  stress = c(labelitass$Stress[13:153]),
  treatment =c(labelitass$Stress[13:153]),
  intensity = c(labelitass$Intensity[13:153]),
  location = c(labelitass$Technique[13:153]),
  population = c(labelitass$Intensity[13:153]),
  group = groups,
  sample_number = c(1:ncol(ProteinsList$Proteins))
)
sample_metadata$intensity <- gsub(pattern = "T3", replacement = "T2", sample_metadata$intensity, fixed = T)
sample_metadata$intensity <- gsub(pattern = "T5", replacement = "T3", sample_metadata$intensity, fixed = T)
sample_metadata$intensity <- gsub(pattern = "T10", replacement = "T4", sample_metadata$intensity, fixed = T)
sample_metadata$stress <- gsub(pattern = "FU", replacement = "Stress", sample_metadata$stress, fixed = T)
sample_metadata$stress <- gsub(pattern = "Heat", replacement = "Stress", sample_metadata$stress, fixed = T)
sample_metadata$stress <- gsub(pattern = "UV", replacement = "Stress", sample_metadata$stress, fixed = T)

data_opts <- get_default_data_options(MOFAgrouped)
model_opts <- get_default_model_options(MOFAgrouped)
model_opts$num_factors <- 10
train_opts <- get_default_training_options(MOFAgrouped)
train_opts$maxiter <- "1000000"
train_opts$convergence_mode <- "slow"
train_opts$drop_factor_threshold <- 0.05
train_opts$seed <- 404
outfileGrouped <- "F:/ATLAS_Pra/Protein_module/0.Process/Figures/MOFA/StressTotal/MOFAStressTotal_log101_V0.hdf5"

MOFAgrouped <- prepare_mofa(MOFAgrouped, data_options = data_opts, model_options = model_opts, training_options = train_opts)
MOFAgrouped.trained <- run_mofa(MOFAgrouped, outfile = outfileGrouped, use_basilisk = F)

## 3.1.3 Plot General

plot_factor_cor(MOFAgrouped.trained)
plot_variance_explained(MOFAgrouped.trained, x = "group", y = "factor")
plot_variance_explained(MOFAgrouped.trained, y = "factor")
plot_variance_explained(MOFAgrouped.trained, y = "factor", plot_total = T)

MOFAgrouped.trained@cache$variance_explained
MOFAgrouped.trained@cache$variance_explained$r2_per_factor$female
MOFAgrouped.trained@cache$variance_explained$r2_per_factor$male

a <- plot_variance_explained(MOFAgrouped.trained, x="view", y="factor")
b <- plot_variance_explained(MOFAgrouped.trained, x="group", y="factor", plot_total = T)[[2]]
gridExtra::grid.arrange(b, a, ncol = 1, nrow = 2)

## 3.1.4 LFs Inspection

# Single

MOFAgrouped.trained@samples_metadata <- sample_metadata

plot_factor(MOFAgrouped.trained, 
            factor = 1,
            color_by = "location",
            shape_by = "stress", dot_size = 3
)

# Multiple

plot_factors(MOFAgrouped.trained, 
             factors = c(1,2,3,4),
             color_by = "group", dot_size = 8, shape_by = "group"
)

plot_factors(MOFAgrouped.trained, 
             factors = c(1,2,3,4),
             color_by = "stress", dot_size = 3, shape_by = "location"
)



#### 3.2 MOFA multigroup TotalStresses ####

## 3.2.1 Format data and import to MOFA

ProtAtlas.Stress <- ProtAtlas[,c(1,15:51)]
AllProteinsList <- list(Proteins = ProtAtlas.Stress)
AllProteinsList$Proteins[,2:ncol(AllProteinsList$Proteins)] <- log10(AllProteinsList$Proteins[,2:ncol(AllProteinsList$Proteins)] + 1)
rownames(AllProteinsList$Proteins) <- AllProteinsList$Proteins$ProtID
AllProteinsList$Proteins <- AllProteinsList$Proteins[,-1]

# limma
pheno <- data.frame(sample = 1:length(15:51), 
                    outcome = sample_metadata$stress,
                    batch = c(rep(1,10), rep(2,12), rep(3,15)))
pheno$outcome <- gsub(pattern = "R", replacement = "C", x = pheno$outcome)
rownames(pheno) <- colnames(AllProteinsList$Proteins)
mm = model.matrix(~as.factor(outcome), data=pheno)
mat <- limma::removeBatchEffect(AllProteinsList$Proteins, batch=pheno$batch, design=mm)
write.table(x = mat ,file = "F:/ATLAS_Pra/Protein_module/0.Process/Figures/MOFA/StressTotal/StressTotal_false_v02/StressTotal_LimmaOutcomeStress.txt", sep = "\t", quote = F, col.names = T, row.names = T)


AllProteinsList.vars <- apply(AllProteinsList$Proteins, 1, var)
AllProteinsList$Proteins <- AllProteinsList$Proteins[names(which(AllProteinsList.vars > 0)),] 
AllProteinsList$Proteins <- as.matrix(AllProteinsList$Proteins)

groups <- c(rep("FU",10), rep("HS",12), rep("UV",15))

MOFAgrouped <- create_mofa(data = AllProteinsList, groups = groups, save_metadata = TRUE)

## 3.2.2 Prepare MOFA and run

plot_data_overview(MOFAgrouped)
sample_metadata <- data.frame(
  sample = colnames(AllProteinsList$Proteins),
  stress = c(rep("C",5), rep("S",5), rep("C",4), rep("S",8), rep("C",3), rep("S",9), rep("R",3)),
  treatment =c(rep("C",5), rep("FU",5), rep("C",4), rep("HS-T1",4), rep("HS-T3",4), rep("C",3), rep("UV-T1",3), rep("UV-T2",3), rep("UV-T3",3), rep("UVR",3)),
  time = c(rep("C",5), rep("FU",5), rep("C",4), rep("T1",4), rep("T3",4), rep("C",3), rep("T1",3), rep("T2",3), rep("T3",3), rep("UVR",3)),
  group = groups,
  sample_number = c(1:37)
)
data_opts <- get_default_data_options(MOFAgrouped)
model_opts <- get_default_model_options(MOFAgrouped)
model_opts$num_factors <- 5
train_opts <- get_default_training_options(MOFAgrouped)
train_opts$maxiter <- "100000"
train_opts$convergence_mode <- "slow"
train_opts$drop_factor_threshold <- 0.005
train_opts$seed <- 400
data_opts$scale_views <- FALSE
outfileGrouped <- "F:/ATLAS_Pra/Protein_module/0.Process/MOFA2grouped_log101_V0_rescue.hdf5"

MOFAgrouped <- prepare_mofa(MOFAgrouped, data_options = data_opts, model_options = model_opts, training_options = train_opts)
MOFAgrouped.trained <- run_mofa(MOFAgrouped, outfile = outfileGrouped)

MOFAgrouped <- prepare_mofa(MOFAgrouped, data_options = data_opts, model_options = model_opts, training_options = train_opts)
MOFAgrouped.trained <- run_mofa(MOFAgrouped, outfile = outfileGrouped)

MOFAgrouped <- load_model("F:/ATLAS_Pra/Protein_module/0.Process/Figures/MOFA/StressTotal/StressTotal_false_v02/MOFA2grouped_StressTotal_LimmaOutcomeStress_log101_V0_400_false.hdf5")
MOFAgrouped.trained <- MOFAgrouped

## 3.2.3 Plot General

plot_factor_cor(MOFAgrouped.trained)
plot_variance_explained(MOFAgrouped.trained, x = "group", y = "factor")
plot_variance_explained(MOFAgrouped.trained, x = "group", y = "factor", plot_total = T)

MOFAgrouped.trained@cache$variance_explained
MOFAgrouped.trained@cache$variance_explained$r2_per_factor$FU
MOFAgrouped.trained@cache$variance_explained$r2_per_factor$HS
MOFAgrouped.trained@cache$variance_explained$r2_per_factor$UV

a <- plot_variance_explained(MOFAgrouped.trained, x="view", y="factor", factors = 1:4)
b <- plot_variance_explained(MOFAgrouped.trained, x="group", y="factor", plot_total = T, factors = 1:4)[[2]]
PlotVariance <- grid.arrange(b, a, ncol = 1, nrow = 2)

## 3.2.4 LFs Inspection

# Single

MOFAgrouped.trained@samples_metadata <- sample_metadata

plot_factor(MOFAgrouped.trained, 
            factor = c(1,2,3),
            color_by = "time",
            shape_by = "stress", dot_size = 4
)

# Multiple

plot_factors(MOFAgrouped.trained, 
             factors = c(1,2,3,4),
             color_by = "stress", dot_size = 8, shape_by = "treatment"
)

plot_factors(MOFAgrouped.trained, 
             factors = c(1,2,3),
             color_by = "stress", dot_size = 3, shape_by = "sex"
)

# Weights of biologically relevant LFs

plot_weights(MOFAgrouped.trained,
             view = "Proteins",
             factor = 1,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

plot_weights(MOFAgrouped.trained,
             view = "Proteins",
             factor = 2,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

plot_weights(MOFAgrouped.trained,
             view = "Proteins",
             factor = 3,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

plot_weights(MOFAgrouped.trained,
             view = "Proteins",
             factor = 4,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

LF1 <- plot_top_weights(MOFAgrouped.trained,
                 view = "Proteins",
                 factor = 1,
                 nfeatures = 10,     # Number of features to highlight
                 scale = T,          # Scale weights from -1 to 1
                 abs = F  )           # Take the absolute value?)

LF2 <- plot_top_weights(MOFAgrouped.trained,
                 view = "Proteins",
                 factor = 2,
                 nfeatures = 10,     # Number of features to highlight
                 scale = T,          # Scale weights from -1 to 1
                 abs = F  )           # Take the absolute value?)

LF3 <- plot_top_weights(MOFAgrouped.trained,
                 view = "Proteins",
                 factor = 3,
                 nfeatures = 10,     # Number of features to highlight
                 scale = T,          # Scale weights from -1 to 1
                 abs = F  )           # Take the absolute value?)

LF4 <- plot_top_weights(MOFAgrouped.trained,
                 view = "Proteins",
                 factor = 4,
                 nfeatures = 10,     # Number of features to highlight
                 scale = T,          # Scale weights from -1 to 1
                 abs = F  )           # Take the absolute value?)


gridExtra::grid.arrange(LF1, LF2,LF3 ,LF4, ncol = 2, nrow = 2)



PROTAS.annotation$Symbol[grep("all-04-326945", PROTAS.annotation$ProtID)]

# Enrichment Analyses #
# Positive and neg splitted
# Mercator bins, PS age, gene family founder PS

# Mercator

Bin_db

# Gene Age

PROTIS <- read.delim(file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/Protein_module/ProtAtlas.Analyses.txt", header = T, sep = "\t")
Ages <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/4.GenEra/Output/Pinus_radiata_proteome/3347_gene_ages.tsv", header = T, sep = "\t")
Ages$X.gene <- do.call(rbind, strsplit(x = Ages$X.gene, split = ".", fixed = T))[,1] 
PROTIS$Phylostratum <- NA
PROTIS <- PROTIS[PROTIS$ProtID %in% Ages$X.gene,]

for(i in 1:nrow(PROTIS)){
  age <- unique(Ages$rank[which(Ages$X.gene == PROTIS$ProtID[i])])
  if(length(age) == 1){
    PROTIS$Phylostratum[i] <- age
  }else if(length(age) > 1){
    PROTIS$Phylostratum[i] <- min(age)
  }else if(length(age) == 0){
    stop("Something went wrong boy ... \n")
  }
}

PROTIS <- PROTIS[,c(ncol(PROTIS),1:ncol(PROTIS)-1)]
colnames(PROTIS)[2] <- "GeneID"
PROTIS <- PROTIS[,-3]
PROTIS <- PROTIS[!is.na(PROTIS$Phylostratum),]

PROTIS$Phylostratum <- gsub("15", "12", PROTIS$Phylostratum, fixed = T) 
PROTIS$Phylostratum <- gsub("13", "12", PROTIS$Phylostratum, fixed = T) 
PROTIS$Phylostratum <- gsub("14", "12", PROTIS$Phylostratum, fixed = T) 

PS_db <- list()
PS_names <- levels(as.factor(PROTIS$Phylostratum))

for(i in PS_names){
  PS_db[[i]] <- PROTIS$GeneID[PROTIS$Phylostratum == i]
}

names(PS_db) <- paste0("PS", names(PS_db))
PS_db

# Gene family founder

PROTIS <- read.delim(file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/Protein_module/ProtAtlas.Analyses.txt", header = T, sep = "\t")
Ages <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/4.GenEra/Output/Pinus_radiata_proteome/3347_founder_events.tsv", header = T, sep = "\t")
Ages$rank <- gsub("15", "12", Ages$rank, fixed = T) 
Ages$rank <- gsub("13", "12", Ages$rank, fixed = T) 
Ages$rank <- gsub("14", "12", Ages$rank, fixed = T) 
Ages$FamilyID <- paste0("Founder",rownames(Ages),".", Ages$rank, ".",Ages$family_size)
Ages <- tidyr::separate_longer_delim(Ages, X.gene_family, delim = ",")
Ages$X.gene <- do.call(rbind, strsplit(x = Ages$X.gene, split = ".", fixed = T))[,1] 
PROTIS$Phylostratum <- NA
PROTIS <- PROTIS[PROTIS$ProtID %in% Ages$X.gene,]

for(i in 1:nrow(PROTIS)){
  age <- unique(Ages$rank[which(Ages$X.gene == PROTIS$ProtID[i])])
  if(length(age) == 1){
    PROTIS$Phylostratum[i] <- age
  }else if(length(age) > 1){
    PROTIS$Phylostratum[i] <- min(age)
  }else if(length(age) == 0){
    stop("Something went wrong boy ... \n")
  }
}

PROTIS <- PROTIS[,c(ncol(PROTIS),1:ncol(PROTIS)-1)]
colnames(PROTIS)[2] <- "GeneID"
PROTIS <- PROTIS[,-3]
PROTIS <- PROTIS[!is.na(PROTIS$Phylostratum),]

PSF_db <- list()
PSF_names <- levels(as.factor(PROTIS$Phylostratum))

for(i in PSF_names){
  PSF_db[[i]] <- PROTIS$GeneID[PROTIS$Phylostratum == i]
}

names(PSF_db) <- paste0("PS", names(PSF_db))
PSF_db

# Compute enrichments; LF1-4; Pos and neg; Mercator, Gene Age and Family foundation

Bin.paths <- list_to_matrix(Bin_db)
PS.paths <- list_to_matrix(PS_db)
PSF.paths <- list_to_matrix(PSF_db)


enrichment.parametric.Bin.all <- run_enrichment(MOFAgrouped.trained,
                                               view = "Proteins", factors = c(1,2,3,4),
                                               feature.sets = t(Bin.paths),
                                               sign = "all", set.statistic = "rank.sum",
                                               # statistical.test = "permutation",
                                               min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                               #nperm = 5000
)


enrichment.parametric.PS.all <- run_enrichment(MOFAgrouped.trained,
                                                view = "Proteins", factors = c(1,2,3,4),
                                                feature.sets = t(PS.paths),
                                                sign = "all", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, alpha = 0.1
                                                #nperm = 5000
)

enrichment.parametric.PSF.all <- run_enrichment(MOFAgrouped.trained,
                                                view = "Proteins", factors = c(1,2,3,4),
                                                feature.sets = t(PSF.paths),
                                                sign = "all", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, alpha = 0.1
                                                #nperm = 5000
)

enrichment.parametric.Bin.pos <- run_enrichment(MOFAgrouped.trained,
                                                view = "Proteins", factors = c(1,2,3,4),
                                                feature.sets = t(Bin.paths),
                                                sign = "pos", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                #nperm = 5000
)


enrichment.parametric.PS.pos <- run_enrichment(MOFAgrouped.trained,
                                               view = "Proteins", factors = c(1,2,3,4),
                                               feature.sets = t(PS.paths),
                                               sign = "pos", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                               # statistical.test = "permutation",
                                               min.size = 5, alpha = 0.1
                                               #nperm = 5000
)

enrichment.parametric.PSF.pos <- run_enrichment(MOFAgrouped.trained,
                                                view = "Proteins", factors = c(1,2,3,4),
                                                feature.sets = t(PSF.paths),
                                                sign = "pos", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, alpha = 0.1
                                                #nperm = 5000
)

enrichment.parametric.Bin.neg <- run_enrichment(MOFAgrouped.trained,
                                                view = "Proteins", factors = c(1,2,3,4),
                                                feature.sets = t(Bin.paths),
                                                sign = "neg", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                #nperm = 5000
)


enrichment.parametric.PS.neg <- run_enrichment(MOFAgrouped.trained,
                                               view = "Proteins", factors = c(1,2,3,4),
                                               feature.sets = t(PS.paths),
                                               sign = "neg", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                               # statistical.test = "permutation",
                                               min.size = 5, alpha = 0.1
                                               #nperm = 5000
)

enrichment.parametric.PSF.neg <- run_enrichment(MOFAgrouped.trained,
                                                view = "Proteins", factors = c(1,2,3,4),
                                                feature.sets = t(PSF.paths),
                                                sign = "neg", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, alpha = 0.1
                                                #nperm = 5000
)

## plot enrich

order <- sort(rownames(enrichment.parametric.Bin.all$pval.adj))

Bin.df <- t(enrichment.parametric.Bin.all$pval.adj)[,order]

for(i in 1:nrow(Bin.df)){
  for(j in 1:ncol(Bin.df)){
    if(Bin.df[i,j] > 0.1 ){
      Bin.df[i,j] <- 1
    }
  }
}

Bin.df <- -log10(Bin.df)

for(i in 1:nrow(Bin.df)){
  for(j in 1:ncol(Bin.df)){
    if(Bin.df[i,j] > 10 ){
      Bin.df[i,j] <- 10
    }
  }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Greens"))(500)

ht.Bin <- Heatmap(t(Bin.df), name = "Bins.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
                       )

order <- sort(rownames(enrichment.parametric.PS.all$pval.adj))

PS.df <- t(enrichment.parametric.PS.all$pval.adj)[,order]

for(i in 1:nrow(PS.df)){
  for(j in 1:ncol(PS.df)){
    if(PS.df[i,j] > 0.1 ){
      PS.df[i,j] <- 1
    }
  }
}

PS.df <- -log10(PS.df)

for(i in 1:nrow(PS.df)){
  for(j in 1:ncol(PS.df)){
    if(PS.df[i,j] > 10 ){
      PS.df[i,j] <- 10
    }
  }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Purples"))(100)

ht.PS <- Heatmap(t(PS.df), name = "PS.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

order <- sort(rownames(enrichment.parametric.PSF.all$pval.adj))

PSF.df <- t(enrichment.parametric.PSF.all$pval.adj)[,order]

for(i in 1:nrow(PSF.df)){
  for(j in 1:ncol(PSF.df)){
    if(PSF.df[i,j] > 0.1 ){
      PSF.df[i,j] <- 1
    }
  }
}

PSF.df <- -log10(PSF.df)

for(i in 1:nrow(PSF.df)){
  for(j in 1:ncol(PSF.df)){
    if(PSF.df[i,j] > 10 ){
      PSF.df[i,j] <- 10
    }
  }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100)

ht.PSF <- Heatmap(t(PSF.df), name = "PSF.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

draw(ht.Bin %v% ht.PS %v% ht.PSF, ht_gap = unit(1, "cm"))

# nVenn or Upset same gene family: not that useful (all Stress Total Proteins comes from same 803 founder events)
factors <- get_factors(object =MOFAgrouped.trained, as.data.frame = T)
weights <- get_weights(object = MOFAgrouped.trained, as.data.frame = T)

VennFamily <- Ages[,c(6,5)]
LF1.50 <- names(sort(abs(deframe(weights[which(weights$factor == "Factor1"), c(1,3)])), decreasing = T)[1:500])
LF1.50 <-  VennFamily$FamilyID[match(x =  LF1.50,table = VennFamily$X.gene)]
LF2.50 <- names(sort(abs(deframe(weights[which(weights$factor == "Factor2"), c(1,3)])), decreasing = T)[1:500])
LF2.50 <-  VennFamily$FamilyID[match(x =  LF2.50,table = VennFamily$X.gene)]
LF3.50 <- names(sort(abs(deframe(weights[which(weights$factor == "Factor3"), c(1,3)])), decreasing = T)[1:500])
LF3.50 <-  VennFamily$FamilyID[match(x =  LF3.50,table = VennFamily$X.gene)]
LF4.50 <- names(sort(abs(deframe(weights[which(weights$factor == "Factor4"), c(1,3)])), decreasing = T)[1:500])
LF4.50 <-  VennFamily$FamilyID[match(x =  LF4.50,table = VennFamily$X.gene)]
dataList <- list(LF1 = LF1.50, LF2 = LF2.50, LF3 = LF3.50, LF4 = LF4.50)
#myV <- nVennR::plotVenn(data, nCycles = 20000, opacity = 0.2, borderWidth = 3, systemShow = T, fontScale = 2, outFile = paste0(output.path, "nVennR.svg"))
upset(fromList(dataList), sets = names(dataList), order.by = "freq", 
      keep.order = TRUE, empty.intersections = F, mb.ratio = c(0.6,0.4), sets.bar.color = "navajowhite3")

#### 3.2 MOFA multigroup Stress Compartiments (only UV and HS) ####

## 3.2.1 Format data and import to MOFA

ProtAtlas.Comp <- ProtAtlas[,c(1,25:155)]
AllProteinsList <- list(Proteins = ProtAtlas.Comp)
AllProteinsList$Proteins[,2:ncol(AllProteinsList$Proteins)] <- log10(AllProteinsList$Proteins[,2:ncol(AllProteinsList$Proteins)] + 1)
rownames(AllProteinsList$Proteins) <- AllProteinsList$Proteins$ProtID
AllProteinsList$Proteins <- AllProteinsList$Proteins[,-1]
pheno <- data.frame(sample = 1:length(25:155), 
                    #outcome = sample_metadata$stress,
                    batch = c(rep(1, 12), rep(2,15), rep(3,28), rep(4,12), rep(5,32), rep(6,32)))
rownames(pheno) <- colnames(AllProteinsList$Proteins)
mod = model.matrix(~as.factor(outcome), data=pheno)
mod = model.matrix(~1, data=pheno)
mod0 = model.matrix(~1,data=pheno)
pheno$outcome <- gsub(pattern = "R", replacement = "C", x = pheno$outcome)
# combat not useful
combat_edata = ComBat(dat=AllProteinsList$Proteins, batch=pheno$batch, mod=mod, par.prior=TRUE, prior.plots = FALSE, mean.only=FALSE)
write.table(x = combat_edata ,file = "F:/ATLAS_Pra/Protein_module/0.Process/ProtAtlas.Analyses.Comp.batch.txt", sep = "\t", quote = F, col.names = T, row.names = T)
# limma
pheno <- data.frame(sample = 1:length(25:155), 
                    outcome = paste0(sample_metadata$location, "_", sample_metadata$stress),
                    batch = c(rep(1, 12), rep(2,15), rep(3,28), rep(4,12), rep(5,32), rep(6,32)))
pheno$outcome <- gsub(pattern = "_R", replacement = "_C", pheno$outcome)
rownames(pheno) <- colnames(AllProteinsList$Proteins)
mm = model.matrix(~as.factor(outcome), data=pheno)
mm = model.matrix(~1, data=pheno)
mod0 = model.matrix(~1,data=pheno)
mat <- limma::removeBatchEffect(AllProteinsList$Proteins, batch=pheno$batch, design=mm)
write.table(x = mat ,file = "F:/ATLAS_Pra/Protein_module/0.Process/Figures/MOFA/StressComp/ProtAtlas.Analyses.Comp.LimmaOutcomeLocationStress.batch.txt", sep = "\t", quote = F, col.names = T, row.names = T)

AllProteinsList.vars.HS <- apply(AllProteinsList$Proteins[,c(grep("HS",colnames(AllProteinsList$Proteins)))], 1, var)
AllProteinsList.vars.HS <- names(which(AllProteinsList.vars.HS > 0.2))
AllProteinsList.vars.UV <- apply(AllProteinsList$Proteins[,c(grep("UV",colnames(AllProteinsList$Proteins)))], 1, var)
AllProteinsList.vars.UV <- names(which(AllProteinsList.vars.UV > 0.2))
AllProteinsList.vars <- unique(AllProteinsList.vars.HS, AllProteinsList.vars.UV)
AllProteinsList$Proteins <- AllProteinsList$Proteins[AllProteinsList.vars,]
AllProteinsList$Proteins <- as.matrix(AllProteinsList$Proteins)

AllProteinsList.vars <- apply(AllProteinsList$Proteins, 1, var)
AllProteinsList$Proteins <- AllProteinsList$Proteins[names(which(AllProteinsList.vars > 0.2)),] 
AllProteinsList$Proteins <- as.matrix(AllProteinsList$Proteins)


groups = c(rep("Total", 27), rep("Nucleus", 40), rep("Chloroplast", 64))
groups <- c(rep("HS",12), rep("UV",15), rep("HS", 28), rep("UV", 12), rep("HS", 32), rep("UV", 32))

MOFAgrouped <- create_mofa(data = AllProteinsList, groups = groups, save_metadata = TRUE)

## 3.2.2 Prepare MOFA and run

plot_data_overview(MOFAgrouped)
sample_metadata <- data.frame(
  sample = colnames(AllProteinsList$Proteins),
  stress = c(rep("C", 4), rep("S",8), rep("C",3), rep("S",9), rep("R",3),
             rep("C", 4), rep("S", 16), rep("R", 8),
             rep("C", 4), rep("S", 8), 
             rep("C", 8), rep("S", 24), 
             rep("C", 8), rep("S", 24)),
  treatment = c(rep("C", 4), rep("HS",8), rep("C",3), rep("UV",9), rep("UVR",3),
                rep("C", 4), rep("HS", 16), rep("HSR", 8),
                rep("C", 4), rep("UV", 8), 
                rep("C", 8), rep("HS", 24), 
                rep("C", 8), rep("UV", 24)),
  time = c(rep("T0", 4), rep("T1", 4), rep("T2", 4),
           rep("T0",3), rep("T1", 3), rep("T2", 3), rep("T3", 3) , rep("R",3),
           rep("T0", 4), rep("T1", 4), rep("T2", 4), rep("T3", 4), rep("T4", 4), rep("R", 8),
           rep("T0", 4), rep("T1", 4), rep("T3", 4), 
           rep("T0", 8), rep("T1", 8), rep("T2", 8), rep("T3", 8), 
           rep("T0", 8),  rep("T1", 8), rep("T2", 8), rep("T3", 8)),
  group = groups,
  location = c(rep("Total", 27), rep("Nucleus", 40), rep("Chloroplast", 64)),
  sample_number = c(1:ncol(AllProteinsList$Proteins))
)
sample_metadata <- data.frame(
  sample = colnames(AllProteinsList$Proteins),
  stress = c(
             #rep("C", 4), rep("S",8), rep("C",3), rep("S",9), rep("R",3),
             rep("C", 4), rep("S", 16), rep("R", 8),
             rep("C", 4), rep("S", 8) 
             #rep("C", 8), rep("S", 24), 
             #rep("C", 8), rep("S", 24)
             ),
  treatment = c(
             #rep("C", 4), rep("HS",8), rep("C",3), rep("UV",9), rep("UVR",3),
               rep("C", 4), rep("HS", 16), rep("HSR", 8),
                rep("C", 4), rep("UV", 8) 
                #rep("C", 8), rep("HS", 24), 
                #rep("C", 8), rep("UV", 24)
               ),
  time = c(
    #rep("T0", 4), rep("T1", 4), rep("T2", 4),
           #rep("T0",3), rep("T1", 3), rep("T2", 3), rep("T3", 3) , rep("R",3),
           rep("T0", 4), rep("T1", 4), rep("T2", 4), rep("T3", 4), rep("T4", 4), rep("R", 8),
           rep("T0", 4), rep("T1", 4), rep("T3", 4) 
           #rep("T0", 8), rep("T1", 8), rep("T2", 8), rep("T3", 8), 
           #rep("T0", 8),  rep("T1", 8), rep("T2", 8), rep("T3", 8)
           ),
  group = groups,
  #location = c(rep("Total", 27), rep("Nucleus", 40), rep("Chloroplast", 64)),
  sample_number = c(1:ncol(AllProteinsList$Proteins))
)
data_opts <- get_default_data_options(MOFAgrouped)
model_opts <- get_default_model_options(MOFAgrouped)
model_opts$num_factors <- 30
train_opts <- get_default_training_options(MOFAgrouped)
train_opts$maxiter <- "100000"
train_opts$convergence_mode <- "slow"
train_opts$drop_factor_threshold <- 0.005
train_opts$seed <- 400
data_opts$scale_views <- FALSE
outfileGrouped <- "F:/ATLAS_Pra/Protein_module/0.Process/Figures/MOFA/StressComp/MOFA2grouped_StressComp_log101_V02_400_false.hdf5"

MOFAgrouped <- prepare_mofa(MOFAgrouped, data_options = data_opts, model_options = model_opts, training_options = train_opts)
MOFAgrouped.trained <- run_mofa(MOFAgrouped, outfile = outfileGrouped)

MOFAgrouped <- prepare_mofa(MOFAgrouped, data_options = data_opts, model_options = model_opts, training_options = train_opts)
MOFAgrouped.trained <- run_mofa(MOFAgrouped, outfile = outfileGrouped)

MOFAgrouped <- load_model("F:/ATLAS_Pra/Protein_module/0.Process/Figures/MOFA/StressComp/woBatch/MOFA2grouped_StressComp_30LFs_log101_V02_400_false.hdf5")
MOFAgrouped.trained <- MOFAgrouped

## 3.2.3 Plot General

plot_factor_cor(MOFAgrouped.trained)
plot_variance_explained(MOFAgrouped.trained, x = "group", y = "factor")
plot_variance_explained(MOFAgrouped.trained, x = "group", y = "factor", plot_total = T)

MOFAgrouped.trained@cache$variance_explained
MOFAgrouped.trained@cache$variance_explained$r2_per_factor$FU
MOFAgrouped.trained@cache$variance_explained$r2_per_factor$HS
MOFAgrouped.trained@cache$variance_explained$r2_per_factor$UV

a <- plot_variance_explained(MOFAgrouped.trained, x="view", y="factor", factors = c(1:11))
b <- plot_variance_explained(MOFAgrouped.trained, x="group", y="factor", plot_total = T, factors = c(1:11))[[2]]
PlotVariance <- grid.arrange(b, a, ncol = 1, nrow = 2)

## 3.2.4 LFs Inspection

# Single

MOFAgrouped.trained@samples_metadata <- sample_metadata

plot_factor(MOFAgrouped.trained, 
            factor = c(1:11),
            color_by = "time",
            shape_by = "location", dot_size = 4
)

# Multiple: 2 location between both stress shared, UV location unique, HS nuclues;chloroplasto and UV only chloroplast

plot_factors(MOFAgrouped.trained, 
             factors = c(2, 3),
             color_by = "time", dot_size = 8, shape_by = "location"
)

plot_factors(MOFAgrouped.trained, 
             factors = c(4, 5),
             color_by = "time", dot_size = 1, shape_by = "location"
)

plot_factors(MOFAgrouped.trained, 
             factors = c(7, 9),
             color_by = "time", dot_size = 8, shape_by = "location"
)

# Weights of biologically relevant LFs

plot_weights(MOFAgrouped.trained,
             view = "Proteins",
             factor = c(7,9),
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

plot_weights(MOFAgrouped.trained,
             view = "Proteins",
             factor = 2,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

plot_weights(MOFAgrouped.trained,
             view = "Proteins",
             factor = 3,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

plot_weights(MOFAgrouped.trained,
             view = "Proteins",
             factor = 4,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

LF1<- plot_top_weights(MOFAgrouped.trained,
                        view = "Proteins",
                        factor = 1,
                        nfeatures = 10,     # Number of features to highlight
                        scale = T,          # Scale weights from -1 to 1
                        abs = F  )           # Take the absolute value?)

LF2 <- plot_top_weights(MOFAgrouped.trained,
                        view = "Proteins",
                        factor = 2,
                        nfeatures = 10,     # Number of features to highlight
                        scale = T,          # Scale weights from -1 to 1
                        abs = F  )           # Take the absolute value?)

LF3 <- plot_top_weights(MOFAgrouped.trained,
                        view = "Proteins",
                        factor = 4,
                        nfeatures = 10,     # Number of features to highlight
                        scale = T,          # Scale weights from -1 to 1
                        abs = F  )           # Take the absolute value?)

LF4 <- plot_top_weights(MOFAgrouped.trained,
                        view = "Proteins",
                        factor = 5,
                        nfeatures = 10,     # Number of features to highlight
                        scale = T,          # Scale weights from -1 to 1
                        abs = F  )           # Take the absolute value?)


gridExtra::grid.arrange(LF1, LF2,LF3 ,LF4, ncol = 2, nrow = 2)



PROTAS.annotation$Symbol[grep("all-04−2099895", PROTAS.annotation$ProtID)]

# Enrichment Analyses #
# Positive and neg splitted
# Mercator bins, PS age, gene family founder PS

# Mercator

Bin_db

# Gene Age

PROTIS <- read.delim(file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/Protein_module/ProtAtlas.Analyses.txt", header = T, sep = "\t")
Ages <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/4.GenEra/Output/Pinus_radiata_proteome/3347_gene_ages.tsv", header = T, sep = "\t")
Ages$X.gene <- do.call(rbind, strsplit(x = Ages$X.gene, split = ".", fixed = T))[,1] 
PROTIS$Phylostratum <- NA
PROTIS <- PROTIS[PROTIS$ProtID %in% Ages$X.gene,]

for(i in 1:nrow(PROTIS)){
  age <- unique(Ages$rank[which(Ages$X.gene == PROTIS$ProtID[i])])
  if(length(age) == 1){
    PROTIS$Phylostratum[i] <- age
  }else if(length(age) > 1){
    PROTIS$Phylostratum[i] <- min(age)
  }else if(length(age) == 0){
    stop("Something went wrong boy ... \n")
  }
}

PROTIS <- PROTIS[,c(ncol(PROTIS),1:ncol(PROTIS)-1)]
colnames(PROTIS)[2] <- "GeneID"
PROTIS <- PROTIS[,-3]
PROTIS <- PROTIS[!is.na(PROTIS$Phylostratum),]

PROTIS$Phylostratum <- gsub("15", "12", PROTIS$Phylostratum, fixed = T) 
PROTIS$Phylostratum <- gsub("13", "12", PROTIS$Phylostratum, fixed = T) 
PROTIS$Phylostratum <- gsub("14", "12", PROTIS$Phylostratum, fixed = T) 

PS_db <- list()
PS_names <- levels(as.factor(PROTIS$Phylostratum))

for(i in PS_names){
  PS_db[[i]] <- PROTIS$GeneID[PROTIS$Phylostratum == i]
}

names(PS_db) <- paste0("PS", names(PS_db))
PS_db

# Gene family founder

PROTIS <- read.delim(file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/Protein_module/ProtAtlas.Analyses.txt", header = T, sep = "\t")
Ages <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/4.GenEra/Output/Pinus_radiata_proteome/3347_founder_events.tsv", header = T, sep = "\t")
Ages$rank <- gsub("15", "12", Ages$rank, fixed = T) 
Ages$rank <- gsub("13", "12", Ages$rank, fixed = T) 
Ages$rank <- gsub("14", "12", Ages$rank, fixed = T) 
Ages$FamilyID <- paste0("Founder",rownames(Ages),".", Ages$rank, ".",Ages$family_size)
Ages <- tidyr::separate_longer_delim(Ages, X.gene_family, delim = ",")
Ages$X.gene <- do.call(rbind, strsplit(x = Ages$X.gene, split = ".", fixed = T))[,1] 
PROTIS$Phylostratum <- NA
PROTIS <- PROTIS[PROTIS$ProtID %in% Ages$X.gene,]

for(i in 1:nrow(PROTIS)){
  age <- unique(Ages$rank[which(Ages$X.gene == PROTIS$ProtID[i])])
  if(length(age) == 1){
    PROTIS$Phylostratum[i] <- age
  }else if(length(age) > 1){
    PROTIS$Phylostratum[i] <- min(age)
  }else if(length(age) == 0){
    stop("Something went wrong boy ... \n")
  }
}

PROTIS <- PROTIS[,c(ncol(PROTIS),1:ncol(PROTIS)-1)]
colnames(PROTIS)[2] <- "GeneID"
PROTIS <- PROTIS[,-3]
PROTIS <- PROTIS[!is.na(PROTIS$Phylostratum),]

PSF_db <- list()
PSF_names <- levels(as.factor(PROTIS$Phylostratum))

for(i in PSF_names){
  PSF_db[[i]] <- PROTIS$GeneID[PROTIS$Phylostratum == i]
}

names(PSF_db) <- paste0("PS", names(PSF_db))
PSF_db

# Compute enrichments; LF1-4; Pos and neg; Mercator, Gene Age and Family foundation

Bin.paths <- list_to_matrix(Bin_db)
PS.paths <- list_to_matrix(PS_db)
PSF.paths <- list_to_matrix(PSF_db)


enrichment.parametric.Bin.all <- run_enrichment(MOFAgrouped.trained,
                                                view = "Proteins", factors = c(1,2,4,5),
                                                feature.sets = t(Bin.paths),
                                                sign = "all", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                #nperm = 5000
)


enrichment.parametric.PS.all <- run_enrichment(MOFAgrouped.trained,
                                               view = "Proteins", factors = c(1,2,4,5),
                                               feature.sets = t(PS.paths),
                                               sign = "all", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                               # statistical.test = "permutation",
                                               min.size = 5, alpha = 0.1
                                               #nperm = 5000
)

enrichment.parametric.PSF.all <- run_enrichment(MOFAgrouped.trained,
                                                view = "Proteins", factors = c(1,2,4,5),
                                                feature.sets = t(PSF.paths),
                                                sign = "all", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, alpha = 0.1
                                                #nperm = 5000
)

enrichment.parametric.Bin.pos <- run_enrichment(MOFAgrouped.trained,
                                                view = "Proteins", factors = c(1,2,4,5),
                                                feature.sets = t(Bin.paths),
                                                sign = "pos", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                #nperm = 5000
)


enrichment.parametric.PS.pos <- run_enrichment(MOFAgrouped.trained,
                                               view = "Proteins", factors = c(1,2,4,5),
                                               feature.sets = t(PS.paths),
                                               sign = "pos", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                               # statistical.test = "permutation",
                                               min.size = 5, alpha = 0.1
                                               #nperm = 5000
)

enrichment.parametric.PSF.pos <- run_enrichment(MOFAgrouped.trained,
                                                view = "Proteins", factors = c(1,2,4,5),
                                                feature.sets = t(PSF.paths),
                                                sign = "pos", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, alpha = 0.1
                                                #nperm = 5000
)

enrichment.parametric.Bin.neg <- run_enrichment(MOFAgrouped.trained,
                                                view = "Proteins", factors = c(1,2,4,5),
                                                feature.sets = t(Bin.paths),
                                                sign = "neg", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                #nperm = 5000
)


enrichment.parametric.PS.neg <- run_enrichment(MOFAgrouped.trained,
                                               view = "Proteins", factors = c(1,2,4,5),
                                               feature.sets = t(PS.paths),
                                               sign = "neg", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                               # statistical.test = "permutation",
                                               min.size = 5, alpha = 0.1
                                               #nperm = 5000
)

enrichment.parametric.PSF.neg <- run_enrichment(MOFAgrouped.trained,
                                                view = "Proteins", factors = c(1,2,4,5),
                                                feature.sets = t(PSF.paths),
                                                sign = "neg", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, alpha = 0.1
                                                #nperm = 5000
)

## plot enrich

order <- sort(rownames(enrichment.parametric.Bin.all$pval.adj))

Bin.df <- t(enrichment.parametric.Bin.all$pval.adj)[,order]

for(i in 1:nrow(Bin.df)){
  for(j in 1:ncol(Bin.df)){
    if(Bin.df[i,j] > 0.1 ){
      Bin.df[i,j] <- 1
    }
  }
}

Bin.df <- -log10(Bin.df)

for(i in 1:nrow(Bin.df)){
  for(j in 1:ncol(Bin.df)){
    if(Bin.df[i,j] > 10 ){
      Bin.df[i,j] <- 10
    }
  }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Greens"))(500)

ht.Bin <- Heatmap(Bin.df, name = "Bins.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

order <- sort(rownames(enrichment.parametric.PS.all$pval.adj))

PS.df <- t(enrichment.parametric.PS.all$pval.adj)[,order]

for(i in 1:nrow(PS.df)){
  for(j in 1:ncol(PS.df)){
    if(PS.df[i,j] > 0.1 ){
      PS.df[i,j] <- 1
    }
  }
}

PS.df <- -log10(PS.df)

for(i in 1:nrow(PS.df)){
  for(j in 1:ncol(PS.df)){
    if(PS.df[i,j] > 10 ){
      PS.df[i,j] <- 10
    }
  }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Purples"))(100)

ht.PS <- Heatmap(PS.df, name = "PS.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

order <- sort(rownames(enrichment.parametric.PSF.all$pval.adj))

PSF.df <- t(enrichment.parametric.PSF.all$pval.adj)[,order]

for(i in 1:nrow(PSF.df)){
  for(j in 1:ncol(PSF.df)){
    if(PSF.df[i,j] > 0.1 ){
      PSF.df[i,j] <- 1
    }
  }
}

PSF.df <- -log10(PSF.df)

for(i in 1:nrow(PSF.df)){
  for(j in 1:ncol(PSF.df)){
    if(PSF.df[i,j] > 10 ){
      PSF.df[i,j] <- 10
    }
  }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100)

ht.PSF <- Heatmap(PSF.df, name = "PSF.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

draw(ht.Bin + ht.PS + ht.PSF, ht_gap = unit(1, "cm"))
draw(ht.Bin %v% ht.PS %v% ht.PSF, ht_gap = unit(1, "cm"))

# nVenn or Upset same gene family: not that useful (all Stress Total Proteins comes from same 803 founder events)
factors <- get_factors(object =MOFAgrouped.trained, as.data.frame = T)
weights <- get_weights(object = MOFAgrouped.trained, as.data.frame = T)
write.table(factors, "F:/ATLAS_Pra/Protein_module/0.Process/Figures/MOFA/StressComp/FACTORS.txt", col.names = T, sep = "\t", quote = F)
VennFamily <- Ages[,c(6,5)]
LF1.50 <- names(sort(abs(deframe(weights[which(weights$factor == "Factor2"), c(1,3)])), decreasing = T)[1:500])
LF1.50 <-  VennFamily$FamilyID[match(x =  LF1.50,table = VennFamily$X.gene)]
LF2.50 <- names(sort(abs(deframe(weights[which(weights$factor == "Factor6"), c(1,3)])), decreasing = T)[1:500])
LF2.50 <-  VennFamily$FamilyID[match(x =  LF2.50,table = VennFamily$X.gene)]
LF3.50 <- names(sort(abs(deframe(weights[which(weights$factor == "Factor7"), c(1,3)])), decreasing = T)[1:500])
LF3.50 <-  VennFamily$FamilyID[match(x =  LF3.50,table = VennFamily$X.gene)]
LF4.50 <- names(sort(abs(deframe(weights[which(weights$factor == "Factor9"), c(1,3)])), decreasing = T)[1:500])
LF4.50 <-  VennFamily$FamilyID[match(x =  LF4.50,table = VennFamily$X.gene)]
dataList <- list(LF2 = LF1.50, LF6 = LF2.50, LF7 = LF3.50, LF9 = LF4.50)
#myV <- nVennR::plotVenn(data, nCycles = 20000, opacity = 0.2, borderWidth = 3, systemShow = T, fontScale = 2, outFile = paste0(output.path, "nVennR.svg"))
upset(fromList(dataList), sets = names(dataList), order.by = "freq", 
      keep.order = TRUE, empty.intersections = F, mb.ratio = c(0.6,0.4), sets.bar.color = "navajowhite3")


#### 3.3 MOFA multigroup Populations (only chloroplast for both UV and heat stress) ####

## 3.2.1 Format data and import to MOFA

ProtAtlas.Pop <- ProtAtlas[,c(1,92:155)]
AllProteinsList <- list(Proteins = ProtAtlas.Pop)
AllProteinsList$Proteins[,2:ncol(AllProteinsList$Proteins)] <- log10(AllProteinsList$Proteins[,2:ncol(AllProteinsList$Proteins)] + 1)
rownames(AllProteinsList$Proteins) <- AllProteinsList$Proteins$ProtID
AllProteinsList$Proteins <- AllProteinsList$Proteins[,-1]

# combat: not useful
pheno <- data.frame(sample = 1:length(92:155), 
                    outcome = paste0(sample_metadata$population, "_", sample_metadata$stress, "_", sample_metadata$time),
                    batch = c(rep(1,32), rep(2,32)))
rownames(pheno) <- colnames(AllProteinsList$Proteins)
mod = model.matrix(~as.factor(outcome), data=pheno)
mod = model.matrix(~1, data=pheno)
mod0 = model.matrix(~1,data=pheno)

combat_edata = ComBat(dat=AllProteinsList$Proteins, batch=pheno$batch, mod=mod, par.prior=TRUE, prior.plots = FALSE, mean.only=FALSE)
write.table(x = combat_edata ,file = "F:/ATLAS_Pra/Protein_module/0.Process/Figures/MOFA/StressPop/ProtAtlas.Analyses.Pop.OutcomeAll.batch.txt", sep = "\t", quote = F, col.names = T, row.names = T)

# limma
pheno <- data.frame(sample = 1:length(92:155), 
                    outcome = paste0(sample_metadata$population, "_", sample_metadata$stress, "_", sample_metadata$time),
                    batch = c(rep(1,32), rep(2,32)))
rownames(pheno) <- colnames(AllProteinsList$Proteins)
mm = model.matrix(~as.factor(outcome), data=pheno)
mm = model.matrix(~1, data=pheno)
mod0 = model.matrix(~1,data=pheno)
mat <- limma::removeBatchEffect(AllProteinsList$Proteins, batch=pheno$batch, design=mm)
write.table(x = mat ,file = "F:/ATLAS_Pra/Protein_module/0.Process/Figures/MOFA/StressPop/ProtAtlas.Analyses.Pop.LimmaOutcomeAll.batch.txt", sep = "\t", quote = F, col.names = T, row.names = T)


AllProteinsList.vars.HS <- apply(AllProteinsList$Proteins[,c(grep("HS",colnames(AllProteinsList$Proteins)))], 1, var)
AllProteinsList.vars.HS <- names(which(AllProteinsList.vars.HS > 0))
AllProteinsList.vars.UV <- apply(AllProteinsList$Proteins[,c(grep("UV",colnames(AllProteinsList$Proteins)))], 1, var)
AllProteinsList.vars.UV <- names(which(AllProteinsList.vars.UV > 0))
AllProteinsList.vars <- unique(AllProteinsList.vars.HS, AllProteinsList.vars.UV)
AllProteinsList$Proteins <- AllProteinsList$Proteins[AllProteinsList.vars,]
AllProteinsList$Proteins <- as.matrix(AllProteinsList$Proteins)

AllProteinsList.vars.HS <- apply(AllProteinsList$Proteins[,c(grep("PopE", groups))], 1, var)
AllProteinsList.vars.HS <- names(which(AllProteinsList.vars.HS > 0))
AllProteinsList.vars.UV <- apply(AllProteinsList$Proteins[,c(grep("PopT", groups))], 1, var)
AllProteinsList.vars.UV <- names(which(AllProteinsList.vars.UV > 0))
AllProteinsList.vars <- unique(AllProteinsList.vars.HS, AllProteinsList.vars.UV)
AllProteinsList$Proteins <- AllProteinsList$Proteins[AllProteinsList.vars,]
AllProteinsList$Proteins <- as.matrix(AllProteinsList$Proteins)

AllProteinsList.vars <- apply(AllProteinsList$Proteins, 1, var)
AllProteinsList$Proteins <- AllProteinsList$Proteins[names(which(AllProteinsList.vars > 0.2)),] 
AllProteinsList$Proteins <- as.matrix(AllProteinsList$Proteins)


groups = c(rep("PopE", 4), rep("PopT", 4), rep("PopE", 4), rep("PopT", 4), rep("PopE", 4), rep("PopT", 4), rep("PopE", 4), rep("PopT", 4),
           rep("PopE", 4), rep("PopT", 4), rep("PopE", 4), rep("PopT", 4), rep("PopE", 4), rep("PopT", 4), rep("PopE", 4), rep("PopT", 4))
groups <- c(rep("HS",32), rep("UV",32))

MOFAgrouped <- create_mofa(data = AllProteinsList, groups = groups, save_metadata = TRUE)

## 3.2.2 Prepare MOFA and run

plot_data_overview(MOFAgrouped)
sample_metadata <- data.frame(
  sample = colnames(AllProteinsList$Proteins),
  stress = c(
             rep("C", 8), rep("S", 24), 
             rep("C", 8), rep("S", 24)),
  treatment = c(
                rep("C", 8), rep("HS", 24), 
                rep("C", 8), rep("UV", 24)),
  time = c(
           rep("T0", 8), rep("T1", 8), rep("T2", 8), rep("T3", 8), 
           rep("T0", 8),  rep("T1", 8), rep("T2", 8), rep("T3", 8)),
  group = groups,
  population = c(rep("PopE", 4), rep("PopT", 4), rep("PopE", 4), rep("PopT", 4), rep("PopE", 4), rep("PopT", 4), rep("PopE", 4), rep("PopT", 4),
                 rep("PopE", 4), rep("PopT", 4), rep("PopE", 4), rep("PopT", 4), rep("PopE", 4), rep("PopT", 4), rep("PopE", 4), rep("PopT", 4)),
  sample_number = c(1:ncol(AllProteinsList$Proteins))
)

data_opts <- get_default_data_options(MOFAgrouped)
model_opts <- get_default_model_options(MOFAgrouped)
model_opts$num_factors <- 10
train_opts <- get_default_training_options(MOFAgrouped)
train_opts$maxiter <- "100000"
train_opts$convergence_mode <- "slow"
train_opts$drop_factor_threshold <- 0.005
train_opts$seed <- 400
data_opts$scale_views <- FALSE
outfileGrouped <- "F:/ATLAS_Pra/Protein_module/0.Process/Figures/MOFA/StressComp/MOFA2grouped_StressComp_log101_V02_400_false.hdf5"

MOFAgrouped <- prepare_mofa(MOFAgrouped, data_options = data_opts, model_options = model_opts, training_options = train_opts)
MOFAgrouped.trained <- run_mofa(MOFAgrouped, outfile = outfileGrouped)

MOFAgrouped <- prepare_mofa(MOFAgrouped, data_options = data_opts, model_options = model_opts, training_options = train_opts)
MOFAgrouped.trained <- run_mofa(MOFAgrouped, outfile = outfileGrouped)

MOFAgrouped <- load_model("F:/ATLAS_Pra/Protein_module/0.Process/Figures/MOFA/StressPop/MOFA2grouped_StressPop_Popgroups_CombatOutcomePopStressTime_10LFs_log101_V0_400_false.hdf5")
MOFAgrouped.trained <- MOFAgrouped

## 3.2.3 Plot General

plot_factor_cor(MOFAgrouped.trained)
plot_variance_explained(MOFAgrouped.trained, x = "group", y = "factor")
plot_variance_explained(MOFAgrouped.trained, x = "group", y = "factor", plot_total = T)

MOFAgrouped.trained@cache$variance_explained

a <- plot_variance_explained(MOFAgrouped.trained, x="view", y="factor", factors = c(2,3,4,5,6))
b <- plot_variance_explained(MOFAgrouped.trained, x="group", y="factor", plot_total = T, factors = c(2,3,4,5,6))[[2]]
PlotVariance <- grid.arrange(b, a, ncol = 1, nrow = 2)

## 3.2.4 LFs Inspection

# Single

MOFAgrouped.trained@samples_metadata <- sample_metadata

plot_factor(MOFAgrouped.trained, 
            factor = c(1,2,3,4,5,6),
            color_by = "time",
            shape_by = "treatment", dot_size = 4
)

# Multiple: 2 location between both stress shared, UV location unique, HS nuclues;chloroplasto and UV only chloroplast

plot_factors(MOFAgrouped.trained, 
             factors = c(1,2,3,4,5,6),
             color_by = "time", dot_size = 8, shape_by = "population"
)

plot_factors(MOFAgrouped.trained, 
             factors = c(2,5),
             color_by = "time", dot_size = 8, shape_by = "population"
)

plot_factors(MOFAgrouped.trained, 
             factors = c(7, 9),
             color_by = "time", dot_size = 8, shape_by = "location"
)

# Weights of biologically relevant LFs

plot_weights(MOFAgrouped.trained,
             view = "Proteins",
             factor = c(7,9),
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

plot_weights(MOFAgrouped.trained,
             view = "Proteins",
             factor = 2,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

plot_weights(MOFAgrouped.trained,
             view = "Proteins",
             factor = 3,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

plot_weights(MOFAgrouped.trained,
             view = "Proteins",
             factor = 4,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

LF1<- plot_top_weights(MOFAgrouped.trained,
                       view = "Proteins",
                       factor = 2,
                       nfeatures = 10,     # Number of features to highlight
                       scale = T,          # Scale weights from -1 to 1
                       abs = F  )           # Take the absolute value?)

LF2 <- plot_top_weights(MOFAgrouped.trained,
                        view = "Proteins",
                        factor = 3,
                        nfeatures = 10,     # Number of features to highlight
                        scale = T,          # Scale weights from -1 to 1
                        abs = F  )           # Take the absolute value?)

LF3 <- plot_top_weights(MOFAgrouped.trained,
                        view = "Proteins",
                        factor = 5,
                        nfeatures = 10,     # Number of features to highlight
                        scale = T,          # Scale weights from -1 to 1
                        abs = F  )           # Take the absolute value?)

LF4 <- plot_top_weights(MOFAgrouped.trained,
                        view = "Proteins",
                        factor = 6,
                        nfeatures = 10,     # Number of features to highlight
                        scale = T,          # Scale weights from -1 to 1
                        abs = F)           # Take the absolute value?)


gridExtra::grid.arrange(LF1, LF2,LF3 ,LF4, ncol = 2, nrow = 2)



PROTAS.annotation$Symbol[grep("all-04-292070", PROTAS.annotation$ProtID)]

# Enrichment Analyses #
# Positive and neg splitted
# Mercator bins, PS age, gene family founder PS

# Mercator

Bin_db

# Gene Age

PROTIS <- read.delim(file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/Protein_module/ProtAtlas.Analyses.txt", header = T, sep = "\t")
Ages <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/4.GenEra/Output/Pinus_radiata_proteome/3347_gene_ages.tsv", header = T, sep = "\t")
Ages$X.gene <- do.call(rbind, strsplit(x = Ages$X.gene, split = ".", fixed = T))[,1] 
PROTIS$Phylostratum <- NA
PROTIS <- PROTIS[PROTIS$ProtID %in% Ages$X.gene,]

for(i in 1:nrow(PROTIS)){
  age <- unique(Ages$rank[which(Ages$X.gene == PROTIS$ProtID[i])])
  if(length(age) == 1){
    PROTIS$Phylostratum[i] <- age
  }else if(length(age) > 1){
    PROTIS$Phylostratum[i] <- min(age)
  }else if(length(age) == 0){
    stop("Something went wrong boy ... \n")
  }
}

PROTIS <- PROTIS[,c(ncol(PROTIS),1:ncol(PROTIS)-1)]
colnames(PROTIS)[2] <- "GeneID"
PROTIS <- PROTIS[,-3]
PROTIS <- PROTIS[!is.na(PROTIS$Phylostratum),]

PROTIS$Phylostratum <- gsub("15", "12", PROTIS$Phylostratum, fixed = T) 
PROTIS$Phylostratum <- gsub("13", "12", PROTIS$Phylostratum, fixed = T) 
PROTIS$Phylostratum <- gsub("14", "12", PROTIS$Phylostratum, fixed = T) 

PS_db <- list()
PS_names <- levels(as.factor(PROTIS$Phylostratum))

for(i in PS_names){
  PS_db[[i]] <- PROTIS$GeneID[PROTIS$Phylostratum == i]
}

names(PS_db) <- paste0("PS", names(PS_db))
PS_db

# Gene family founder

PROTIS <- read.delim(file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/Protein_module/ProtAtlas.Analyses.txt", header = T, sep = "\t")
Ages <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/4.GenEra/Output/Pinus_radiata_proteome/3347_founder_events.tsv", header = T, sep = "\t")
Ages$rank <- gsub("15", "12", Ages$rank, fixed = T) 
Ages$rank <- gsub("13", "12", Ages$rank, fixed = T) 
Ages$rank <- gsub("14", "12", Ages$rank, fixed = T) 
Ages$FamilyID <- paste0("Founder",rownames(Ages),".", Ages$rank, ".",Ages$family_size)
Ages <- tidyr::separate_longer_delim(Ages, X.gene_family, delim = ",")
Ages$X.gene <- do.call(rbind, strsplit(x = Ages$X.gene, split = ".", fixed = T))[,1] 
PROTIS$Phylostratum <- NA
PROTIS <- PROTIS[PROTIS$ProtID %in% Ages$X.gene,]

for(i in 1:nrow(PROTIS)){
  age <- unique(Ages$rank[which(Ages$X.gene == PROTIS$ProtID[i])])
  if(length(age) == 1){
    PROTIS$Phylostratum[i] <- age
  }else if(length(age) > 1){
    PROTIS$Phylostratum[i] <- min(age)
  }else if(length(age) == 0){
    stop("Something went wrong boy ... \n")
  }
}

PROTIS <- PROTIS[,c(ncol(PROTIS),1:ncol(PROTIS)-1)]
colnames(PROTIS)[2] <- "GeneID"
PROTIS <- PROTIS[,-3]
PROTIS <- PROTIS[!is.na(PROTIS$Phylostratum),]

PSF_db <- list()
PSF_names <- levels(as.factor(PROTIS$Phylostratum))

for(i in PSF_names){
  PSF_db[[i]] <- PROTIS$GeneID[PROTIS$Phylostratum == i]
}

names(PSF_db) <- paste0("PS", names(PSF_db))
PSF_db

# Compute enrichments; LF1-4; Pos and neg; Mercator, Gene Age and Family foundation

Bin.paths <- list_to_matrix(Bin_db)
PS.paths <- list_to_matrix(PS_db)
PSF.paths <- list_to_matrix(PSF_db)


enrichment.parametric.Bin.all <- run_enrichment(MOFAgrouped.trained,
                                                view = "Proteins", factors = c(2,3,5,6),
                                                feature.sets = t(Bin.paths),
                                                sign = "all", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                #nperm = 5000
)


enrichment.parametric.PS.all <- run_enrichment(MOFAgrouped.trained,
                                               view = "Proteins", factors = c(2,3,5,6),
                                               feature.sets = t(PS.paths),
                                               sign = "all", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                               # statistical.test = "permutation",
                                               min.size = 5, alpha = 0.1
                                               #nperm = 5000
)

enrichment.parametric.PSF.all <- run_enrichment(MOFAgrouped.trained,
                                                view = "Proteins", factors = c(2,3,5,6),
                                                feature.sets = t(PSF.paths),
                                                sign = "all", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, alpha = 0.1
                                                #nperm = 5000
)

enrichment.parametric.Bin.pos <- run_enrichment(MOFAgrouped.trained,
                                                view = "Proteins", factors = c(2,3,5,6),
                                                feature.sets = t(Bin.paths),
                                                sign = "pos", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                #nperm = 5000
)


enrichment.parametric.PS.pos <- run_enrichment(MOFAgrouped.trained,
                                               view = "Proteins", factors = c(2,3,5,6),
                                               feature.sets = t(PS.paths),
                                               sign = "pos", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                               # statistical.test = "permutation",
                                               min.size = 5, alpha = 0.1
                                               #nperm = 5000
)

enrichment.parametric.PSF.pos <- run_enrichment(MOFAgrouped.trained,
                                                view = "Proteins", factors = c(2,3,5,6),
                                                feature.sets = t(PSF.paths),
                                                sign = "pos", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, alpha = 0.1
                                                #nperm = 5000
)

enrichment.parametric.Bin.neg <- run_enrichment(MOFAgrouped.trained,
                                                view = "Proteins", factors = c(2,3,5,6),
                                                feature.sets = t(Bin.paths),
                                                sign = "neg", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                #nperm = 5000
)


enrichment.parametric.PS.neg <- run_enrichment(MOFAgrouped.trained,
                                               view = "Proteins", factors = c(2,3,5,6),
                                               feature.sets = t(PS.paths),
                                               sign = "neg", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                               # statistical.test = "permutation",
                                               min.size = 5, alpha = 0.1
                                               #nperm = 5000
)

enrichment.parametric.PSF.neg <- run_enrichment(MOFAgrouped.trained,
                                                view = "Proteins", factors = c(2,3,5,6),
                                                feature.sets = t(PSF.paths),
                                                sign = "neg", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, alpha = 0.1
                                                #nperm = 5000
)

## plot enrich

order <- sort(rownames(enrichment.parametric.Bin.all$pval.adj))

Bin.df <- t(enrichment.parametric.Bin.all$pval.adj)[,order]

for(i in 1:nrow(Bin.df)){
  for(j in 1:ncol(Bin.df)){
    if(Bin.df[i,j] > 0.1 ){
      Bin.df[i,j] <- 1
    }
  }
}

Bin.df <- -log10(Bin.df)

for(i in 1:nrow(Bin.df)){
  for(j in 1:ncol(Bin.df)){
    if(Bin.df[i,j] > 10 ){
      Bin.df[i,j] <- 10
    }
  }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Greens"))(500)

ht.Bin <- Heatmap(t(Bin.df), name = "Bins.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

order <- sort(rownames(enrichment.parametric.PS.all$pval.adj))

PS.df <- t(enrichment.parametric.PS.all$pval.adj)[,order]

for(i in 1:nrow(PS.df)){
  for(j in 1:ncol(PS.df)){
    if(PS.df[i,j] > 0.1 ){
      PS.df[i,j] <- 1
    }
  }
}

PS.df <- -log10(PS.df)

for(i in 1:nrow(PS.df)){
  for(j in 1:ncol(PS.df)){
    if(PS.df[i,j] > 10 ){
      PS.df[i,j] <- 10
    }
  }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Purples"))(100)

ht.PS <- Heatmap(t(PS.df), name = "PS.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

order <- sort(rownames(enrichment.parametric.PSF.all$pval.adj))

PSF.df <- t(enrichment.parametric.PSF.all$pval.adj)[,order]

for(i in 1:nrow(PSF.df)){
  for(j in 1:ncol(PSF.df)){
    if(PSF.df[i,j] > 0.1 ){
      PSF.df[i,j] <- 1
    }
  }
}

PSF.df <- -log10(PSF.df)

for(i in 1:nrow(PSF.df)){
  for(j in 1:ncol(PSF.df)){
    if(PSF.df[i,j] > 10 ){
      PSF.df[i,j] <- 10
    }
  }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100)

ht.PSF <- Heatmap(t(PSF.df), name = "PSF.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

draw(ht.Bin + ht.PS + ht.PSF, ht_gap = unit(1, "cm"))
draw(ht.Bin %v% ht.PS %v% ht.PSF, ht_gap = unit(1, "cm"))

# nVenn or Upset same gene family: not that useful (all Stress Total Proteins comes from same 803 founder events)
factors <- get_factors(object =MOFAgrouped.trained, as.data.frame = T)
weights <- get_weights(object = MOFAgrouped.trained, as.data.frame = T)
write.table(factors, "F:/ATLAS_Pra/Protein_module/0.Process/Figures/MOFA/StressPop/FACTORS.txt", col.names = T, sep = "\t", quote = F)
VennFamily <- Ages[,c(6,5)]
LF1.50 <- names(sort(abs(deframe(weights[which(weights$factor == "Factor2"), c(1,3)])), decreasing = T)[1:150])
LF1.50 <-  VennFamily$FamilyID[match(x =  LF1.50,table = VennFamily$X.gene)]
LF2.50 <- names(sort(abs(deframe(weights[which(weights$factor == "Factor3"), c(1,3)])), decreasing = T)[1:150])
LF2.50 <-  VennFamily$FamilyID[match(x =  LF2.50,table = VennFamily$X.gene)]
LF3.50 <- names(sort(abs(deframe(weights[which(weights$factor == "Factor5"), c(1,3)])), decreasing = T)[1:150])
LF3.50 <-  VennFamily$FamilyID[match(x =  LF3.50,table = VennFamily$X.gene)]
LF4.50 <- names(sort(abs(deframe(weights[which(weights$factor == "Factor6"), c(1,3)])), decreasing = T)[1:150])
LF4.50 <-  VennFamily$FamilyID[match(x =  LF4.50,table = VennFamily$X.gene)]
dataList <- list(LF2 = LF1.50, LF3 = LF2.50, LF5 = LF3.50, LF6 = LF4.50)
#myV <- nVennR::plotVenn(data, nCycles = 20000, opacity = 0.2, borderWidth = 3, systemShow = T, fontScale = 2, outFile = paste0(output.path, "nVennR.svg"))
upset(fromList(dataList), sets = names(dataList), order.by = "freq", 
      keep.order = TRUE, empty.intersections = F, mb.ratio = c(0.6,0.4), sets.bar.color = "navajowhite3")
