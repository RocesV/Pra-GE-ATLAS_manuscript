##########################################
# Pra-GE-ATLAS - Transcriptomic Analyses #
##########################################

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
library(kissDE)
library(nVennR)
library(tximport)
library(DESeq2)
library(readxl)
library(stringr)
library(seqRFLP)
library(GenomicRanges)
library(GenomicFeatures)
library(hexbin)
library(pals)

#### 1. Import metadata, all events and expression ####

setwd("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/")
sample_metadata <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/sample_metadata.txt", header = T, sep = "\t")

## Expression

files <- list.files("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/5.Pinus_taeda/Expression/",recursive=TRUE)
files <- paste0("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/5.Pinus_taeda/Expression/",files)
files <- files[grep("quant.sf", files)]
# for loop sample names
names.list <- list()
for(i in 1:length(files)){
 Paired <- grep(pattern = "Paired_",x =  files[i], fixed = T)
 if(length(Paired) >  0){
   split1 <- strsplit(x = files[i], split = "Paired_", fixed = T)[[1]][2]
   split1 <- gsub(pattern = ".sra/", replacement = "", x = split1)
   split1 <- gsub(pattern = "/", replacement = "", x = split1)
   names.list[[i]] <- sample_metadata$SampleName[which(sample_metadata$SampleID == split1)]
 }else if(length(Paired) == 0){
   split1 <- strsplit(x = files[i], split = ".sra", fixed = T)[[1]][2]
   split1 <- gsub(pattern = "/", replacement = "", x = split1)
   names.list[[i]] <- sample_metadata$SampleName[which(sample_metadata$SampleID == split1)]
 }
}
names(files) <- make.unique(unlist(names.list))
txi.salmon <- tximport(files, type = "salmon", txOut = TRUE, )

# mean technical reps

Raw.counts <- as.data.frame(txi.salmon$counts)
Raw.counts.jic <- Raw.counts

Raw.counts.jic <- Raw.counts.jic[,-grep("HS", colnames(Raw.counts.jic))] # 21 22 23 24 55 56 57 58 59 60 61 62
Raw.counts.jic$Needle_HS_Control_Control <- rowMeans(cbind(Raw.counts$Needle_HS_Control_Control, Raw.counts$Needle_HS_Control_Control.1))
Raw.counts.jic$Needle_HS_Control_Control.1 <- rowMeans(cbind(Raw.counts$Needle_HS_Control_Control.2, Raw.counts$Needle_HS_Control_Control.3))
Raw.counts.jic$Needle_HS_Stress_1 <- rowMeans(cbind(Raw.counts$Needle_HS_Stress_1, Raw.counts$Needle_HS_Stress_1.1))
Raw.counts.jic$Needle_HS_Stress_1.1 <- rowMeans(cbind(Raw.counts$Needle_HS_Stress_1.2, Raw.counts$Needle_HS_Stress_1.3))
Raw.counts.jic$Needle_HS_Stress_2 <- rowMeans(cbind(Raw.counts$Needle_HS_Stress_2, Raw.counts$Needle_HS_Stress_2.1))
Raw.counts.jic$Needle_HS_Stress_2.1 <- rowMeans(cbind(Raw.counts$Needle_HS_Stress_2.2, Raw.counts$Needle_HS_Stress_2.3))
colnames(Raw.counts)
colnames(Raw.counts.jic)

Raw.counts.F <- Raw.counts.jic[,c(1:20, 136:137, 21:50, 138:141, 51:135)]
colnames(Raw.counts.F)

# combat-seq batch effect in raw counts removal

sample_metadata <- sample_metadata[-c(23, 24, 57, 58, 61, 62),]
sample_metadata$SampleName <- make.unique(sample_metadata$SampleName)

Raw.counts.F <- apply(Raw.counts.F, 2, function(x){
  round(x, digits = 0)
})

Raw.counts.F.filt <- Raw.counts.F[apply(Raw.counts.F, 1, function(x) sum(x == 0)) < ncol(Raw.counts.F) / 3, ] # filter genes that are not expressed in less than half of the samples
Raw.counts.F.batch <- ComBat_seq(counts = as.matrix(Raw.counts.F.filt),batch = as.numeric(as.factor(paste0(sample_metada$Batch1, "_", sample_metada$Batch2))) )
Raw.counts.F.batch.int <- apply(Raw.counts.F.batch, 2, as.integer)
which(is.na(Raw.counts.F.batch.int), arr.ind=TRUE)
# low counts filter

keep <- rowSums(Raw.counts.F.batch >= 10) >= 3
dds <- DESeqDataSetFromMatrix(countData = Raw.counts.F.batch,
                              colData = sample_metada,
                              design = ~ Group)
dds <- dds[keep,]
vsd<-vst(dds, blind = TRUE)
write.table(Raw.counts.F.batch[keep,], file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Final_RawCounts_lowFilt_Presence_batchCorrected.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(assay(vsd), file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Final_VST_lowFilt_Presence_batchCorrected.txt", sep = "\t", col.names = T, row.names = T, quote = F)

keep1 <- rowSums(Raw.counts.F.filt >= 10) >= 3
dds1 <- DESeqDataSetFromMatrix(countData = Raw.counts.F.filt,
                              colData = sample_metada,
                              design = ~ Group)
dds1 <- dds1[keep1,]
vsd1<-vst(dds1, blind = TRUE)
write.table(Raw.counts.F.filt[keep1,], file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Batch_NotCorrected/Final_RawCounts_lowFilt_Presence.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(assay(vsd1), file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Batch_NotCorrected/Final_VST_lowFilt_Presence.txt", sep = "\t", col.names = T, row.names = T, quote = F)


## Splicing

countsSplicing <- kissplice2counts("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/5.Pinus_taeda/Splicing/Pra-GE-ATLAS-AS-all_type1_kiss2genome.tsv", pairedEnd = FALSE, k2rg = TRUE, counts = 2)
sample_metada <- read.delim("clipboard", header = T)

# mean technical

splicing.names <- list()
for(i in 1:nrow(sample_metada)){
   if(sample_metada$Tech[i] == "PE"){
      splicing.names[[i]] <- c(paste0(sample_metada$NameID[i], "_F"),paste0(sample_metada$NameID[i], "_R"))
   }else if(sample_metada$Tech[i] == "SE"){
      splicing.names[[i]] <- sample_metada$NameID[i]
   }
}

dim(countsSplicing$countsEvents)
length(unlist(splicing.names))
colnames(countsSplicing$countsEvents) <- c("events.names", "events.length", unlist(splicing.names))
colnames(countsSplicing$psiInfo) <- c("events.names", unlist(splicing.names))

countsSplicing.jic <- countsSplicing
countsSplicing.jic$countsEvents <- countsSplicing.jic$countsEvents[,-grep("HS", colnames(countsSplicing.jic$countsEvents))] # 37 38 39 40 41 42 43 44 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120
countsSplicing.jic$countsEvents$Needle_HS_Control_Control_F <- rowMeans(cbind(countsSplicing$countsEvents$Needle_HS_Control_Control_F, countsSplicing$countsEvents$Needle_HS_Control_Control.1_F))
countsSplicing.jic$countsEvents$Needle_HS_Control_Control_R <- rowMeans(cbind(countsSplicing$countsEvents$Needle_HS_Control_Control_R, countsSplicing$countsEvents$Needle_HS_Control_Control.1_R))
countsSplicing.jic$countsEvents$Needle_HS_Control_Control.1_F <- rowMeans(cbind(countsSplicing$countsEvents$Needle_HS_Control_Control.2_F, countsSplicing$countsEvents$Needle_HS_Control_Control.3_F))
countsSplicing.jic$countsEvents$Needle_HS_Control_Control.1_R <- rowMeans(cbind(countsSplicing$countsEvents$Needle_HS_Control_Control.2_R, countsSplicing$countsEvents$Needle_HS_Control_Control.3_R))

countsSplicing.jic$countsEvents$Needle_HS_Stress_1_F <- rowMeans(cbind(countsSplicing$countsEvents$Needle_HS_Stress_1_F, countsSplicing$countsEvents$Needle_HS_Stress_1.1_F))
countsSplicing.jic$countsEvents$Needle_HS_Stress_1_R <- rowMeans(cbind(countsSplicing$countsEvents$Needle_HS_Stress_1_R, countsSplicing$countsEvents$Needle_HS_Stress_1.1_R))
countsSplicing.jic$countsEvents$Needle_HS_Stress_1.1_F <- rowMeans(cbind(countsSplicing$countsEvents$Needle_HS_Stress_1.2_F, countsSplicing$countsEvents$Needle_HS_Stress_1.3_F))
countsSplicing.jic$countsEvents$Needle_HS_Stress_1.1_R <- rowMeans(cbind(countsSplicing$countsEvents$Needle_HS_Stress_1.2_R, countsSplicing$countsEvents$Needle_HS_Stress_1.3_R))

countsSplicing.jic$countsEvents$Needle_HS_Stress_2_F <- rowMeans(cbind(countsSplicing$countsEvents$Needle_HS_Stress_2_F, countsSplicing$countsEvents$Needle_HS_Stress_2.1_F))
countsSplicing.jic$countsEvents$Needle_HS_Stress_2_R <- rowMeans(cbind(countsSplicing$countsEvents$Needle_HS_Stress_2_R, countsSplicing$countsEvents$Needle_HS_Stress_2.1_R))
countsSplicing.jic$countsEvents$Needle_HS_Stress_2.1_F <- rowMeans(cbind(countsSplicing$countsEvents$Needle_HS_Stress_2.2_F, countsSplicing$countsEvents$Needle_HS_Stress_2.3_F))
countsSplicing.jic$countsEvents$Needle_HS_Stress_2.1_R <- rowMeans(cbind(countsSplicing$countsEvents$Needle_HS_Stress_2.2_R, countsSplicing$countsEvents$Needle_HS_Stress_2.3_R))

colnames(countsSplicing$countsEvents)
colnames(countsSplicing.jic$countsEvents)
countsSplicing.F <- countsSplicing.jic

countsSplicing.F$countsEvents <- countsSplicing.jic$countsEvents[,c(1:36, 205:208, 37:96, 209:216, 97:204)]
colnames(countsSplicing.F$countsEvents)


countsSplicing.jic$psiInfo <- countsSplicing.jic$psiInfo[,-grep("HS", colnames(countsSplicing.jic$psiInfo))] # 37 38 39 40 41 42 43 44 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120
countsSplicing.jic$psiInfo$Needle_HS_Control_Control_F <- rowMeans(cbind(countsSplicing$psiInfo$Needle_HS_Control_Control_F, countsSplicing$psiInfo$Needle_HS_Control_Control.1_F))
countsSplicing.jic$psiInfo$Needle_HS_Control_Control_R <- rowMeans(cbind(countsSplicing$psiInfo$Needle_HS_Control_Control_R, countsSplicing$psiInfo$Needle_HS_Control_Control.1_R))
countsSplicing.jic$psiInfo$Needle_HS_Control_Control.1_F <- rowMeans(cbind(countsSplicing$psiInfo$Needle_HS_Control_Control.2_F, countsSplicing$psiInfo$Needle_HS_Control_Control.3_F))
countsSplicing.jic$psiInfo$Needle_HS_Control_Control.1_R <- rowMeans(cbind(countsSplicing$psiInfo$Needle_HS_Control_Control.2_R, countsSplicing$psiInfo$Needle_HS_Control_Control.3_R))

countsSplicing.jic$psiInfo$Needle_HS_Stress_1_F <- rowMeans(cbind(countsSplicing$psiInfo$Needle_HS_Stress_1_F, countsSplicing$psiInfo$Needle_HS_Stress_1.1_F))
countsSplicing.jic$psiInfo$Needle_HS_Stress_1_R <- rowMeans(cbind(countsSplicing$psiInfo$Needle_HS_Stress_1_R, countsSplicing$psiInfo$Needle_HS_Stress_1.1_R))
countsSplicing.jic$psiInfo$Needle_HS_Stress_1.1_F <- rowMeans(cbind(countsSplicing$psiInfo$Needle_HS_Stress_1.2_F, countsSplicing$psiInfo$Needle_HS_Stress_1.3_F))
countsSplicing.jic$psiInfo$Needle_HS_Stress_1.1_R <- rowMeans(cbind(countsSplicing$psiInfo$Needle_HS_Stress_1.2_R, countsSplicing$psiInfo$Needle_HS_Stress_1.3_R))

countsSplicing.jic$psiInfo$Needle_HS_Stress_2_F <- rowMeans(cbind(countsSplicing$psiInfo$Needle_HS_Stress_2_F, countsSplicing$psiInfo$Needle_HS_Stress_2.1_F))
countsSplicing.jic$psiInfo$Needle_HS_Stress_2_R <- rowMeans(cbind(countsSplicing$psiInfo$Needle_HS_Stress_2_R, countsSplicing$psiInfo$Needle_HS_Stress_2.1_R))
countsSplicing.jic$psiInfo$Needle_HS_Stress_2.1_F <- rowMeans(cbind(countsSplicing$psiInfo$Needle_HS_Stress_2.2_F, countsSplicing$psiInfo$Needle_HS_Stress_2.3_F))
countsSplicing.jic$psiInfo$Needle_HS_Stress_2.1_R <- rowMeans(cbind(countsSplicing$psiInfo$Needle_HS_Stress_2.2_R, countsSplicing$psiInfo$Needle_HS_Stress_2.3_R))

colnames(countsSplicing$psiInfo)
colnames(countsSplicing.jic$psiInfo)

countsSplicing.F$psiInfo <- countsSplicing.jic$psiInfo[,c(1:35, 204:207, 36:95, 208:215, 96:203)]
colnames(countsSplicing.F$psiInfo)

# mean PE

countsSplicing.FPE <- countsSplicing.F
sample_metada <- read.delim("clipboard", header = T)
countsSplicing.FPE$countsEvents <- as.data.frame(matrix(data = NA, nrow = nrow(countsSplicing.F$countsEvents), ncol = length(sample_metada$NameID)+2))
countsSplicing.FPE$psiInfo <- as.data.frame(matrix(data = NA, nrow = nrow(countsSplicing.F$countsEvents), ncol = length(sample_metada$NameID)+2))

lapply(countsSplicing.FPE, dim)

colnames(countsSplicing.FPE$countsEvents) <- c("events.names", "events.length", sample_metada$NameID)
colnames(countsSplicing.FPE$psiInfo) <- c("events.names", sample_metada$NameID)
countsSplicing.FPE$countsEvents$events.names <- countsSplicing.F$countsEvents$events.names
countsSplicing.FPE$countsEvents$events.length <- countsSplicing.F$countsEvents$events.length
countsSplicing.FPE$psiInfo$events.names <- countsSplicing.F$countsEvents$events.names

for(i in 3:ncol(countsSplicing.FPE$countsEvents)){
   NAME <- colnames(countsSplicing.FPE$countsEvents)[i]
   Tech <- sample_metada$Tech[which(sample_metada$NameID == NAME)]
   if(Tech == "PE"){
   Forward <- paste0(NAME, "_F")
   Reverse <- paste0(NAME, "_R")
   countsSplicing.FPE$countsEvents[i] <- round(rowMeans(countsSplicing.F$countsEvents[,c(Forward, Reverse)]), digits = 0)
   countsSplicing.FPE$psiInfo[i-1] <- round(rowMeans(countsSplicing.F$psiInfo[,c(Forward, Reverse)]), digits = 0)
   }else if(Tech == "SE"){
    countsSplicing.FPE$countsEvents[i] <- countsSplicing.F$countsEvents[,NAME]
    countsSplicing.FPE$psiInfo[i-1] <- countsSplicing.F$psiInfo[,NAME]
   }
}
countsSplicing.FPE$psiInfo <- countsSplicing.FPE$psiInfo[,-ncol(countsSplicing.FPE$psiInfo)] 

## Splicing strict filtering 
# >= 10 counts and at least >= 4 counts in >= 2 samples for each isoform
# 1) Filter by expression

keep <- (rowSums(countsSplicing.FPE$countsEvents[3:ncol(countsSplicing.FPE$countsEvents)])  >= 10) & (rowSums(countsSplicing.FPE$countsEvents[3:ncol(countsSplicing.FPE$countsEvents)] >= 4) >= 2)
Events <- names(which((table(countsSplicing.FPE$countsEvents$events.names[keep]) == 2) == TRUE))

countsSplicing.FPE$countsEvents <- countsSplicing.FPE$countsEvents[countsSplicing.FPE$countsEvents$events.names %in% Events,]
countsSplicing.FPE$psiInfo <- countsSplicing.FPE$psiInfo[countsSplicing.FPE$psiInfo$events.names %in% Events,]
write.table(countsSplicing.FPE$countsEvents, "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_NotCorrected/countsEvents_PEandTE.txt", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(countsSplicing.FPE$psiInfo, "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithoutGroup/psiInfo_PEandTE.txt", sep = "\t", quote = F, row.names = T, col.names = T)

# Remove batch effect PE and SE / study from raw counts

sample_metada <- read.delim("clipboard", header = T)

counts <- countsSplicing.FPE$countsEvents[3:ncol( countsSplicing.FPE$countsEvents)]
rownames(counts) <- make.unique(countsSplicing.FPE$countsEvents$events.names) # large and .1 small
test <- ComBat_seq(counts = as.matrix(counts),batch = as.numeric(as.factor(paste0(sample_metada$Batch1, "_", sample_metada$Batch2))), group = paste0(sample_metada$Tissue, "_" ,sample_metada$Treatment))
test <- ComBat_seq(counts = as.matrix(counts),batch = as.numeric(as.factor(paste0(sample_metada$Batch1, "_", sample_metada$Batch2))))
countsSplicing.FPE$countsEvents[3:ncol( countsSplicing.FPE$countsEvents)] <- test
write.table(countsSplicing.FPE$countsEvents, "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithoutGroup/countsEvents_PEandTE_batchCorrectedTechStudy.txt", sep = "\t", quote = F, row.names = T, col.names = T)
write.table(countsSplicing.FPE$psiInfo, "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithoutGroup/psiInfo_PEandTE_batchCorrectedTechStudy.txt", sep = "\t", quote = F, row.names = T, col.names = T)

# read normalization

matrix <- countsSplicing.FPE$countsEvents[,c(1, 3:ncol( countsSplicing.FPE$countsEvents))]
matrix$events.names <- make.unique(matrix$events.names)
rownames(matrix) <- matrix$events.names
matrix <- matrix[,-1]
matrix <- as.matrix(matrix)
dds2 <- DESeqDataSetFromMatrix(matrix, sample_metada, design = ~ Group)
dds2 <- estimateSizeFactors(dds2) 
matrix.counts <- counts(dds2, normalized=TRUE)
countsSplicing.FPE$countsEvents[3:ncol(countsSplicing.FPE$countsEvents)] <- matrix.counts
write.table(countsSplicing.FPE$countsEvents, "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_NotCorrected/countsEvents_lowFilter_PEandTE_Normalized.txt", sep = "\t", quote = F, row.names = T, col.names = T)

# PSI computation

PSI.df <- as.data.frame(matrix(data = NA, nrow = length(unique(countsSplicing.FPE$countsEvents$events.names)), ncol = nrow(sample_metada)+1))
colnames(PSI.df) <- colnames(countsSplicing.FPE$countsEvents)[-2]
PSI.df$events.names <- unique(countsSplicing.FPE$countsEvents$events.names)
rows <- 1:nrow(countsSplicing.FPE$countsEvents)
pairs <- rows %% 2
PSI.df[2:ncol(PSI.df)] <- apply(countsSplicing.FPE$countsEvents[3:ncol(countsSplicing.FPE$countsEvents)], MARGIN = 2, function(x){
   PSI <- x[pairs == 1]/(x[pairs == 1] + x[pairs == 0])
   PSI
})
write.table(PSI.df, "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_NotCorrected/PSI_lowFilter_PEandTE_Normalized.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# FINAL EVENTS: PSI between 0.1-0.9 (spliced) in at least 10% of samples and/or PSI range (max-min) >= 0.25
# potential repeat for unused processing sets
PSI.df <- read.table("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/0.PSI_lowFilter_PEandTE_batchCorrectedTechStudy_Normalized.txt", header = T)

retain.list <- list()

for(i in 1:nrow(PSI.df)){
   row.info <- as.vector(PSI.df[i, 2:ncol(PSI.df)])
   row.info <- row.info[!is.na(row.info)]
   range <- max(unlist(row.info)) - min(unlist(row.info))
   presence <- length(which(unlist(row.info) >= 0.1 & unlist(row.info) <= 0.9)) >= 0.1*(ncol(PSI.df)-1)
   retain.list[[i]] <- FALSE
   if(range >= 0.25 | presence){
      retain.list[[i]] <- TRUE
   }
}

Final.EventIDs <- PSI.df$events.names[unlist(retain.list)]
write.table(as.data.frame(Final.EventIDs),"G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/Final_EventIDs.txt", quote = F, col.names = T, row.names = F, sep = "\t")
Final_PSI.df <- PSI.df[PSI.df$events.names %in% Final.EventIDs,]
write.table(Final_PSI.df,"G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/Final_PSI_lowFilter_PEandTE_Normalized_PresenceRange.txt", quote = F, col.names = T, row.names = F, sep = "\t")
countsSplicing.FPE$countsEvents <- read.table("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/0.countsEvents_lowFilter_PEandTE_batchCorrectedTechStudy_Normalized.txt", header = T)
countsSplicing.FPE$countsEvents <- countsSplicing.FPE$countsEvents[countsSplicing.FPE$countsEvents$events.names %in% Final.EventIDs,]
write.table(countsSplicing.FPE$countsEvents,"G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/Final_countsEvents_lowFilter_PEandTE_Normalized_PresenceRange.txt", quote = F, col.names = T, row.names = F, sep = "\t")
countsSplicing.FPE$psiInfo <- read.table("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/0.psiInfo_lowFilter_PEandTE_batchCorrectedTechStudy.txt", header = T)
countsSplicing.FPE$psiInfo <- countsSplicing.FPE$psiInfo[,-ncol(countsSplicing.FPE$psiInfo)]
countsSplicing.FPE$psiInfo <- countsSplicing.FPE$psiInfo[countsSplicing.FPE$psiInfo$events.names %in% Final.EventIDs,]
write.table(countsSplicing.FPE$psiInfo,"G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/Final_psiInfo_PEandTE_batchCorrectedTechniqueStudy_PresenceRange.txt", quote = F, col.names = T, row.names = F, sep = "\t")
Rawcounts <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/0.countsEvents_lowFilter_PEandTE_batchCorrectedTechStudy.txt", header = T, sep = "\t")
Rawcounts <- Rawcounts[Rawcounts$events.names %in% Final.EventIDs,]
write.table(Rawcounts,"G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/Final_countsEvents_lowFilter_PEandTE_batchCorrectedTechStudy_PresenceRange.txt", quote = F, col.names = T, row.names = F, sep = "\t")

#### 2. Expand annotation to the nearest gene in a range | Generate annotation database per taeda gene and splicing event ####

## Gene-level annotation for expression

Gene.annotation <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/5.Pinus_taeda/Annotation/Pta_Annotation_F.txt", header = T)

## Event-level annotation for splicing

# EVENTID; START, END, SCAFF, STRAND, GeneID, EVENTLENGTH, Frameshift, CDS, Type, Group (NAs till filled)

Event.annotation <- data.frame(EVENTID = Final.EventIDs, START = NA, END = NA, SCAFF = NA, STRAND = NA, SEQ = NA, GENEID = NA,
                               EVENTLENGTH = NA, Frameshift = NA, CDS = NA, PTC = NA, Type = NA, CoreSet = NA, GROUP = NA)

k2rg <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/5.Pinus_taeda/Splicing/Pra-GE-ATLAS-AS-all_type1_kiss2genome.tsv", header = T, sep = "\t")
sequences <- seqinr::read.fasta("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/5.Pinus_taeda/Splicing/Pra-GE-ATLAS-AS-all_type1.fa", whole.header = F, as.string = T, forceDNAtolower = F)
names.tmp <- do.call(rbind, strsplit(names(sequences), split = "|", fixed = T))[,c(1,2)]
names.tmp <- as.data.frame(names.tmp)
names.tmp$Final <- NA
names.tmp$Final[(1:nrow(names.tmp) %% 2) == 1] <- paste0(names.tmp$V1[(1:nrow(names.tmp) %% 2) == 1], "|", names.tmp$V2[(1:nrow(names.tmp) %% 2) == 1], ".U")
names.tmp$Final[(1:nrow(names.tmp) %% 2) == 0] <- paste0(names.tmp$V1[(1:nrow(names.tmp) %% 2) == 0], "|", names.tmp$V2[(1:nrow(names.tmp) %% 2) == 0], ".L")
names(sequences) <- names.tmp$Final

Svariable <- data.frame(kissID = names(sequences), Sequences = unlist(sequences), Lengths = unlist(lapply(sequences, nchar)))

# helper functions

compare.DNA <- function(x,y){
   ########################
   # compare.DNA: search differences between two sequences
   ########################
   # x - sequence of nts
   # y - sequence of nts
   ########################
   sapply(seq(length(x)),
          function(i){
             x[i] == y[i]
          }
   )
}
Variants.seq.diff <- function(matrix1, kmer, export = TRUE, path) {
   ########################
   # Variants.seq.diff: KisSplice extension. Finds the variable sequence (S) between two isoforms. 
   ########################
   # matrix1 - matrix or dataframe (n events, ), first column kissIDs, second column Sequences and third column lengths.
   # kmer - kmer used in KisSplice
   # export - logical. Export S variable seq as fasta to indicated path
   ########################
   Lower <- matrix1[(1:nrow(matrix1) %% 2) == 0,]
   Upper <- matrix1[(1:nrow(matrix1) %% 2) == 1,]
   matrix1 <- data.frame(kissID = Upper[,1], Upper.Seq = Upper[,2], Lower.seq = Lower[,2], Length = Upper[,3] - Lower[,3])
   lower <- as.matrix(matrix1[,3])
   upper <- as.matrix(matrix1[,2])
   list.lower <- list()
   list.upper <- list()
   list.length <- list()
   for (i in seq_along(lower)) {
      seq_lower <- unlist(strsplit(lower[i,], split = ""))
      seq_upper <- unlist(strsplit(upper[i,], split = ""))
      list.lower[[i]] <- seq_lower
      list.upper[[i]] <- seq_upper
      list.length[[i]] <- length(seq_upper) - length(seq_lower) 
   }
   
   list.Svariable <- list()
   for (k in 1:length(list.lower)) {
      TrueFalse <- compare.DNA(list.upper[[k]], list.lower[[k]])
      S0 <- grep("FALSE", TrueFalse)[1]
      #if(S0 != (kmer + 1)){print(stop("\n Error01: wrong kmer"))}
      Svariable <- list.upper[[k]][S0: (S0+ list.length[[k]] -1)]
      list.Svariable[[k]] <- paste(Svariable,collapse = "")
   }
   df <- data.frame(KissID = matrix1[,1],
                    Lower_seq = lower[,1], 
                    Upper_seq = upper[,1],
                    Length_Variable = unlist(list.length, use.names = FALSE),
                    Variable_seq = unlist(list.Svariable, use.names = FALSE))
   if(export){
      dataframe2fas(df[,c(1,5)], file = paste0(path, "allSparts.fasta"))
   }
   return(df)
}

# extract variable sequences

Variant.seqs.df <- Variants.seq.diff(matrix1 = Svariable, kmer = 51, export = F)
# quality check
which(((Svariable$Lengths[(1:nrow(Svariable) %% 2) == 1] - Svariable$Lengths[(1:nrow(Svariable) %% 2) == 0]) == nchar(Variant.seqs.df$Variable_seq)) == FALSE)
Variant.seqs.df$KissID <- do.call(rbind, strsplit(Variant.seqs.df$KissID, split = ".", fixed = T))[,1] 

# Main annotation loop

which((Event.annotation$EVENTID %in% k2rg$X16.Event_name) == FALSE)
StopCodons <- c("TAG", "TAA", "TGA")

for(i in 1:nrow(Event.annotation)){
   
   hit <- which(k2rg$X16.Event_name == Event.annotation$EVENTID[i])
   
   POS <- strsplit(x = k2rg$X3.Chromosome_and_genomic_position[hit], split = ":", fixed = T)[[1]]
   
   # Start
   
   Event.annotation$START[i] <- strsplit(POS[2], "-", fixed = T)[[1]][1]
   
   # End
   
   Event.annotation$END[i] <- strsplit(POS[2], "-", fixed = T)[[1]][2]
   
   # Scaff
   
   Event.annotation$SCAFF[i] <- POS[1]
   
   # Strand
   
   Event.annotation$STRAND[i] <- k2rg$X4.Strand[hit] 
   
   # GeneID
   
   Event.annotation$GENEID[i] <- k2rg$X.1.Gene_Id[hit] 
   
   # Event length
   
   Event.annotation$EVENTLENGTH[i] <- k2rg$X6.Variable_part_length[hit]
   
   # Frameshift
   
   Event.annotation$Frameshift[i] <- k2rg$X7.Frameshift_.[hit]
   
   # CDS
   
   Event.annotation$CDS[i] <- k2rg$X8.CDS_.[hit]
   
   # EVENT TYPE
   
   Event.annotation$Type[i] <- k2rg$X5.Event_type[hit]
   
   # Variable sequence
   
   hit2 <- which(Variant.seqs.df$KissID == Event.annotation$EVENTID[i])
   Event.annotation$SEQ[i] <- Variant.seqs.df$Variable_seq[hit2]
   
   # PTC
   
   PTC.list <- list()
   for(j in 1:length(StopCodons)){
      PTC.list[[j]] <- length(grep(StopCodons[j], Event.annotation$SEQ[i]))
   }
   res <- sum(unlist(PTC.list))
   if(res > 0){
      Event.annotation$PTC[i] <- "True"   
   }else if(res == 0){
      Event.annotation$PTC[i] <- "False"
   }
   
}

## Expand to closest gene in a 1500 window | re locate all CDS events just in case by coordinates in a range (200-300 nt) and check proportion

# expand to closest genes
TxDb.Pta <- GenomicFeatures::makeTxDbFromGFF(file = "E:/ATLAS_Pra/Pita.2_01_PraGEATLASF.gtf", format = "gtf")
QueryEvents <- Event.annotation[is.na(Event.annotation$GENEID) & !is.na(Event.annotation$START), c(4,2,3,5,1, 12)]
colnames(QueryEvents) <- c("chr", "start", "end","strand", "EventID", "Type")
QueryRanges <- makeGRangesFromDataFrame(QueryEvents, keep.extra.columns = T)
genes <- genes(TxDb.Pta)
Hits <- distanceToNearest(QueryRanges, genes,  ignore.strand=FALSE)
rownames(Event.annotation) <- Event.annotation$EVENTID
rowids <- QueryRanges@elementMetadata$EventID[Hits[Hits@elementMetadata$distance <= 3000]@from]
newids <- genes$gene_id[Hits[Hits@elementMetadata$distance <= 3000]@to]
Event.annotation[rowids,7] <- newids

# CDS hit
table(Event.annotation$CDS)
QueryEvents <- Event.annotation[!is.na(Event.annotation$START), c(4,2,3,5,1, 12)]
colnames(QueryEvents) <- c("chr", "start", "end","strand", "EventID", "Type")
QueryEvents$strand <- gsub("+,", "", QueryEvents$strand, fixed = T)
QueryEvents$strand <- gsub("-,", "", QueryEvents$strand, fixed = T)
QueryRanges <- makeGRangesFromDataFrame(QueryEvents, keep.extra.columns = T)
cdss <- cds(TxDb.Pta)
Hits <- distanceToNearest(QueryRanges, cdss,  ignore.strand=FALSE)
rownames(Event.annotation) <- Event.annotation$EVENTID
rowids <- QueryRanges@elementMetadata$EventID[Hits[Hits@elementMetadata$distance == 0]@from]
Event.annotation[rowids,10] <- "True"
table(Event.annotation$CDS)
which(!is.na(Event.annotation$CDS) & is.na(Event.annotation$GENEID)) # check


# parse major types (AS events classified without geneID; genes not annotated)

Event.annotation.Stric <- Event.annotation
table(Event.annotation.Stric$Type)
Event.annotation.Stric$Type[which(Event.annotation.Stric$Type == "altA")] <- "AltAD"
Event.annotation.Stric$Type[which(Event.annotation.Stric$Type == "altAD")] <- "AltAD"
Event.annotation.Stric$Type[which(Event.annotation.Stric$Type == "altD")] <- "AltAD"
Event.annotation.Stric$Type[which(Event.annotation.Stric$Type == "-")] <- "AS"
Event.annotation.Stric$Type[which(Event.annotation.Stric$Type == "insertion")] <- "indel"
Event.annotation.Stric$Type[which(Event.annotation.Stric$Type == "deletion")] <- "indel"
Event.annotation.Stric$Type[which(Event.annotation.Stric$Type == "ES_altA")] <- "ES"
Event.annotation.Stric$Type[which(Event.annotation.Stric$Type == "ES_altAD")] <- "ES"
Event.annotation.Stric$Type[which(Event.annotation.Stric$Type == "ES_altD")] <- "ES"
Event.annotation.Stric$Type[which(Event.annotation.Stric$Type == "ES_MULTI")] <- "ES"
Event.annotation.Stric$Type[which(Event.annotation.Stric$Type == "ES_MULTI_altA")] <- "ES"
Event.annotation.Stric$Type[which(Event.annotation.Stric$Type == "ES_MULTI_altAD")] <- "ES"
Event.annotation.Stric$Type[which(Event.annotation.Stric$Type == "ES_MULTI_altD")] <- "ES"
Event.annotation <- Event.annotation.Stric
table(Event.annotation.Stric$Type)

# each row for each event and gene (separate_longer_delim like multiple columns )

Event.annotation.split <- separate_longer_delim(Event.annotation, cols = c(STRAND, GENEID), delim = ",") # 6462 genes with structure variation
table(Event.annotation.split$Type)

# export final event annotation table

colnames(Gene.annotation)
colnames(Event.annotation)
Event.annotation.split$Mercator4.Bin <- NA
Event.annotation.split$Mercator4.Desc <- NA
Event.annotation.split$EggNOG.Symbol <- NA
Event.annotation.split$EggNOG.GOm <- NA
Event.annotation.split$EggNOG.GOc <- NA
Event.annotation.split$EggNOG.GOb <- NA
Event.annotation.split$Phylostratum <- NA
Event.annotation.split$GeneFamilyFounder.ID <- NA
Event.annotation.split$GeneFamilyFounder.Age <- NA

for(i in 1:nrow(Event.annotation.split)){
   if(is.na(Event.annotation.split$GENEID[i])){next}
   hit <- which(Gene.annotation$GeneID == Event.annotation.split$GENEID[i])
   Event.annotation.split$Mercator4.Bin[i] <- Gene.annotation$Mercator4.Bin[hit]
   Event.annotation.split$Mercator4.Desc[i] <- Gene.annotation$Mercator4.Desc[hit]
   Event.annotation.split$EggNOG.Symbol[i] <- Gene.annotation$EggNOG.Symbol[hit]
   Event.annotation.split$EggNOG.GOm[i] <- Gene.annotation$EggNOG.GOm[hit]
   Event.annotation.split$EggNOG.GOc[i] <- Gene.annotation$EggNOG.GOc[hit]
   Event.annotation.split$EggNOG.GOb[i] <- Gene.annotation$EggNOG.GOb[hit]
   Event.annotation.split$Phylostratum[i] <- Gene.annotation$PhyloStratum[hit]
   Event.annotation.split$GeneFamilyFounder.ID[i] <- Gene.annotation$GeneFamilyFounder.ID[hit]
   Event.annotation.split$GeneFamilyFounder.Age[i] <- Gene.annotation$GeneFamilyFounder.Age[hit]
}

write.table(Event.annotation.split, "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/Event_annotation.txt", sep = "\t", quote = F, col.names = T, row.names = F)

#### 3. Descriptives ####

#### UMAP (Expression vst; Splicing PSI)
## Use this to select BatchCorrected without group or raw

labelitass <- list(Tissues = sample_metada$Tissue, Group = sample_metada$Group, Stress = sample_metada$Stress, Treatment =  sample_metada$Treatment, Intensity = sample_metada$Time, Genotype = sample_metada$Genotype, Technique = sample_metada$Batch1, Study = sample_metada$Batch2,
                   Treatment.Intensity = paste0(sample_metada$Treatment, "_", sample_metada$Time), Tissue.Stress = paste0(sample_metada$Tissue, "_", sample_metada$Stress))
labelitass$Intensity <- as.character(labelitass$Intensity)
labelitass$Study <- as.character(labelitass$Study)
ncolors <- lapply(labelitass, function(x) length(levels(as.factor(x))))
fcolors <- c(RColorBrewer::brewer.pal(12, "Paired") ,"black",  "peachpuff3", "lightgoldenrod")
color.labelitass <- list(Tissues = c("purple", "darkgreen", "salmon3", "gold", "lightsteelblue3"),
                         Group = c("purple", "darkgreen", "salmon3", "gold", "lightsteelblue3"),
                         Stress = c("gold", "purple"), 
                         Treament = fcolors[1:6],
                         Intensity = c("gray68", "gray48", "gray28", "gray18", "gray8"),
                         Genotype = c("purple", "black", "gold"),
                         Technique = c("purple", "gold"),
                         Study = fcolors[1:7],
                         Treatment.Intensity = fcolors[1:14],
                         Tissue.Stress = fcolors[1:7])

# Expression

set.seed(404)
rawinputs <- list(NoBatch = t(read.table("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Batch_NotCorrected/Final_VST_lowFilt_Presence.txt", header = T)),
                  Batch = t(read.table("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Batch_Corrected/Final_VST_lowFilt_Presence_batchCorrected.txt", header = T)))
inputs <- lapply(t(rawinputs), function(x) umap(x, n_components = 5, n_threads = 7))
inputs.pca <- lapply(t(rawinputs), function(x) umap(x, n_components = 5, n_threads = 7, pca = 50))
inputs <- lapply(inputs, function(x){ 
   colnames(x) <- paste0("UMAP", 1:ncol(x))
   x
})

inputs.pca <- lapply(inputs.pca, function(x){ 
   colnames(x) <- paste0("UMAP", 1:ncol(x))
   x
})
INPUTS.F <- list(WI = inputs.pca, WO = inputs)
write.table(INPUTS.F$WO[[1]], file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/UMAP_Expression_NoBatch.txt", sep = "\t", col.names = T, row.names = T)

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

# Splicing

set.seed(404)
rawinputs <- list(NoBatch = read.table("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_NotCorrected/Final_PSI_lowFilter_PEandTE_Normalized_PresenceRange.txt", header = T),
                  Batch = read.table("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithoutGroup/Final_PSI_lowFilter_PEandTE_batchCorrectedTechStudy_Normalized_PresenceRange.txt", header = T))
rawinputs <- lapply(rawinputs, function(x){
   rownames(x) <- x[,1]
   x <- t(x[,-1])
   x[is.na(x)] <- 0.5
   x
})

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
write.table(INPUTS.F$WI[[1]], file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/UMAP/Splicing/UMAP_Splicing_NoBatch_PCA.txt", sep = "\t", col.names = T, row.names = T)

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


#### 3.0 Define sets ####
# DiffSplicing with batch corrected raw counts and DESeq2 with raw counts and including batch in design

### Splicing ###

## PanAs + ASNR
# Presence at least in 20%/total samples of samples with sufficient coverage (not NaN)
# 0.1 < PSI < 0.9 in 70% of the samples with sufficient read coverage (may be 60%)

PanAS.list <- list()
for(i in 1:nrow(Final_PSI.df)){
   PanAS.list[[i]] <- FALSE
   row.info <- as.vector(Final_PSI.df[i, 2:ncol(Final_PSI.df)])
   row.info <- row.info[!is.na(row.info)]
   # Presence in 20%  of samples with sufficient read coverage (not NaN)
   presence <- length(row.info) >= 0.2*(ncol(Final_PSI.df)-1)
   # 0.1 < PSI < 0.9 in at least 70% of the samples with sufficient read coverage
   splicing <- length(which(unlist(row.info) >= 0.1 & unlist(row.info) <= 0.9)) >= 0.7*(length(row.info))
   if(splicing & presence){
      PanAS.list[[i]] <- TRUE
   }
}

PanAS.events <- Final.EventIDs[unlist(PanAS.list)]
length(PanAS.events)
Final_PSI.df[Final_PSI.df$events.names %in% PanAS.events,]
write.table(Final_PSI.df[Final_PSI.df$events.names %in% PanAS.events,1], file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/PanAS_EventsID.txt", col.names = T, quote = F, row.names = F)

## TissuesAS 
# have at least two replicates without NaN for each tissue evaluated (Needle, Buds, Vascular tissue); Megagametophyte excluded
# Absolute Delta PSI >  25
# Using all the samples under Tissues group + control samples of other groups
## AS-NR Tissues 
# have at least two replicates without NaN for each tissue evaluated (Needle, Buds, Vascular tissue); Megagametophyte excluded
# Absolute Delta PSI < 5
# mean PSI tissue 0.1 < PSI < 0.9 (Ensure that is alternatively spliced)
## Genome Tissues
# Same as above but without PSI filtering

TissuesSamples <- c(1:16, 18:52)
Tissues.Final_PSI.df <- Final_PSI.df[,c(1,TissuesSamples+1)]
colnames(Tissues.Final_PSI.df) <- c("events.names",rep("Bud", 16), rep("Needle", 27), rep("Vascular", 8))

TissuesAS.list <- list()
TissuesASNR.list <- list()
TissuesGenome.list <- list()

for(i in 1:nrow(Tissues.Final_PSI.df)){
   for(j in unique(colnames(Tissues.Final_PSI.df)[-1])){
      
   # initial set
   
   TissuesAS.list[[j]][[i]] <- FALSE
   TissuesASNR.list[[j]][[i]] <- FALSE
   TissuesGenome.list[[j]][[i]] <- FALSE
      
   # extract information
   
   Source <- grep(j, colnames(Tissues.Final_PSI.df))
   Target <- grep(j, colnames(Tissues.Final_PSI.df), invert = T)[-1]
   
   row.info.source <- as.vector(Tissues.Final_PSI.df[i, Source])
   row.info.source <- row.info.source[!is.na(row.info.source)]
   if(length(row.info.source) == 0){ next }
   
   row.info.target <- as.vector(Tissues.Final_PSI.df[i, Target])
   row.info.target <- row.info.target[!is.na(row.info.target)]
   if(length(row.info.target) == 0){ next }
   
   # presence: at least two replicates with enough coverage in the three group of tissues
   
   presence <- length(which((table(do.call(rbind, strsplit(names(c(unlist(row.info.source), unlist(row.info.target))), split = ".", fixed = T))[,1]) >= 2) == TRUE)) == 3
   
   if(!presence){
      next
   }else if(presence){
      TissuesGenome.list[[j]][[i]] <- TRUE   
      }
   
   
   
   # mean PSI and difference and absolute delta
   
   AbsPSI <- abs(mean(unlist(row.info.source)) - mean(unlist(row.info.target)))
   
   # mean PSI per tissue for TissueASNR other condition set
   
   MeanPSIperTissue <- c(unlist(row.info.source), unlist(row.info.target))
   names(MeanPSIperTissue) <- do.call(rbind, strsplit(names(MeanPSIperTissue), split = ".", fixed = T))[,1]
   agg <- aggregate(MeanPSIperTissue, by= list(names = names(MeanPSIperTissue)), FUN=mean)
   splicing <- length(which(agg$x >= 0.1 & agg$x <= 0.9)) > 2
   
   # Subset AS events
   
   if(AbsPSI >= 0.25){
      TissuesAS.list[[j]][[i]] <- TRUE
   }else if(AbsPSI <= 0.05 & splicing){
      TissuesASNR.list[[j]][[i]] <- TRUE
   }
   
   
   } # tissues loop
} # row loop

TissuesAS.events <- unique(c(Tissues.Final_PSI.df$events.names[unlist(TissuesAS.list$Bud)],Tissues.Final_PSI.df$events.names[unlist(TissuesAS.list$Needle)], Tissues.Final_PSI.df$events.names[unlist(TissuesAS.list$Vascular)]))
TissuesASNR.events <- unique(c(Tissues.Final_PSI.df$events.names[unlist(TissuesASNR.list$Bud)],Tissues.Final_PSI.df$events.names[unlist(TissuesASNR.list$Needle)], Tissues.Final_PSI.df$events.names[unlist(TissuesASNR.list$Vascular)]))
TissuesGenome.events <- unique(c(Tissues.Final_PSI.df$events.names[unlist(TissuesGenome.list$Bud)],Tissues.Final_PSI.df$events.names[unlist(TissuesGenome.list$Needle)], Tissues.Final_PSI.df$events.names[unlist(TissuesGenome.list$Vascular)]))

length(TissuesAS.events)
length(TissuesASNR.events)
length(TissuesGenome.events)

Final_PSI.df[Final_PSI.df$events.names %in% TissuesAS.events,]
Final_PSI.df[Final_PSI.df$events.names %in% TissuesASNR.events,]
write.table(Final_PSI.df[Final_PSI.df$events.names %in% TissuesAS.events,1], file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/TissuesAS_EventsID.txt", col.names = T, quote = F, row.names = F)
write.table(Final_PSI.df[Final_PSI.df$events.names %in% TissuesASNR.events,1], file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/TissuesASNR_EventsID.txt", col.names = T, quote = F, row.names = F)
write.table(Final_PSI.df[Final_PSI.df$events.names %in% TissuesGenome.events,1], file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/TissuesGenome_EventsID.txt", col.names = T, quote = F, row.names = F)

## StressAS 
# have at least two replicates without NaN for each control and stress matched samples (PH, DO, HS, FU, FU)
# Absolute Delta PSI > 15 in the same direction in at least two of the five stress experiments (pos and neg Delta PSI)
## AS-NR Stress 
# have at least two replicates without NaN for each control and stress matched samples (PH, DO, HS, FU, FU)
# Absolute Delta PSI < 5
# at least one sample with PSI between 0.1 < PSI < 0.9 (Ensure that is alternatively spliced)
## Genome Stress
# Same as above but without PSI filtering

StressSamples <- c(19,20,22,23,25:36, 54:142)
StressColumns <- c(1,StressSamples)
Stress.Final_PSI.df <- Final_PSI.df[,StressColumns]
colnames(Stress.Final_PSI.df) <- c("events.names",
                                   rep("DO_Control", 2),
                                   rep("HS_Control", 2),
                                   rep("PH_Control", 12),
                                   rep("HS_Stress", 4),
                                   rep("DO_Stress", 10),
                                   rep("PH_Stress", 6),
                                   rep("FU1_Control", 3),
                                   rep("FU2_Control", 18),
                                   rep("FU1_Stress", 4),
                                   rep("FU2_Stress", 44)
                                   )# study, stress/control

StressAS.list <- list()
StressASNR.list <- list()
StressGenome.list <- list()

for(i in 1:nrow(Stress.Final_PSI.df)){
   for(j in c("DO", "HS", "PH", "FU1", "FU2")){
      
      # initial set
      
      StressAS.list[[j]][[i]] <- 0
      StressASNR.list[[j]][[i]] <- FALSE
      StressGenome.list[[j]][[i]] <- FALSE
      
      # extract information
      
      Source <- grep(paste0(j, "_Stress"), colnames(Stress.Final_PSI.df))
      Target <- grep(paste0(j, "_Control"), colnames(Stress.Final_PSI.df))
      
      row.info.source <- as.vector(Stress.Final_PSI.df[i, Source])
      row.info.source <- row.info.source[!is.na(row.info.source)]
      if(length(row.info.source) == 0){ 
         StressAS.list[[j]][[i]] <- 0
         next 
         }
      
      row.info.target <- as.vector(Stress.Final_PSI.df[i, Target])
      row.info.target <- row.info.target[!is.na(row.info.target)]
      if(length(row.info.target) == 0){ 
         StressAS.list[[j]][[i]] <- 0
         next
         }
      
      # presence: at least two replicates with enough coverage in each group studied
      
      presence <- length(which((table(do.call(rbind, strsplit(names(c(unlist(row.info.source), unlist(row.info.target))), split = ".", fixed = T))[,1]) >= 2) == TRUE)) == 2
      
      if(!presence){
         StressAS.list[[j]][[i]] <- 0
         next
      }else if(presence){
         StressGenome.list[[j]][[i]] <- TRUE   
      }
      
      # Delta PSI for future check with same direction and |PSI| > 0.15
      
      PSIn <- mean(unlist(row.info.source)) - mean(unlist(row.info.target))
      
      # at least one replicate (control or stress) with 0.1 < PSI < 0.9 (ensure alternative spliced)
      
      PSIper <- c(unlist(row.info.source), unlist(row.info.target))
      splicing <- length(which(PSIper >= 0.1 & PSIper <= 0.9)) >= 1
      
      # Subset AS events
      
      StressAS.list[[j]][[i]] <- PSIn # future subset with with |PSI| > 0.15 same direction in at least two experiments
      if(abs(PSIn) <= 0.05 & splicing){
         StressASNR.list[[j]][[i]] <- TRUE
      }
      
      
   } # stress loop
} # row loop

# >= 0.15 absolute and at least two replicates with same sign
StressAS.list.df <- data.frame(DO = unlist(StressAS.list$DO), HS = unlist(StressAS.list$HS),
                               PH = unlist(StressAS.list$PH), FU1 = unlist(StressAS.list$FU1),
                               FU2 = unlist(StressAS.list$FU2))

StressAS.list.df <- t(StressAS.list.df)
StressAS.list2 <- list()
for(i in 1:ncol(StressAS.list.df)){
   val <- StressAS.list.df[,i]
   hits <- which(abs(val) >= 0.15)
   if(length(hits) <= 1){
      StressAS.list2[[i]] <- FALSE
      next
   }else if(length(hits) > 1){
         con <- val[hits]
      if(length(which(con < 0)) >= 2 | length(which(con > 0)) >= 2){
         StressAS.list2[[i]] <- TRUE
      }else {
         StressAS.list2[[i]] <- FALSE
      }
   }
}

StressAS.events <- Stress.Final_PSI.df$events.names[unlist(StressAS.list2)]
StressASNR.events <- unique(c(Stress.Final_PSI.df$events.names[unlist(StressASNR.list$DO)],
                              Stress.Final_PSI.df$events.names[unlist(StressASNR.list$HS)],
                              Stress.Final_PSI.df$events.names[unlist(StressASNR.list$PH)],
                              Stress.Final_PSI.df$events.names[unlist(StressASNR.list$FU1)],
                              Stress.Final_PSI.df$events.names[unlist(StressASNR.list$FU2)]))
StressGenome.events <- unique(c(Stress.Final_PSI.df$events.names[unlist(StressGenome.list$DO)],
                                Stress.Final_PSI.df$events.names[unlist(StressGenome.list$HS)],
                                Stress.Final_PSI.df$events.names[unlist(StressGenome.list$PH)],
                                Stress.Final_PSI.df$events.names[unlist(StressGenome.list$FU1)],
                                Stress.Final_PSI.df$events.names[unlist(StressGenome.list$FU2)]))

length(StressAS.events)
length(StressASNR.events)
length(StressGenome.events)

write.table(Final_PSI.df[Final_PSI.df$events.names %in% StressAS.events,1], file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/StressAS_EventsID.txt", col.names = T, quote = F, row.names = F)
write.table(Final_PSI.df[Final_PSI.df$events.names %in% StressASNR.events,1], file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/StressASNR_EventsID.txt", col.names = T, quote = F, row.names = F)
write.table(Final_PSI.df[Final_PSI.df$events.names %in% StressGenome.events,1], file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/StressGenome_EventsID.txt", col.names = T, quote = F, row.names = F)

## Final ASNR & Genome set
# Tissues and Stress ASNR and GenomeSet intersection

ASNR <- intersect(TissuesASNR.events, StressASNR.events)
Genome <- intersect(TissuesGenome.events, StressGenome.events)

write.table(Final_PSI.df[Final_PSI.df$events.names %in% ASNR,1], file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/ASNR_EventsID.txt", col.names = T, quote = F, row.names = F)
write.table(Final_PSI.df[Final_PSI.df$events.names %in% Genome,1], file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/Genome_EventsID.txt", col.names = T, quote = F, row.names = F)

### Expression ###

# normalized counts for filtering and raw counts for expression testing

Raw.counts.expression <- read.table("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Batch_Corrected/Final_RawCounts_lowFilt_Presence_batchCorrected.txt", header = T) 
Raw.counts.expression <- as.matrix(Raw.counts.expression)
dds3 <- DESeqDataSetFromMatrix(Raw.counts.expression, sample_metada, design = ~ Group)
dds3 <- estimateSizeFactors(dds3) 
matrix.counts.expression <- counts(dds3, normalized=TRUE)

## Genome Universe

write.table(data.frame(GenomeUniverse = rownames(matrix.counts.expression)), "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Batch_Corrected/GenomeUniverse_GeneID.txt", quote = F, sep = "\t", col.names = T, row.names = F)

## Pan Expressed genes (greater than 20 normalized counts in at least 70 % of the samples)

PanGE.genes <- list()
for(i in 1:nrow(matrix.counts.expression)){
   PanGE.genes[[i]] <- length(which((matrix.counts.expression[i,] >= 20) == TRUE)) >= ncol(matrix.counts.expression)*0.7
}

PanGE.genes <- rownames(matrix.counts.expression)[unlist(PanGE.genes)]

write.table(data.frame(PanGE.genes = PanGE.genes), "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Batch_Corrected/PanGE_GeneID.txt", quote = F, sep = "\t", col.names = T, row.names = F)

## Tissues
# median per tissue type normalized count >= 5 in at least one tissue
# logFC >= 3 (only +) as tissue specific

TissuesSamples <- c(1:16, 18:52)
matrix.counts.expression.tissues <- matrix.counts.expression[,TissuesSamples]
Raw.counts.expression.tissues <- Raw.counts.expression[,TissuesSamples]
sample_metada.tissues <- sample_metada[TissuesSamples,]
sample_metada.tissues$Bud <- c(rep("Bud", 16), rep("Other", 27), rep("Other", 8))
sample_metada.tissues$Needle <- c(rep("Other", 16), rep("Needle", 27), rep("Other", 8))
sample_metada.tissues$Vascular <- c(rep("Other", 16), rep("Other", 27), rep("Vascular", 8))
class(sample_metada.tissues$Bud)
sample_metada.tissues$Bud <- as.factor(sample_metada.tissues$Bud)# to factor
class(sample_metada.tissues$Bud)
class(sample_metada.tissues$Needle)
sample_metada.tissues$Needle <- as.factor(sample_metada.tissues$Needle)# to factor
class(sample_metada.tissues$Needle)
class(sample_metada.tissues$Vascular)
sample_metada.tissues$Vascular <- as.factor(sample_metada.tissues$Vascular)# to factor
class(sample_metada.tissues$Vascular)
dds.tissues <- DESeqDataSetFromMatrix(Raw.counts.expression.tissues, sample_metada.tissues, design = ~ Tissue)

TissuesGE.list <- list()
colnames(matrix.counts.expression.tissues) <- gsub("Xylem_", "Vascular_", colnames(matrix.counts.expression.tissues), fixed = T)
colnames(matrix.counts.expression.tissues) <- gsub("Phloem_", "Vascular_", colnames(matrix.counts.expression.tissues), fixed = T)
sample_metada.tissues$Tissue <- gsub("Xylem", "Vascular", sample_metada.tissues$Tissue, fixed = T)
sample_metada.tissues$Tissue <- gsub("Phloem", "Vascular", sample_metada.tissues$Tissue, fixed = T)

for(i in 1:nrow(matrix.counts.expression.tissues)){
   for(j in c("Bud", "Needle", "Vascular")){
      
      # initial set
      
      TissuesGE.list[[j]][[i]] <- TRUE
      
      # extract information
      
      Source <- grep(j, colnames(matrix.counts.expression.tissues))
      Target <- grep(j, colnames(matrix.counts.expression.tissues), invert = T)
      
      row.info.source <- as.vector(matrix.counts.expression.tissues[i, Source])
      names(row.info.source) <- sample_metada.tissues$Tissue[Source]
      row.info.target <- as.vector(matrix.counts.expression.tissues[i, Target])
      names(row.info.target) <- sample_metada.tissues$Tissue[Target]
      
      # median per tissue type normalized count >= 5 in at least one tissue
      
      presence.test <- aggregate(c(row.info.source, row.info.target), by= list(tissue = sample_metada.tissues$Tissue), median)
      presence <- length(which(presence.test$x >= 5)) >= 1
        
      if(!presence){
         TissuesGE.list[[j]][[i]] <- FALSE
         next
      }
      
   } # tissues loop
} # row loop

# Differential expression and at least logFC >= +- 1.584963 (|FC| >= 3) (tissue specific corresponding to the tissue evaluated versus the rest)

TissuesLFC.list <- list()

for(j in c("Bud", "Needle", "Vascular")){
   dds.tissues <- DESeqDataSetFromMatrix(Raw.counts.expression.tissues[unlist(TissuesGE.list[j]),], sample_metada.tissues, design = ~ Tissue)
   design(dds.tissues) <- formula(paste("~ ", j, sep = " "))
   dds.tissues <- DESeq(dds.tissues)
   res <- results(dds.tissues, contrast=c(j, j,"Other"))
   TissuesLFC.list[[j]] <- c(rownames(res)[res$log2FoldChange >= 1.584963 & res$padj <= 0.05], rownames(res)[res$log2FoldChange <= -1.584963  & res$padj <= 0.05])
}

TissuesGE.genes <- unique(c(TissuesLFC.list$Bud, TissuesLFC.list$Needle, TissuesLFC.list$Vascular))

length(TissuesGE.genes)

write.table(data.frame(TissuesGE.genes = TissuesGE.genes), file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Batch_Corrected/TissuesGE_GenesID.txt", col.names = T, quote = F, row.names = F)

## Stress
# normalized counts >= 5 in at least two samples from same stress experiment
# logFC >= 2 same sign as tissue specific

StressSamples <- c(19,20,22,23,25:35, 54:142)-1
matrix.counts.expression.stress <- matrix.counts.expression[,StressSamples]
Raw.counts.expression.stress <- Raw.counts.expression[,StressSamples]
sample_metada.stress <- sample_metada[StressSamples,]
sample_metada.stress$diff <- paste0(sample_metada.stress$Group, "_", sample_metada.stress$Stress)
sample_metada.stress$diff[c(36:38, 57:60)] <- c(rep("FU1_Control", 3), rep("FU1_Stress", 4))
sample_metada.stress$diff <- gsub(pattern = "FU_", replacement = "FU2_", x = sample_metada.stress$diff, fixed = T)
class(sample_metada.stress$diff)
sample_metada.stress$diff <- as.factor(sample_metada.stress$diff)# to factor
class(sample_metada.stress$diff)
dds.stress <- DESeqDataSetFromMatrix(Raw.counts.expression.stress, sample_metada.stress, design = ~ diff)

StressGE.list <- list()
colnames(matrix.counts.expression.stress) <- sample_metada.stress$diff

for(i in 1:nrow(matrix.counts.expression.tissues)){
   for(j in c("DO", "HS", "PH","FU1", "FU2")){
      
      # initial set
      
      StressGE.list[[j]][[i]] <- TRUE
      
      # extract information
      
      Source <- grep(paste0(j, "_Control"), colnames(matrix.counts.expression.stress))
      Target <- grep(paste0(j, "_Stress"), colnames(matrix.counts.expression.stress))
      
      row.info.source <- as.vector(matrix.counts.expression.stress[i, Source])
      names(row.info.source) <- sample_metada.stress$diff[Source]
      row.info.target <- as.vector(matrix.counts.expression.stress[i, Target])
      names(row.info.target) <- sample_metada.stress$diff[Target]
      
      # >= 5 normalized counts in two samples (control or stress)
      
      presence <- length(which((c(row.info.source, row.info.target) >= 5) == TRUE)) >= 2
      
      if(!presence){
         StressGE.list[[j]][[i]] <- FALSE
         next
      }
      
   } # tissues loop
} # row loop

# Differential expression and at least logFC >= +- 1 (|FC| >= 2) (same sign in at least 2/5 experiments)

StressLFC.list <- list()
for(j in c("DO", "HS", "PH", "FU1", "FU2")){
   dds.stress <- DESeqDataSetFromMatrix(Raw.counts.expression.stress[unlist(StressGE.list[j]),], sample_metada.stress, design = ~ diff)
   dds.stress <- DESeq(dds.stress)
   res <- results(dds.stress, contrast=c("diff", paste0(j, "_Stress"), paste0(j, "_Control")))
   StressLFC.list[[j]] <- res
}

StressGE.q <- list()
StressGE.q[["POS"]] <- c(rownames(StressLFC.list$DO)[which(StressLFC.list$DO$log2FoldChange >= 1 & StressLFC.list$DO$padj <= 0.05)],
                           rownames(StressLFC.list$HS)[which(StressLFC.list$HS$log2FoldChange >= 1 & StressLFC.list$HS$padj <= 0.05)],
                           rownames(StressLFC.list$PH)[which(StressLFC.list$PH$log2FoldChange >= 1 & StressLFC.list$PH$padj <= 0.05)],
                           rownames(StressLFC.list$FU1)[which(StressLFC.list$FU1$log2FoldChange >= 1 & StressLFC.list$FU1$padj <= 0.05)],
                           rownames(StressLFC.list$FU2)[which(StressLFC.list$FU2$log2FoldChange >= 1 & StressLFC.list$FU2$padj <= 0.05)]
                           )
StressGE.q[["NEG"]] <- c(rownames(StressLFC.list$DO)[which(StressLFC.list$DO$log2FoldChange <= -1 & StressLFC.list$DO$padj <= 0.05)],
                             rownames(StressLFC.list$HS)[which(StressLFC.list$HS$log2FoldChange <= -1 & StressLFC.list$HS$padj <= 0.05)],
                             rownames(StressLFC.list$PH)[which(StressLFC.list$PH$log2FoldChange <= -1 & StressLFC.list$PH$padj <= 0.05)],
                             rownames(StressLFC.list$FU1)[which(StressLFC.list$FU1$log2FoldChange <= -1 & StressLFC.list$FU1$padj <= 0.05)],
                             rownames(StressLFC.list$FU2)[which(StressLFC.list$FU2$log2FoldChange <= -1 & StressLFC.list$FU2$padj <= 0.05)]
)


StressGE.genes <- unique(c(names(table(StressGE.q$POS))[which(table(StressGE.q$POS) >= 2)],names(table(StressGE.q$NEG))[which(table(StressGE.q$NEG) >= 2)]))

length(StressGE.genes)

write.table(data.frame(StressGE.genes = StressGE.genes), file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Batch_Corrected/StressGE_GenesID.txt", col.names = T, quote = F, row.names = F)


#### 3.1 Event type and proportion of each AS type in each core set ####

## Event type per core set

df <- data.frame(CoreSet = c(rep("StressAS", 4), rep("TissuesAS", 4), rep("PanAS", 4), rep("ASNR", 4), rep("Genome", 4)),
                 Events = c(rep(c("AltAD", "Unk-AS", "ES", "IR"), 5)),
                 Percentage = NA)

table(Event.annotation$Type[Event.annotation$EVENTID %in% StressAS.events])[-4]/sum(table(Event.annotation$Type[Event.annotation$EVENTID %in% StressAS.events])[-4])
table(Event.annotation$Type[Event.annotation$EVENTID %in% TissuesAS.events])[-4]/sum(table(Event.annotation$Type[Event.annotation$EVENTID %in% TissuesAS.events])[-4])
table(Event.annotation$Type[Event.annotation$EVENTID %in% PanAS.events])[-4]/sum(table(Event.annotation$Type[Event.annotation$EVENTID %in% PanAS.events])[-4])
table(Event.annotation$Type[Event.annotation$EVENTID %in% ASNR])[-4]/sum(table(Event.annotation$Type[Event.annotation$EVENTID %in% ASNR])[-4])
table(Event.annotation$Type[Event.annotation$EVENTID %in% Genome])[-4]/sum(table(Event.annotation$Type[Event.annotation$EVENTID %in% Genome])[-4])

df$Percentage <- c(table(Event.annotation$Type[Event.annotation$EVENTID %in% StressAS.events])[-4]/sum(table(Event.annotation$Type[Event.annotation$EVENTID %in% StressAS.events])[-4]),
                   table(Event.annotation$Type[Event.annotation$EVENTID %in% TissuesAS.events])[-4]/sum(table(Event.annotation$Type[Event.annotation$EVENTID %in% TissuesAS.events])[-4]),
                   table(Event.annotation$Type[Event.annotation$EVENTID %in% PanAS.events])[-4]/sum(table(Event.annotation$Type[Event.annotation$EVENTID %in% PanAS.events])[-4]),
                   table(Event.annotation$Type[Event.annotation$EVENTID %in% ASNR])[-4]/sum(table(Event.annotation$Type[Event.annotation$EVENTID %in% ASNR])[-4]),
                   table(Event.annotation$Type[Event.annotation$EVENTID %in% Genome])[-4]/sum(table(Event.annotation$Type[Event.annotation$EVENTID %in% Genome])[-4])
                   )

df$CoreSet <- factor(df$CoreSet, levels = c("Genome", "ASNR", "PanAS", "TissuesAS", "StressAS"))
pal_island <- c("#54C7FC","#F0E68C","#FF83FA","#FFB6C1", "#9AFF9A", "#FF6347")
pal_annotation<-c("#104E8B","#7A378B","#8B3A3A","#778899","#548B54","#CD853F","#8B5F65","#CD9B1D")
brewer.pal(n = 8, name = "Dark2")
pal_annotation <- c("#1B9E77", "#D95F02", "#7570B3", "#A6761D")
plot1 <- ggplot(df,aes(x=CoreSet,y=Percentage,fill=Events,color=Events,width=0.75),size=0.01)+
   geom_bar(stat="identity",color="black",size=0.1)+theme_classic()+labs(x=NULL,y="Percentage")+
   scale_fill_manual(values=pal_annotation)+scale_colour_manual(values=rep("#000000",11))+
   theme(legend.key = element_rect(color="gray",size=0.1),legend.key.size=unit(0.2,"cm"),
         legend.position="right",legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.3, unit="cm"),
         axis.text.x = element_text(color="black",size=5,face="bold"),axis.text.y = element_text(color="black",size=5.5),
         axis.title.x = element_text(size=6,face="bold"),axis.title.y = element_text(size=6.5,face="bold"),
         legend.title = element_text(size=6.5,face="bold"),legend.text = element_text(size=5.5),
         panel.border = element_rect(colour = "black", fill=NA, size=0.2,linetype=3))

# 1 number of items in the sample that are classified as success
# 2 number of items in the universe that are classified as success
# 3 universe size - universe success
# 4 sample size

## ASNR

table(Event.annotation$Type[Event.annotation$EVENTID %in% ASNR])[-4]
sum(table(Event.annotation$Type[Event.annotation$EVENTID %in% ASNR])[-4])
table(Event.annotation$Type[Event.annotation$EVENTID %in% Genome])[-4]
sum(table(Event.annotation$Type[Event.annotation$EVENTID %in% Genome])[-4])

phyper(1438-1,6159,16495-6159, 3838, lower.tail = FALSE) # ASNR vs Genome || IR: 0.43
phyper(313-1,1513,16495-1513, 3838, lower.tail = FALSE) # ASNR vs Genome || ES: 0.99
phyper(1263-1,4556,16495-4556, 3838, lower.tail = FALSE) # ASNR vs Genome || AltAD: 8.401931e-17
phyper(824-1,4267,16495-4267, 3838, lower.tail = FALSE) # ASNR vs Genome || Unk-AS: 1

## PanAS

table(Event.annotation$Type[Event.annotation$EVENTID %in% PanAS.events])[-4]
sum(table(Event.annotation$Type[Event.annotation$EVENTID %in% PanAS.events])[-4])
table(Event.annotation$Type[Event.annotation$EVENTID %in% Genome])[-4]
sum(table(Event.annotation$Type[Event.annotation$EVENTID %in% Genome])[-4])

phyper(780-1,6159,16495-6159, 2082, lower.tail = FALSE) # PanAS vs Genome || IR: 0.45
phyper(168-1,1513,16495-1513, 2082, lower.tail = FALSE) # PanAS vs Genome || ES: 0.97
phyper(688-1,4556,16495-4556, 2082, lower.tail = FALSE) # PanAS vs Genome || AltAD: 3.165927e-09
phyper(446-1,4267,16495-4267, 2082, lower.tail = FALSE) # PanAS vs Genome || Unk-AS: 0.99999

## StressAS

table(Event.annotation$Type[Event.annotation$EVENTID %in% StressAS.events])[-4]
sum(table(Event.annotation$Type[Event.annotation$EVENTID %in% StressAS.events])[-4])
table(Event.annotation$Type[Event.annotation$EVENTID %in% Genome])[-4]
sum(table(Event.annotation$Type[Event.annotation$EVENTID %in% Genome])[-4])

phyper(834-1,6159,16495-6159, 2179, lower.tail = FALSE) # StressAS vs Genome || IR: 0.17
phyper(188-1,1513,16495-1513, 2179, lower.tail = FALSE) # StressAS vs Genome || ES: 0.97
phyper(566-1,4556,16495-4556, 2179, lower.tail = FALSE) # StressAS vs Genome || AltAD: 0.9697
phyper(591-1,4267,16495-4267, 2179, lower.tail = FALSE) # StressAS vs Genome || Unk-AS: 0.07

## Tissues AS 

table(Event.annotation$Type[Event.annotation$EVENTID %in% TissuesAS.events])[-4]
sum(table(Event.annotation$Type[Event.annotation$EVENTID %in% TissuesAS.events])[-4])
table(Event.annotation$Type[Event.annotation$EVENTID %in% Genome])[-4]
sum(table(Event.annotation$Type[Event.annotation$EVENTID %in% Genome])[-4])

phyper(922-1,6159,16495-6159, 2501, lower.tail = FALSE) # StressAS vs Genome || IR: 0.7
phyper(164-1,1513,16495-1513, 2501, lower.tail = FALSE) # StressAS vs Genome || ES: 0.99
phyper(536-1,4556,16495-4556, 2501, lower.tail = FALSE) # StressAS vs Genome || AltAD: 1
phyper(879-1,4267,16495-4267, 2501, lower.tail = FALSE) # StressAS vs Genome || Unk-AS: 2.609771e-29

df$Significance <- c("", "", "", "",
                     "", "*", "", "",
                     "*", "", "", "",
                     "*", "", "", "", 
                     "", "", "", "")

plot1 + annotate("text", x = 2, y = 0.80, label = "*", size = 15) +
   annotate("text", x = 3, y = 0.80, label = "*", size = 15) +
   annotate("text", x = 4, y = 0.15, label = "*", size = 15) 

write.table(df, "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Barplots/EventTypesPerCoreSet.txt", sep = "\t", quote = F, row.names = F, col.names = T)


## Functional impact for each Event type per core set
# CDS or not, CDS-noframeshift, CDS-frameshift-ptc, CDS-frameshifht-clean
# not CDS, CDS-disrupting, CDS-PTC, CDS-change

df2.IR <- data.frame(CoreSet = c(rep("StressAS", 4), rep("TissuesAS", 4), rep("PanAS", 4), rep("ASNR", 4), rep("Genome", 4)),
                 Events = c(rep(c("Not-CDS", "CDS-disrupt", "CDS-PTC", "CDS-change"), 5)),
                 Percentage = NA, pvalue = NA)

df2.ES <- data.frame(CoreSet = c(rep("StressAS", 4), rep("TissuesAS", 4), rep("PanAS", 4), rep("ASNR", 4), rep("Genome", 4)),
                     Events = c(rep(c("Not-CDS", "CDS-disrupt", "CDS-PTC", "CDS-change"), 5)),
                     Percentage = NA, pvalue = NA)

df2.AltAD <- data.frame(CoreSet = c(rep("StressAS", 4), rep("TissuesAS", 4), rep("PanAS", 4), rep("ASNR", 4), rep("Genome", 4)),
                     Events = c(rep(c("Not-CDS", "CDS-disrupt", "CDS-PTC", "CDS-change"), 5)),
                     Percentage = NA, pvalue = NA)

df2.AS <- data.frame(CoreSet = c(rep("StressAS", 4), rep("TissuesAS", 4), rep("PanAS", 4), rep("ASNR", 4), rep("Genome", 4)),
                     Events = c(rep(c("Not-CDS", "CDS-disrupt", "CDS-PTC", "CDS-change"), 5)),
                     Percentage = NA, pvalue = NA)

Fimpact <- list(ES = df2.ES, IR = df2.IR, AltAD = df2.AltAD, AS  = df2.AS)
Event.annotation.Fimpact <- Event.annotation
Event.annotation.Fimpact$CDS[is.na(Event.annotation.Fimpact$CDS)] <- "False"
Event.annotation.Fimpact$Fimpact <- NA
Event.annotation.Fimpact$Fimpact[which(Event.annotation.Fimpact$CDS == "False")] <- "Not-CDS"
Event.annotation.Fimpact$Fimpact[which(Event.annotation.Fimpact$CDS == "True" & Event.annotation.Fimpact$Frameshift == "False")] <- "CDS-disrupt"
Event.annotation.Fimpact$Fimpact[which(Event.annotation.Fimpact$CDS == "True" & Event.annotation.Fimpact$Frameshift == "True" & Event.annotation.Fimpact$PTC == "True")] <- "CDS-PTC"
Event.annotation.Fimpact$Fimpact[which(Event.annotation.Fimpact$CDS == "True" & Event.annotation.Fimpact$Frameshift == "True" & Event.annotation.Fimpact$PTC == "False")] <- "CDS-change"
table(Event.annotation.Fimpact$Fimpact)

for(i in names(Fimpact)){
   
   StressAS.sub <- table(Event.annotation.Fimpact$Fimpact[Event.annotation.Fimpact$EVENTID %in% StressAS.events & Event.annotation$Type == i])/sum(table(Event.annotation.Fimpact$Fimpact[Event.annotation.Fimpact$EVENTID %in% StressAS.events & Event.annotation$Type == i]))
   TissuesAS.sub <- table(Event.annotation.Fimpact$Fimpact[Event.annotation.Fimpact$EVENTID %in% TissuesAS.events & Event.annotation$Type == i])/sum(table(Event.annotation.Fimpact$Fimpact[Event.annotation.Fimpact$EVENTID %in% TissuesAS.events & Event.annotation$Type == i]))
   PanAS.sub <- table(Event.annotation.Fimpact$Fimpact[Event.annotation.Fimpact$EVENTID %in% PanAS.events & Event.annotation$Type == i])/sum(table(Event.annotation.Fimpact$Fimpact[Event.annotation.Fimpact$EVENTID %in% PanAS.events & Event.annotation$Type == i]))
   ASNR.sub <- table(Event.annotation.Fimpact$Fimpact[Event.annotation.Fimpact$EVENTID %in% ASNR & Event.annotation$Type == i])/sum(table(Event.annotation.Fimpact$Fimpact[Event.annotation.Fimpact$EVENTID %in% ASNR & Event.annotation$Type == i]))
   Genome.sub  <- table(Event.annotation.Fimpact$Fimpact[Event.annotation.Fimpact$EVENTID %in% Genome & Event.annotation$Type == i])/sum(table(Event.annotation.Fimpact$Fimpact[Event.annotation.Fimpact$EVENTID %in% Genome & Event.annotation$Type == i]))

   totalStressAS <- sum(table(Event.annotation.Fimpact$Fimpact[Event.annotation.Fimpact$EVENTID %in% StressAS.events & Event.annotation$Type == i]))
   totalTissuesAS <- sum(table(Event.annotation.Fimpact$Fimpact[Event.annotation.Fimpact$EVENTID %in% TissuesAS.events & Event.annotation$Type == i]))
   totalPanAS <- sum(table(Event.annotation.Fimpact$Fimpact[Event.annotation.Fimpact$EVENTID %in% PanAS.events & Event.annotation$Type == i]))
   totalASNR <- sum(table(Event.annotation.Fimpact$Fimpact[Event.annotation.Fimpact$EVENTID %in% ASNR & Event.annotation$Type == i]))
   totalGenome <- sum(table(Event.annotation.Fimpact$Fimpact[Event.annotation.Fimpact$EVENTID %in% Genome & Event.annotation$Type == i]))
   
   Fimpact[[i]]$Type <- rep(i, nrow(Fimpact[[i]])) 
   
   ## Genome 
   
   if(length(Genome.sub[names(Genome.sub) %in% "Not-CDS"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"] <- Genome.sub[names(Genome.sub) %in% "Not-CDS"]
   }
   
   
   if(length(Genome.sub[names(Genome.sub) %in% "CDS-PTC"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"] <- Genome.sub[names(Genome.sub) %in% "CDS-PTC"]
   }
   
   if(length(Genome.sub[names(Genome.sub) %in% "CDS-disrupt"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"] <- Genome.sub[names(Genome.sub) %in% "CDS-disrupt"]
   }
   
   if(length(Genome.sub[names(Genome.sub) %in% "CDS-change"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"] <- Genome.sub[names(Genome.sub) %in% "CDS-change"]
   }
   
   
   ## StressAS
   
   if(length(StressAS.sub[names(StressAS.sub) %in% "Not-CDS"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "StressAS" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "StressAS" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"] <- StressAS.sub[names(StressAS.sub) %in% "Not-CDS"]
   }
  
   Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "StressAS" & Fimpact[[i]]$Events == "Not-CDS"), "pvalue"] <- phyper(Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "StressAS" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"]*totalStressAS-1,
                                                                                                                  Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"]*totalGenome,
                                                                                                                  totalGenome-Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"]*totalGenome,
                                                                                                                  totalStressAS, lower.tail = FALSE)
   
   
   if(length(StressAS.sub[names(StressAS.sub) %in% "CDS-PTC"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "StressAS" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "StressAS" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"] <- StressAS.sub[names(StressAS.sub) %in% "CDS-PTC"]
   }
   
   Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "StressAS" & Fimpact[[i]]$Events == "CDS-PTC"), "pvalue"] <- phyper(Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "StressAS" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"]*totalStressAS-1,
                                                                                                                  Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"]*totalGenome,
                                                                                                                  totalGenome-Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"]*totalGenome,
                                                                                                                  totalStressAS, lower.tail = FALSE)
   
   if(length(StressAS.sub[names(StressAS.sub) %in% "CDS-disrupt"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "StressAS" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "StressAS" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"] <- StressAS.sub[names(StressAS.sub) %in% "CDS-disrupt"]
   }
   
   Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "StressAS" & Fimpact[[i]]$Events == "CDS-disrupt"), "pvalue"] <- phyper(Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "StressAS" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"]*totalStressAS-1,
                                                                                                                  Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"]*totalGenome,
                                                                                                                  totalGenome-Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"]*totalGenome,
                                                                                                                  totalStressAS, lower.tail = FALSE)
   
   if(length(StressAS.sub[names(StressAS.sub) %in% "CDS-change"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "StressAS" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "StressAS" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"] <- StressAS.sub[names(StressAS.sub) %in% "CDS-change"]
   }
   
   Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "StressAS" & Fimpact[[i]]$Events == "CDS-change"), "pvalue"] <- phyper(Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "StressAS" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"]*totalStressAS-1,
                                                                                                                  Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"]*totalGenome,
                                                                                                                  totalGenome-Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"]*totalGenome,
                                                                                                                  totalStressAS, lower.tail = FALSE)
   
   ## TissuesAS
   
   if(length(TissuesAS.sub[names(TissuesAS.sub) %in% "Not-CDS"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "TissuesAS" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "TissuesAS" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"] <- TissuesAS.sub[names(TissuesAS.sub) %in% "Not-CDS"]
   }
   
   Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "TissuesAS" & Fimpact[[i]]$Events == "Not-CDS"), "pvalue"] <- phyper(Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "TissuesAS" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"]*totalTissuesAS-1,
                                                                                                                  Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"]*totalGenome,
                                                                                                                  totalGenome-Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"]*totalGenome,
                                                                                                                  totalTissuesAS, lower.tail = FALSE)
   
   
   if(length(TissuesAS.sub[names(TissuesAS.sub) %in% "CDS-PTC"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "TissuesAS" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "TissuesAS" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"] <- TissuesAS.sub[names(TissuesAS.sub) %in% "CDS-PTC"]
   }
   
   Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "TissuesAS" & Fimpact[[i]]$Events == "CDS-PTC"), "pvalue"] <- phyper(Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "TissuesAS" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"]*totalTissuesAS-1,
                                                                                                                  Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"]*totalGenome,
                                                                                                                  totalGenome-Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"]*totalGenome,
                                                                                                                  totalTissuesAS, lower.tail = FALSE)
   
   if(length(TissuesAS.sub[names(TissuesAS.sub) %in% "CDS-disrupt"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "TissuesAS" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "TissuesAS" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"] <- TissuesAS.sub[names(TissuesAS.sub) %in% "CDS-disrupt"]
   }
   
   Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "TissuesAS" & Fimpact[[i]]$Events == "CDS-disrupt"), "pvalue"] <- phyper(Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "TissuesAS" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"]*totalTissuesAS-1,
                                                                                                                      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"]*totalGenome,
                                                                                                                      totalGenome-Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"]*totalGenome,
                                                                                                                      totalTissuesAS, lower.tail = FALSE)
   
   if(length(TissuesAS.sub[names(TissuesAS.sub) %in% "CDS-change"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "TissuesAS" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "TissuesAS" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"] <- TissuesAS.sub[names(TissuesAS.sub) %in% "CDS-change"]
   }
   
   Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "TissuesAS" & Fimpact[[i]]$Events == "CDS-change"), "pvalue"] <- phyper(Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "TissuesAS" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"]*totalTissuesAS-1,
                                                                                                                     Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"]*totalGenome,
                                                                                                                     totalGenome-Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"]*totalGenome,
                                                                                                                     totalTissuesAS, lower.tail = FALSE)
   
   ## PanAS
   
   if(length(PanAS.sub[names(PanAS.sub) %in% "Not-CDS"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "PanAS" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "PanAS" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"] <- PanAS.sub[names(PanAS.sub) %in% "Not-CDS"]
   }
   
   Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "PanAS" & Fimpact[[i]]$Events == "Not-CDS"), "pvalue"] <- phyper(Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "PanAS" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"]*totalPanAS-1,
                                                                                                                  Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"]*totalGenome,
                                                                                                                  totalGenome-Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"]*totalGenome,
                                                                                                                  totalPanAS, lower.tail = FALSE)
   
   
   if(length(PanAS.sub[names(PanAS.sub) %in% "CDS-PTC"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "PanAS" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "PanAS" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"] <- PanAS.sub[names(PanAS.sub) %in% "CDS-PTC"]
   }
   
   Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "PanAS" & Fimpact[[i]]$Events == "CDS-PTC"), "pvalue"] <- phyper(Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "PanAS" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"]*totalPanAS-1,
                                                                                                                  Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"]*totalGenome,
                                                                                                                  totalGenome-Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"]*totalGenome,
                                                                                                                  totalPanAS, lower.tail = FALSE)
   
   if(length(PanAS.sub[names(PanAS.sub) %in% "CDS-disrupt"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "PanAS" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "PanAS" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"] <- PanAS.sub[names(PanAS.sub) %in% "CDS-disrupt"]
   }
   
   Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "PanAS" & Fimpact[[i]]$Events == "CDS-disrupt"), "pvalue"] <- phyper(Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "PanAS" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"]*totalPanAS-1,
                                                                                                                      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"]*totalGenome,
                                                                                                                      totalGenome-Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"]*totalGenome,
                                                                                                                      totalPanAS, lower.tail = FALSE)
   
   if(length(PanAS.sub[names(PanAS.sub) %in% "CDS-change"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "PanAS" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "PanAS" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"] <- PanAS.sub[names(PanAS.sub) %in% "CDS-change"]
   }
   
   Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "PanAS" & Fimpact[[i]]$Events == "CDS-change"), "pvalue"] <- phyper(Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "PanAS" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"]*totalPanAS-1,
                                                                                                                     Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"]*totalGenome,
                                                                                                                     totalGenome-Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"]*totalGenome,
                                                                                                                     totalPanAS, lower.tail = FALSE)
   
   ## ASNR
   
   if(length(ASNR.sub[names(ASNR.sub) %in% "Not-CDS"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "ASNR" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "ASNR" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"] <- ASNR.sub[names(ASNR.sub) %in% "Not-CDS"]
   }
   
   Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "ASNR" & Fimpact[[i]]$Events == "Not-CDS"), "pvalue"] <- phyper(Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "ASNR" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"]*totalASNR-1,
                                                                                                                  Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"]*totalGenome,
                                                                                                                  totalGenome-Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "Not-CDS"), "Percentage"]*totalGenome,
                                                                                                                  totalASNR, lower.tail = FALSE)
   
   if(length(ASNR.sub[names(ASNR.sub) %in% "CDS-PTC"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "ASNR" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "ASNR" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"] <- ASNR.sub[names(ASNR.sub) %in% "CDS-PTC"]
   }
   
   Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "ASNR" & Fimpact[[i]]$Events == "CDS-PTC"), "pvalue"] <- phyper(Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "ASNR" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"]*totalASNR-1,
                                                                                                                  Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"]*totalGenome,
                                                                                                                  totalGenome-Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-PTC"), "Percentage"]*totalGenome,
                                                                                                                  totalASNR, lower.tail = FALSE)
   
   if(length(ASNR.sub[names(ASNR.sub) %in% "CDS-disrupt"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "ASNR" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "ASNR" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"] <- ASNR.sub[names(ASNR.sub) %in% "CDS-disrupt"]
   }
   
   Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "ASNR" & Fimpact[[i]]$Events == "CDS-disrupt"), "pvalue"] <- phyper(Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "ASNR" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"]*totalASNR-1,
                                                                                                                      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"]*totalGenome,
                                                                                                                      totalGenome-Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-disrupt"), "Percentage"]*totalGenome,
                                                                                                                      totalASNR, lower.tail = FALSE)
   
   if(length(ASNR.sub[names(ASNR.sub) %in% "CDS-change"]) == 0){
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "ASNR" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"] <- 0
   }else {
      Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "ASNR" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"] <- ASNR.sub[names(ASNR.sub) %in% "CDS-change"]
   }
   
   Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "ASNR" & Fimpact[[i]]$Events == "CDS-change"), "pvalue"] <- phyper(Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "ASNR" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"]*totalASNR-1,
                                                                                                                     Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"]*totalGenome,
                                                                                                                     totalGenome-Fimpact[[i]][which(Fimpact[[i]]$CoreSet == "Genome" & Fimpact[[i]]$Events == "CDS-change"), "Percentage"]*totalGenome,
                                                                                                                     totalASNR, lower.tail = FALSE)
   
   
}

require(ggsci)
require(scales)
show_col(pal_npg("nrc")(4))
pal_cds <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF")

plot2 <- ggplot(Fimpact$ES,aes(x=CoreSet,y=Percentage,fill=Events,color=Events,width=0.75),size=0.01)+
   geom_bar(stat="identity",color="black",size=0.1)+theme_classic()+labs(x=NULL,y="Percentage")+
   scale_fill_manual(values=pal_cds)+scale_colour_manual(values=rep("#000000",11))+
   theme(legend.key = element_rect(color="gray",size=0.1),legend.key.size=unit(0.2,"cm"),
         legend.position="right",legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.3, unit="cm"),
         axis.text.x = element_text(color="black",size=5,face="bold"),axis.text.y = element_text(color="black",size=5.5),
         axis.title.x = element_text(size=6,face="bold"),axis.title.y = element_text(size=6.5,face="bold"),
         legend.title = element_text(size=6.5,face="bold"),legend.text = element_text(size=5.5),
         panel.border = element_rect(colour = "black", fill=NA, size=0.2,linetype=3)) +
   annotate("text", x = 4, y = 0.285, label = "*", size = 15) + annotate("text", x = 5, y = 0.32, label = "*", size = 15) + coord_flip()

plot3 <- ggplot(Fimpact$IR,aes(x=CoreSet,y=Percentage,fill=Events,color=Events,width=0.75),size=0.01)+
   geom_bar(stat="identity",color="black",size=0.1)+theme_classic()+labs(x=NULL,y="Percentage")+
   scale_fill_manual(values=pal_cds)+scale_colour_manual(values=rep("#000000",11))+
   theme(legend.key = element_rect(color="gray",size=0.1),legend.key.size=unit(0.2,"cm"),
         legend.position="right",legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.3, unit="cm"),
         axis.text.x = element_text(color="black",size=5,face="bold"),axis.text.y = element_text(color="black",size=5.5),
         axis.title.x = element_text(size=6,face="bold"),axis.title.y = element_text(size=6.5,face="bold"),
         legend.title = element_text(size=6.5,face="bold"),legend.text = element_text(size=5.5),
         panel.border = element_rect(colour = "black", fill=NA, size=0.2,linetype=3)) +
   annotate("text", x = 4, y = 0.581/2, label = "*", size = 15) + annotate("text", x = 5, y = 0.577/2, label = "*", size = 15) + 
   annotate("text", x = 3, y = 0.9, label = "*", size = 15) + coord_flip()

plot4 <- ggplot(Fimpact$AltAD,aes(x=CoreSet,y=Percentage,fill=Events,color=Events,width=0.75),size=0.01)+
   geom_bar(stat="identity",color="black",size=0.1)+theme_classic()+labs(x=NULL,y="Percentage")+
   scale_fill_manual(values=pal_cds)+scale_colour_manual(values=rep("#000000",11))+
   theme(legend.key = element_rect(color="gray",size=0.1),legend.key.size=unit(0.2,"cm"),
         legend.position="right",legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.3, unit="cm"),
         axis.text.x = element_text(color="black",size=5,face="bold"),axis.text.y = element_text(color="black",size=5.5),
         axis.title.x = element_text(size=6,face="bold"),axis.title.y = element_text(size=6.5,face="bold"),
         legend.title = element_text(size=6.5,face="bold"),legend.text = element_text(size=5.5),
         panel.border = element_rect(colour = "black", fill=NA, size=0.2,linetype=3)) +
   annotate("text", x = 4, y = 0.31, label = "*", size = 15) + annotate("text", x = 5, y = 0.325, label = "*", size = 15) + 
   annotate("text", x = 3, y = 0.8, label = "*", size = 15) + annotate("text", x = 1, y = 0.79, label = "*", size = 15) + coord_flip()

write.table(do.call(rbind, Fimpact), "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Barplots/CDSimpact.txt", sep = "\t", col.names = T, row.names = F, quote = F)

#### 3.2 Intersection between AS core sets & GE-AS sets ####

AScore.list <- list(PanAS = PanAS.events, TissuesAS = TissuesAS.events, StressAS = StressAS.events)
nVennR::plotVenn(sets = AScore.list, opacity = 0.2, borderWidth = 3, systemShow = F, fontScale = 2, outFile = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Venns/AScore_nVennR.svg")

GEcore.list <- list(PanGE = PanGE.genes, TissuesGE = TissuesGE.genes, StressGE = StressGE.genes)
nVennR::plotVenn(sets = GEcore.list, opacity = 0.2, borderWidth = 3, systemShow = F, fontScale = 2, outFile = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Venns/GEcore_nVennR.svg")

ASGEcore.list <- lapply(AScore.list, function(x){
   tt <- Event.annotation.split$GENEID[Event.annotation.split$EVENTID %in% x]
   tt <- unique(tt[!is.na(tt)])
   tt
})

plotVenn(sets = list(PanAS = ASGEcore.list$PanAS, PanGE = GEcore.list$PanGE), opacity = 0.2, borderWidth = 3, systemShow = F, fontScale = 2, outFile = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Venns/PanASGE_nVennR.svg")
plotVenn(sets = list(TissuesAS = ASGEcore.list$TissuesAS, TissuesGE = GEcore.list$TissuesGE), opacity = 0.2, borderWidth = 3, systemShow = F, fontScale = 2, outFile = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Venns/TissuesASGE_nVennR.svg")
plotVenn(sets = list(StressAS = ASGEcore.list$StressAS, StressGE = GEcore.list$StressGE), opacity = 0.2, borderWidth = 3, systemShow = F, fontScale = 2, outFile = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Venns/StressASGE_nVennR.svg")

### 3.3 Variation plot ###

reference <- read.table("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Variation/Stress_vs_Tissues-input_table.tab", header = T)
Final_PSI.df <- read.table("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/Final_PSI_lowFilter_PEandTE_batchCorrectedTechStudy_Normalized_PresenceRange.txt", header = T)
View(Final_PSI.df)

Variation <- data.frame(EventID = Final_PSI.df$events.names, Range_tissues = NA, Max_dPSI_stress = NA, N_tis = NA, N_stress = NA)

### For each AS event:
## tissues at least two replicates for each tissue with sufficient read coverage (not NA)
## stress at least two replicates for stress and control with sufficient read coverage in at least 3/5 experiments (not NA) 
## PSI range across tissues (Range_tis): difference between median PSI values of the tissue types with the highest and lowest values
## Max_dPSI_stress: the highest absolute dPSI value among the stress experiments

TissuesSamples <- c(1:16, 18:52)
StressSamples <- c(19,20,22,23,25:36, 54:142)

for(i in 1:nrow(Variation)){
   
   # initial row
   
   row.info <- as.vector(Final_PSI.df[which(Final_PSI.df$events.names == Variation$EventID[i]),])
   
   ## Tissues
   
   # tissues row
   
   row.info.tissues <- row.info[TissuesSamples+1]
   names(row.info.tissues) <- c(rep("Bud", 16), rep("Needle", 27), rep("Vascular", 8))
   row.info.tissues <- unlist(row.info.tissues)[!is.na(unlist(row.info.tissues))]
   
   # presence filter
   
   N_tis <- length(which((table(names(row.info.tissues)) >= 2) == TRUE))
   presence <- N_tis == 3
   if(!presence){
      Variation$Range_tissues[i] <- NA
      Variation$N_tis[i] <- N_tis
   }else if(presence){
      
      # PSI range across tissues
      
      Variation$N_tis[i] <- N_tis
      agg <- aggregate(row.info.tissues, by=list(Tissues = names(row.info.tissues)), median)
      Variation$Range_tissues[i] <- max(agg$x) - min(agg$x)
   }
   
   ## Stress
   
   # Stress row
   
   row.info.stress <- row.info[StressSamples]
   names(row.info.stress) <- c(rep("DO_Control", 2),
                               rep("HS_Control", 2),
                               rep("PH_Control", 12),
                               rep("HS_Stress", 4),
                               rep("DO_Stress", 10),
                               rep("PH_Stress", 6),
                               rep("FU1_Control", 3),
                               rep("FU2_Control", 18),
                               rep("FU1_Stress", 4),
                               rep("FU2_Stress", 44))
   row.info.stress <- unlist(row.info.stress)[!is.na(unlist(row.info.stress))]
   
   # presence and dPSI computation
   
   presence.list <- list()
   dPSI.list <- list()
   
   for(j in c("DO", "HS", "PH", "FU1", "FU2")){
      
      experiment.sub <- row.info.stress[grep(j, names(row.info.stress))]
      presence.list[[j]] <- length(which((table(names(experiment.sub)) >= 2) == TRUE)) == 2
      if(!presence.list[[j]]){
         dPSI.list[[j]] <- NA
      }else if(presence.list[[j]]){
         agg2 <- aggregate(experiment.sub, by=list(Stress = names(experiment.sub)), mean)
         dPSI.list[[j]] <- abs(agg2$x[which(agg2$Stress == paste0(j, "_Stress"))] - agg2$x[which(agg2$Stress == paste0(j, "_Control"))])
      }
      
   } # for loop stress experiments
   
   # presence filter and dPSI assignment
   
   N_stress <- length(which((unlist(presence.list)) == TRUE))
   presence <- N_stress == 5
   if(!presence){
    Variation$N_stress[i] <- N_stress
    Variation$Max_dPSI_stress[i] <- NA
   }else if(presence){
      Variation$N_stress[i] <- N_stress
      Variation$Max_dPSI_stress[i] <- max(unlist(dPSI.list)[!is.na(unlist(dPSI.list))])
   }
   
} # for loop AS events

View(Variation)

Variation$Range_tissues <- Variation$Range_tissues*100
Variation$Max_dPSI_stress <- Variation$Max_dPSI_stress*100

Variation$GlobalPSI <- Variation$Range_tissues + Variation$Max_dPSI_stress
Variation$T_prop <- (Variation$Range_tissues/Variation$GlobalPSI)*100
Variation$S_prop <- (Variation$Max_dPSI_stress/Variation$GlobalPSI)*100

i <- 10
j <- 3

Tissues <- na.omit(as.vector(Variation[Variation$GlobalPSI>i & Variation$N_tis>=j & Variation$N_stress>=j,7]))
Stress <- na.omit(as.vector(Variation[Variation$GlobalPSI>i & Variation$N_tis>=j & Variation$N_stress>=j,8]))

binPra_S <- (hexbin(Stress,Tissues,xbins = 20))
plot (binPra_S, main="" , colramp=coolwarm)

# other species mix abiotic with biotic

Table <- reference 
table(Table$N_tis[which(Table$Species == "Hsa")])
table(Table$N_abiotic[which(Table$Species == "Hsa")])
table(Table$N_biotic[which(Table$Species == "Hsa")])
Table <- Table[which(Table$N_tis == 5 & Table$N_abiotic == 5 & Table$N_biotic == 5), ]
Table$Max_dPSI_stress <- apply(Table[,3:4], 1 , max)

Table$GlobalPSI <- Table$Range_tis + Table$Max_dPSI_stress
Table$T_prop <- (Table$Range_tis/Table$GlobalPSI)*100
Table$S_prop <- (Table$Max_dPSI_stress/Table$GlobalPSI)*100

Ath_T <- na.omit(as.vector(Table[Table$Species=="Ath" & Table$GlobalPSI>i ,11]))
Ath_S <- na.omit(as.vector(Table[Table$Species=="Ath" & Table$GlobalPSI>i ,12]))

Cel_T <- na.omit(as.vector(Table[Table$Species=="Cel" & Table$GlobalPSI>i ,11]))
Cel_S <- na.omit(as.vector(Table[Table$Species=="Cel" & Table$GlobalPSI>i ,12]))

Dme_T <- na.omit(as.vector(Table[Table$Species=="Dme" & Table$GlobalPSI>i ,11]))
Dme_S <- na.omit(as.vector(Table[Table$Species=="Dme" & Table$GlobalPSI>i ,12]))

Hsa_T <- na.omit(as.vector(Table[Table$Species=="Hsa" & Table$GlobalPSI>i ,11]))
Hsa_S <- na.omit(as.vector(Table[Table$Species=="Hsa" & Table$GlobalPSI>i ,12]))

binAth_S <- (hexbin(Ath_S, Ath_T, xbins = 20))
plot (binAth_S, main="", colramp=coolwarm)

binCel_S <- (hexbin(Cel_S, Cel_T, xbins = 20))
plot (binCel_S, main="", colramp=coolwarm)

binDme_S <- (hexbin(Dme_S, Dme_T, xbins = 20))
plot (binDme_S, main="", colramp=coolwarm)

binHsa_S <- (hexbin(Hsa_S, Hsa_T, xbins = 20))
plot (binHsa_S, main="", colramp=coolwarm)

Ath_B <- na.omit(as.vector(Table[Table$Species=="Ath" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,12]))
Cel_T <- na.omit(as.vector(Table[Table$Species=="Cel" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,10]))
Cel_A <- na.omit(as.vector(Table[Table$Species=="Cel" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,11]))
Cel_B <- na.omit(as.vector(Table[Table$Species=="Cel" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,12]))
Dme_T <- na.omit(as.vector(Table[Table$Species=="Dme" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,10]))
Dme_A <- na.omit(as.vector(Table[Table$Species=="Dme" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,11]))
Dme_B <- na.omit(as.vector(Table[Table$Species=="Dme" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,12]))
Hsa_T <- na.omit(as.vector(Table[Table$Species=="Hsa" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,10]))
Hsa_A <- na.omit(as.vector(Table[Table$Species=="Hsa" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,11]))
Hsa_B <- na.omit(as.vector(Table[Table$Species=="Hsa" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,12]))

write.table(Table, "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Variation/ReferenceTable_Variation.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(Variation, "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Variation/Pra_Variation.txt", col.names = T, row.names = F, quote = F, sep = "\t")


#### 3.3 Prepare for matt & variation ####

## prepare taeda genome for matt | 
## NOT TAKING INTO ACCOUNT AS TYPE: exons and introns difference for each core set; if similar exons and introns difference for Up and Down in all sets
## exon features could affect IR etc
## compare exons and introns for each core set and event type (AltAD; ES; IR; Splicing)
## For multiple; multiple genes are taken into account

## Compare Exons for ES; Compare Introns for IR; Compare exons for AltAD (probably expand to both intron and exons for both)
## Genome vs (ASNR; PanAS; StressAS; TissuesAS) (compared to genome)
## split per core set and group

Event.annotation.split <- read.table("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/Event_annotation.txt")
PanAS <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/PanAS_EventsID.txt")
TissuesAS <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/TissuesAS_EventsID.txt")
StressAS <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/StressAS_EventsID.txt")
ASNR  <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/ASNR_EventsID.txt")
GenomeAS  <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/Genome_EventsID.txt")

## IR

Event.annotation.split.IR <- Event.annotation.split
Event.annotation.split.IR <- Event.annotation.split.IR[which(Event.annotation.split.IR$Type == "IR"),]
for(i in 1:nrow(Event.annotation.split.IR)){
   Event.id <- Event.annotation.split.IR$EVENTID[i]
   character.list <- list()
   
   Pan <- length(which(PanAS$x == Event.id)) == 0
   if(Pan){
      character.list[["PanAS"]] <- NA
   }else if(!Pan){
      character.list[["PanAS"]] <- "PanAS"
   }
   
   Stress <- length(which(StressAS$x == Event.id)) == 0
   if(Stress){
      character.list[["StressAS"]] <- NA
   }else if(!Stress){
      character.list[["StressAS"]] <- "StressAS"
   }
   
   Tissues <- length(which(TissuesAS$x == Event.id)) == 0
   if(Tissues){
      character.list[["TissuesAS"]] <- NA
   }else if(!Tissues){
      character.list[["TissuesAS"]] <- "TissuesAS"
   }
   
   ASNRs <- length(which(ASNR$x == Event.id)) == 0
   if(ASNRs){
      character.list[["ASNR"]] <- NA
   }else if(!Tissues){
      character.list[["ASNR"]] <- "ASNR"
   }
   
   
   Genome <- length(which(GenomeAS$x == Event.id)) == 0
   if(Genome){
      character.list[["Genome"]] <- NA
   }else if(!Genome){
      character.list[["Genome"]] <- "Genome"
   }
   
   group <- unlist(character.list)[!is.na(unlist(character.list))]
   if(length(group) == 0){
      Event.annotation.split.IR$GROUP[i] <- NA
   }else if(length(group) > 0){
      Event.annotation.split.IR$GROUP[i] <- paste(group, collapse = ",")
   }
   
}

table(Event.annotation.split.IR$GROUP)
Event.annotation.split.IR <- separate_longer_delim(Event.annotation.split.IR, cols = "GROUP", delim = ",")
table(Event.annotation.split.IR$GROUP)
Event.annotation.split.IR <- Event.annotation.split.IR[!(is.na(Event.annotation.split.IR$GENEID)),]
Event.annotation.split.IR <- Event.annotation.split.IR[!(is.na(Event.annotation.split.IR$GROUP)),]
table(Event.annotation.split.IR$GROUP)

Event.annotation.split.IR$EVENTID <- make.unique(Event.annotation.split.IR$EVENTID)
Matt_IR <- Event.annotation.split.IR[,c(1:7,14)]
write.table(Matt_IR, "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Matt_compare/Matt_IR.tab", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(data.frame(unique(Matt_IR$SCAFF)), "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Matt_compare/IR_Seqnames.txt", sep = "\t", row.names = F, col.names = F, quote = F)
Matt_IR_results <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Matt_compare/vsIRs/overview_OnlySig.txt", header = T)
Matt_IR_results.filtered <- Matt_IR_results[which(Matt_IR_results$PVAL_MANN_WHITNEY_U_TEST_Genome_VS_PanAS <= 0.05 | Matt_IR_results$PVAL_MANN_WHITNEY_U_TEST_Genome_VS_StressAS <= 0.05 | Matt_IR_results$PVAL_MANN_WHITNEY_U_TEST_Genome_VS_TissuesAS <= 0.05 | Matt_IR_results$PVAL_MANN_WHITNEY_U_TEST_ASNR_VS_Genome <= 0.05),]
View(Matt_IR_results.filtered
     )
write.table(x = Matt_IR_results.filtered, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Matt_compare/vsIRs/Matt_IR_results_filtered.txt", quote = F, sep = "\t", col.names = T, row.names = F)

## ES

Event.annotation.split.ES <- Event.annotation.split
Event.annotation.split.ES <- Event.annotation.split.ES[which(Event.annotation.split.ES$Type == "ES"),]
for(i in 1:nrow(Event.annotation.split.ES)){
   Event.id <- Event.annotation.split.ES$EVENTID[i]
   character.list <- list()
   
   Pan <- length(which(PanAS$x == Event.id)) == 0
   if(Pan){
      character.list[["PanAS"]] <- NA
   }else if(!Pan){
      character.list[["PanAS"]] <- "PanAS"
   }
   
   Stress <- length(which(StressAS$x == Event.id)) == 0
   if(Stress){
      character.list[["StressAS"]] <- NA
   }else if(!Stress){
      character.list[["StressAS"]] <- "StressAS"
   }
   
   Tissues <- length(which(TissuesAS$x == Event.id)) == 0
   if(Tissues){
      character.list[["TissuesAS"]] <- NA
   }else if(!Tissues){
      character.list[["TissuesAS"]] <- "TissuesAS"
   }
   
   ASNRs <- length(which(ASNR$x == Event.id)) == 0
   if(ASNRs){
      character.list[["ASNR"]] <- NA
   }else if(!Tissues){
      character.list[["ASNR"]] <- "ASNR"
   }
   
   
   Genome <- length(which(GenomeAS$x == Event.id)) == 0
   if(Genome){
      character.list[["Genome"]] <- NA
   }else if(!Genome){
      character.list[["Genome"]] <- "Genome"
   }
   
   group <- unlist(character.list)[!is.na(unlist(character.list))]
   if(length(group) == 0){
      Event.annotation.split.ES$GROUP[i] <- NA
   }else if(length(group) > 0){
      Event.annotation.split.ES$GROUP[i] <- paste(group, collapse = ",")
   }
   
}

table(Event.annotation.split.ES$GROUP)
Event.annotation.split.ES <- separate_longer_delim(Event.annotation.split.ES, cols = "GROUP", delim = ",")
table(Event.annotation.split.ES$GROUP)
Event.annotation.split.ES <- Event.annotation.split.ES[!(is.na(Event.annotation.split.ES$GENEID)),]
Event.annotation.split.ES <- Event.annotation.split.ES[!(is.na(Event.annotation.split.ES$GROUP)),]
table(Event.annotation.split.ES$GROUP)

Event.annotation.split.ES$EVENTID <- make.unique(Event.annotation.split.ES$EVENTID)
Matt_ES <- Event.annotation.split.ES[,c(1:7,14)]
write.table(Matt_ES, "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Matt_compare/vsES/Matt_ES.tab", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(data.frame(unique(Matt_ES$SCAFF)), "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Matt_compare/vsES/ES_Seqnames.txt", sep = "\t", row.names = F, col.names = F, quote = F)
Matt_ES_results <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Matt_compare/vsES/overview_of_features_and_comparisons.tab", header = T)
Matt_ES_results <- Matt_ES_results[,c(1,22,23, 30:35)]
Matt_ES_results.filtered <- Matt_ES_results[which(Matt_ES_results$PVAL_MANN_WHITNEY_U_TEST_Genome_VS_PanAS <= 0.05 | Matt_ES_results$PVAL_MANN_WHITNEY_U_TEST_Genome_VS_StressAS <= 0.05 | Matt_ES_results$PVAL_MANN_WHITNEY_U_TEST_Genome_VS_TissuesAS <= 0.05 | Matt_ES_results$PVAL_MANN_WHITNEY_U_TEST_ASNR_VS_Genome <= 0.05),]
View(Matt_ES_results.filtered)
write.table(x = Matt_ES_results.filtered, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Matt_compare/vsES/Matt_ES_results_filtered.txt", quote = F, sep = "\t", col.names = T, row.names = F)

## AltAD

Event.annotation.split.AltAD <- Event.annotation.split
Event.annotation.split.AltAD <- Event.annotation.split.AltAD[which(Event.annotation.split.AltAD$Type == "AltAD"),]
for(i in 1:nrow(Event.annotation.split.AltAD)){
   Event.id <- Event.annotation.split.AltAD$EVENTID[i]
   character.list <- list()
   
   Pan <- length(which(PanAS$x == Event.id)) == 0
   if(Pan){
      character.list[["PanAS"]] <- NA
   }else if(!Pan){
      character.list[["PanAS"]] <- "PanAS"
   }
   
   Stress <- length(which(StressAS$x == Event.id)) == 0
   if(Stress){
      character.list[["StressAS"]] <- NA
   }else if(!Stress){
      character.list[["StressAS"]] <- "StressAS"
   }
   
   Tissues <- length(which(TissuesAS$x == Event.id)) == 0
   if(Tissues){
      character.list[["TissuesAS"]] <- NA
   }else if(!Tissues){
      character.list[["TissuesAS"]] <- "TissuesAS"
   }
   
   ASNRs <- length(which(ASNR$x == Event.id)) == 0
   if(ASNRs){
      character.list[["ASNR"]] <- NA
   }else if(!Tissues){
      character.list[["ASNR"]] <- "ASNR"
   }
   
   
   Genome <- length(which(GenomeAS$x == Event.id)) == 0
   if(Genome){
      character.list[["Genome"]] <- NA
   }else if(!Genome){
      character.list[["Genome"]] <- "Genome"
   }
   
   group <- unlist(character.list)[!is.na(unlist(character.list))]
   if(length(group) == 0){
      Event.annotation.split.AltAD$GROUP[i] <- NA
   }else if(length(group) > 0){
      Event.annotation.split.AltAD$GROUP[i] <- paste(group, collapse = ",")
   }
   
}

table(Event.annotation.split.AltAD$GROUP)
Event.annotation.split.AltAD <- separate_longer_delim(Event.annotation.split.AltAD, cols = "GROUP", delim = ",")
table(Event.annotation.split.AltAD$GROUP)
Event.annotation.split.AltAD <- Event.annotation.split.AltAD[!(is.na(Event.annotation.split.AltAD$GENEID)),]
Event.annotation.split.AltAD <- Event.annotation.split.AltAD[!(is.na(Event.annotation.split.AltAD$GROUP)),]
table(Event.annotation.split.AltAD$GROUP)

Event.annotation.split.AltAD$EVENTID <- make.unique(Event.annotation.split.AltAD$EVENTID)
Matt_AltAD <- Event.annotation.split.AltAD[,c(1:7,14)]
write.table(Matt_AltAD, "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Matt_compare/vsAltAD/Matt_AltAD.tab", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(data.frame(unique(Matt_AltAD$SCAFF)), "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Matt_compare/vsAltAD/AltAD_Seqnames.txt", sep = "\t", row.names = F, col.names = F, quote = F)
Matt_AltAD_results <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Matt_compare/vsAltAD/overview_of_features_and_comparisons.tab", header = T)
Matt_AltAD_results <- Matt_AltAD_results[,c(1,22,23, 30:35)]
Matt_AltAD_results.filtered <- Matt_AltAD_results[which(Matt_AltAD_results$PVAL_MANN_WHITNEY_U_TEST_Genome_VS_PanAS <= 0.05 | Matt_AltAD_results$PVAL_MANN_WHITNEY_U_TEST_Genome_VS_StressAS <= 0.05 | Matt_AltAD_results$PVAL_MANN_WHITNEY_U_TEST_Genome_VS_TissuesAS <= 0.05 | Matt_AltAD_results$PVAL_MANN_WHITNEY_U_TEST_ASNR_VS_Genome <= 0.05),]
View(Matt_AltAD_results.filtered)
write.table(x = Matt_AltAD_results.filtered, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/Matt_compare/vsAltAD/Matt_AltAD_results_filtered.txt", quote = F, sep = "\t", col.names = T, row.names = F)

#### 3.4 WGCNA analyses ####
## Only for expression

options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 15)
output.path <- "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/WGCNA/"
setwd(output.path)

## Data Prepare and checks

data0 <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Batch_Corrected/Final_VST_lowFilt_Presence_batchCorrected.txt", header = T, sep = "\t")
data0 <- t(data0)

sample_metada <- read.delim("clipboard", header = T)

sample_metada$Tissue <- gsub(pattern = "Phloem", replacement = "Vascular", x = sample_metada$Tissue)
sample_metada$Tissue <- gsub(pattern = "Xylem", replacement = "Vascular", x = sample_metada$Tissue)
labelitass <- list(SampleID = sample_metada$SampleNameUnique,Tissues = sample_metada$Tissue, Group = sample_metada$Group, Stress = sample_metada$Stress, Treatment =  sample_metada$Treatment, Intensity = sample_metada$Time, Genotype = sample_metada$Genotype, Technique = sample_metada$Batch1, Study = sample_metada$Batch2,
                   Treatment.Intensity = paste0(sample_metada$Treatment, "_", sample_metada$Time), Tissue.Stress = paste0(sample_metada$Tissue, "_", sample_metada$Stress), Tissues.Stress.Treatment.Intensity = paste0(sample_metada$Tissue, "_",sample_metada$Stress, "_", sample_metada$Treatment, "_", sample_metada$Time) )

labelitass$Intensity <- as.character(labelitass$Intensity)
labelitass$Study <- as.character(labelitass$Study)
sample0.metadata <- as.data.frame(do.call(cbind ,labelitass))
View(sample0.metadata)

length(match(sample0.metadata$SampleID, rownames(data0))) == nrow(sample0.metadata)

# check outliers: NO OUTLIERS; log10 best transformation

goodSamplesGenes(datExpr = data0) # ALL TRUE/OK

sampleTree = hclust(dist(data0), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
# filter Needle_Control_Control.6.8.9.10 
rownames(data0)
dim(data0)
data0 <- data0[-c(40,42,43,44),]
rownames(data0)
dim(data0)

## Choose soft threshold parameter

powers = c(c(1:100), seq(from = 22, to=100, by=2))
sft = pickSoftThreshold(data0, powerVector = powers, verbose = 5, dataIsExpr = TRUE, networkType = "signed hybrid") # blockSize = 30 

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

softPower.hybrid <- 9

## Turning data into topological overlap matrix

adjacency.hybrid = adjacency(datExpr = data0, power = softPower.hybrid, type = "signed hybrid", corFnc = "bicor")
TOM.hybrid = TOMsimilarity(adjacency.hybrid)
dissTOM.hybrid = 1 - TOM.hybrid

## Grouping genes in modules

geneTree.hybrid = hclust(as.dist(dissTOM.hybrid), method = "average")

dynamicMods.hybrid = cutreeDynamic(dendro = geneTree.hybrid, distM = dissTOM.hybrid, deepSplit = 2, 
                                   pamRespectsDendro = FALSE, minClusterSize = 30)

table(dynamicMods.hybrid)
length(table(dynamicMods.hybrid)) 

dynamicColors.hybrid = labels2colors(dynamicMods.hybrid)
table(dynamicColors.hybrid)

plotDendroAndColors(geneTree.hybrid, dynamicColors.hybrid, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")


## Merge modules

MEList.hybrid = moduleEigengenes(data0, colors = dynamicColors.hybrid)
MEs.hybrid = MEList.hybrid$eigengenes
MEDiss.hybrid = 1-cor(MEs.hybrid)
METree.hybrid = hclust(as.dist(MEDiss.hybrid), method = "average")
plot(METree.hybrid, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres.hybrid = 0.3
abline(h=MEDissThres.hybrid, col = "red")
merge.hybrid = mergeCloseModules(data0, dynamicColors.hybrid, cutHeight = MEDissThres.hybrid, verbose = 3) 
mergedColors.hybrid = merge.hybrid$colors  
mergedMEs.hybrid = merge.hybrid$newMEs  
plotDendroAndColors(geneTree.hybrid, cbind(dynamicColors.hybrid, mergedColors.hybrid), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05) 
write.table(merge.hybrid$oldMEs,file="Hybrid.oldMEs.txt", quote = F, sep = "\t")
write.table(merge.hybrid$newMEs,file="Hybrid.newMEs.txt", quote = F, sep = "\t")

moduleColors.hybrid = mergedColors.hybrid

colorOrder = c("grey", standardColors(50))

moduleLabels.hybrid = match(moduleColors.hybrid, colorOrder)-1
MEs.hybrid = mergedMEs.hybrid

## Associating modules and phenotypes

sample0.metadata$Tissues.Stress.Treatment <- paste0(sample0.metadata$Tissues, "_", sample0.metadata$Stress, "_", sample0.metadata$Treatment)
datTrait <- binarizeCategoricalColumns(sample0.metadata$Tissues.Stress.Treatment.Intensity[-c(40,42,43,44)], includePairwise = F, includeLevelVsAll = T, dropFirstLevelVsAll = F, minCount = 1)
datTrait <- binarizeCategoricalColumns(sample0.metadata$Tissues.Stress.Treatment[-c(40,42,43,44)], includePairwise = F, includeLevelVsAll = T, dropFirstLevelVsAll = F, minCount = 1)

nGenes = ncol(data0)
nSamples = nrow(data0)

MEs0.hybrid = moduleEigengenes(data0, moduleColors.hybrid)$eigengenes
MEs.hybrid = orderMEs(MEs0.hybrid)
moduleTraitCor.hybrid = cor(MEs.hybrid, datTrait, use = "p")
moduleTraitPvalue.hybrid = corPvalueStudent(moduleTraitCor.hybrid, nSamples)
BH=0.05/ncol(moduleTraitCor.hybrid)

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
rownames(moduleTraitCor.hybrid) <- c("M06(1064)", "M04(2833)", "M19(74)", "M05(1360)", "M16(234)", "M15(238)", "M07(635)", "M18(107)", "M09(566)", "M17(137)", "M12(277)",
                                     "M03(3475)", "M13(277)", "M08(598)", "M10(545)", "M20(57)", "M14(271)", "M11(349)", "M02(3683)", "M01(5879)")
colnames(MEs.hybrid) <- rownames(moduleTraitCor.hybrid)

colors <- sort(table(mergedColors.hybrid), decreasing = T)
for(i in 1:20){
   if(i < 10){
      mergedColors.hybrid[which(mergedColors.hybrid == names(colors)[i])] <- paste0("M0", i)   
   }else if(i >= 10){
      mergedColors.hybrid[which(mergedColors.hybrid == names(colors)[i])] <- paste0("M", i)
   }
}

write.table(moduleTraitCor.hybrid, file = "moduleTraitCor.hybrid.filtered.renamed.woIntensity.txt", quote = F, sep = "\t", row.names = T, col.names = T)
write.table(MEs.hybrid, file = "MEs.hybrid.filtered.renamed.woIntensity.txt", quote = F, sep = "\t", row.names = T, col.names = T)
ModuleInfo.WGCNA <- data.frame(GeneID = rownames(data0), Module = mergedColors.hybrid)
write.table(ModuleInfo.WGCNA, "ModuleInfoWGCNA.txt", quote = F, sep = "\t", row.names = F, col.names = T)

## plotting final TraitCor heatmap

col_trait <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

ht.Trait.hor <- Heatmap(t(moduleTraitCor.hybrid), name = "Trait.cor", col = col_fun, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F,
                        right_annotation = rowAnnotation(Tissues = c(rep("Bud", 8), rep("Needle", 5), rep("Vascular", 1)), 
                                                         Stress = c(rep("Control", 4), rep("Fusarium", 4), "Control", "Dothistroma", "Heat", "Heat", "Phytophthora", "Control"), 
                                                         Intensity = c("T0", "T1", "T3", "T4", "T1", "T2", "T3", "T4", "T0", "T2", "T1", "T2", "T1", "T0"),
                                                         col = list(Tissues = c("Vascular"="darkgoldenrod", "Bud"="darkmagenta", "Needle"="darkgreen" ),
                                                                    Stress = c("Control"="#00A087FF", "Fusarium"="#7E6148FF", "Phytophthora"="#3C5488FF", "Heat"="#E64B35FF", "Dothistroma"="darkgoldenrod"),
                                                                    Intensity = c("T0"="gray88", "T1"="gray68", "T2"="gray48", "T3"="gray28", "T4"="gray18", "T5"="gray8")) 
                        ),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                           if(t(moduleTraitCor.hybrid)[i, j] > 0 | t(moduleTraitCor.hybrid)[i, j] < 0)
                              grid.text(sprintf("%.1f", t(moduleTraitCor.hybrid)[i, j]), x, y, gp = gpar(fontsize = 10))
                        })
ht.Trait.ver <- Heatmap(moduleTraitCor.hybrid, name = "Trait.cor", col = col_fun, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F,
                        bottom_annotation = HeatmapAnnotation(Tissues = c(rep("Bud", 8), rep("Needle", 5), rep("Vascular", 1)), 
                                                         Stress = c(rep("Control", 4), rep("Fusarium", 4), "Control", "Dothistroma", "Heat", "Heat", "Phytophthora", "Control"), 
                                                         Intensity = c("T0", "T1", "T3", "T4", "T1", "T2", "T3", "T4", "T0", "T2", "T1", "T2", "T1", "T0"),
                                                         col = list(Tissues = c("Vascular"="darkgoldenrod", "Bud"="darkmagenta", "Needle"="darkgreen" ),
                                                                    Stress = c("Control"="#00A087FF", "Fusarium"="#7E6148FF", "Phytophthora"="#3C5488FF", "Heat"="#E64B35FF", "Dothistroma"="darkgoldenrod"),
                                                                    Intensity = c("T0"="gray88", "T1"="gray68", "T2"="gray48", "T3"="gray28", "T4"="gray18", "T5"="gray8")) 
                        ),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                           if(moduleTraitCor.hybrid[i, j] > 0 | moduleTraitCor.hybrid[i, j] < 0)
                              grid.text(sprintf("%.1f", moduleTraitCor.hybrid[i, j]), x, y, gp = gpar(fontsize = 10))
                        })

ht.Trait.hor <- Heatmap(t(moduleTraitCor.hybrid), name = "Trait.cor", col = col_fun, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F,
                        right_annotation = rowAnnotation(Tissues = c(rep("Bud", 3), rep("Needle", 4), rep("Vascular", 1)), 
                                                         Stress = c(rep("Control", 2), rep("Fusarium", 1), "Control", "Dothistroma", "Heat", "Phytophthora", "Control"), 
                                                         
                                                         col = list(Tissues = c("Vascular"="darkgoldenrod", "Bud"="darkmagenta", "Needle"="darkgreen" ),
                                                                    Stress = c("Control"="#00A087FF", "Fusarium"="#7E6148FF", "Phytophthora"="#3C5488FF", "Heat"="#E64B35FF", "Dothistroma"="darkgoldenrod")
                                                                    ) 
                        ),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                           if(t(moduleTraitCor.hybrid)[i, j] > 0 | t(moduleTraitCor.hybrid)[i, j] < 0)
                              grid.text(sprintf("%.1f", t(moduleTraitCor.hybrid)[i, j]), x, y, gp = gpar(fontsize = 10))
                        })
ht.Trait.ver <- Heatmap(moduleTraitCor.hybrid, name = "Trait.cor", col = col_fun, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F,
                        bottom_annotation = HeatmapAnnotation(Tissues = c(rep("Bud", 3), rep("Needle", 4), rep("Vascular", 1)), 
                                                              Stress = c(rep("Control", 2), rep("Fusarium", 1), "Control", "Dothistroma", "Heat", "Phytophthora", "Control"), 
                                                              
                                                              col = list(Tissues = c("Vascular"="darkgoldenrod", "Bud"="darkmagenta", "Needle"="darkgreen" ),
                                                                         Stress = c("Control"="#00A087FF", "Fusarium"="#7E6148FF", "Phytophthora"="#3C5488FF", "Heat"="#E64B35FF", "Dothistroma"="darkgoldenrod")
                                                                         ) 
                        ),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                           if(moduleTraitCor.hybrid[i, j] > 0 | moduleTraitCor.hybrid[i, j] < 0)
                              grid.text(sprintf("%.1f", moduleTraitCor.hybrid[i, j]), x, y, gp = gpar(fontsize = 10))
                        })

## Enrichment based on module membership

MEs.hybrid <- read.delim("Hybrid.newMEs.txt")
geneModuleMembership = as.data.frame(cor(data0, MEs.hybrid, use = "p"))

# create expression Bin_db

Gene.annotation <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Gene_annotation.txt")
Gen.annotation.bin <- separate_longer_delim(Gene.annotation, Mercator4.Bin, delim = ",")
Bin_db <- list()
bin.names <- levels(as.factor(Gen.annotation.bin$Mercator4.Bin))
for(i in bin.names){
   Bin_db[[i]] <- Gen.annotation.bin$GeneID[Gen.annotation.bin$Mercator4.Bin == i]
}
names(Bin_db) <- c("PS", "Redox hom", "Phytohormone act", "Chromatin org", "Cell division", "DNA dmg response", "RNA biosynthesis", "RNA processing", 
                   "Protein biosynthesis", "Protein modification", "Protein hom",
                   "Cellular respiration", "Cytoskeleton org", "Cell wall org", "Vessicle trafficking", "Protein translocation", "Solute transport",
                   "Nutrient uptake", "External stimuli resp", "Multi-process reg", "Plant reproduction", "Plant organogen",
                   "Carbohydrate met", "Clade-specific met", "Unknown", 
                   "Aminoacid met", "Lipid met", "Enzyme classification", "Nucleotide met", "Coenzyme met", "Polyamine met", "Secondary met")


# compute module membership enrichments

output.path <- "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/WGCNA/"
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


Heatmap.Modules.Q <- data.frame(Factor = Membership.res$M06.1064.$pathway, 
                                "M06(1064)" = NA, "M04(2833)" = NA, "M19(74)" = NA, "M05(1360)" = NA, "M16(234)" = NA,
                                "M15(238)" = NA, "M07(635)" = NA, "M18(107)" = NA, "M09(566)" = NA, "M17(137)" = NA,
                                "M12(277)" = NA, "M03(3475)" = NA, "M13(277)" = NA, "M08(598)" = NA, "M10(545)" = NA,
                                "M20(57)" = NA, "M14(271)" = NA, "M11(349)" = NA, "M02(3683)" = NA, "M01(5879)" = NA
                                )
Heatmap.Modules.S <- data.frame(Factor = Membership.res$M06.1064.$pathway, 
                                "M06(1064)" = NA, "M04(2833)" = NA, "M19(74)" = NA, "M05(1360)" = NA, "M16(234)" = NA,
                                "M15(238)" = NA, "M07(635)" = NA, "M18(107)" = NA, "M09(566)" = NA, "M17(137)" = NA,
                                "M12(277)" = NA, "M03(3475)" = NA, "M13(277)" = NA, "M08(598)" = NA, "M10(545)" = NA,
                                "M20(57)" = NA, "M14(271)" = NA, "M11(349)" = NA, "M02(3683)" = NA, "M01(5879)" = NA)

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

col_significance <- colorRamp2(c(0,1), c("white", "dodgerblue4"))
col_quantity <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

ht.Q.Modules <- Heatmap(Heatmap.Modules.Q[,2:ncol(Heatmap.Modules.Q)], name = "Modules.NES", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F,
                        top_annotation = HeatmapAnnotation(Type = c(rep("Quantity", 20)), col=list(Type = c("Quantity" = "#374E55FF"))))
ht.S.Modules <- Heatmap(Heatmap.Modules.S[,2:ncol(Heatmap.Modules.S)], name = "Modules.S", col = col_significance, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F,
                        top_annotation = HeatmapAnnotation(Type = c(rep("Significance", 20)), col=list(Type = c("Significance" = "#DF8F44FF"))))

ht.Q.Modules + ht.S.Modules

#### 3.5 Functional Enrichment sets ORA ####

## Splicing

# create splicing Bin_db

Event.annotation.bin <- separate_longer_delim(Event.annotation.split, Mercator4.Bin, delim = ",")
View(Event.annotation.bin)
Event.annotation.bin <- Event.annotation.bin[!is.na(Event.annotation.bin$Mercator4.Bin),]
Bin_db_splicing <- list()
bin.names.splicing <- levels(as.factor(Event.annotation.bin$Mercator4.Bin))
for(i in bin.names.splicing){
   Bin_db_splicing[[i]] <- Event.annotation.bin$GENEID[Event.annotation.bin$Mercator4.Bin == i]
}
names(Bin_db_splicing) <- c("PS", "Redox hom", "Phytohormone act", "Chromatin org", "Cell division", "DNA dmg response", "RNA biosynthesis", "RNA processing", 
                   "Protein biosynthesis", "Protein modification", "Protein hom",
                   "Cellular respiration", "Cytoskeleton org", "Cell wall org", "Vessicle trafficking", "Protein translocation", "Solute transport",
                   "Nutrient uptake", "External stimuli resp", "Multi-process reg", "Plant reproduction", "Plant organogen",
                   "Carbohydrate met", "Clade-specific met", "Unknown", 
                   "Aminoacid met", "Lipid met", "Enzyme classification", "Nucleotide met", "Coenzyme met", "Polyamine met", "Secondary met")


# Compute enrichments

Splicing <- list(PanAS = unique(Event.annotation.bin$GENEID[Event.annotation.bin$EVENTID %in% PanAS$x]), 
                 TissuesAS = unique(Event.annotation.bin$GENEID[Event.annotation.bin$EVENTID %in% TissuesAS$x]), 
                 StressAS = unique(Event.annotation.bin$GENEID[Event.annotation.bin$EVENTID %in% StressAS$x]))
Universe <- unique(unlist(Bin_db_splicing, use.names = F))
Splicing <- lapply(Splicing, function(x){
   x[x %in% Universe]
})
Splicing.res <- list()
for(i in names(Splicing)){
   data <- Splicing[[i]]
   #data$ID <- rownames(data)
   #data <- data[,c("ID", "logFC")]
   #data <- deframe(data)
   fGseaRES <- fora(pathways = Bin_db_splicing, genes = data, universe = Universe, minSize = 5, maxSize = 1000)
   write.table(fGseaRES[,-ncol(fGseaRES)], file = paste0(output.path, i, "_fora.txt"), col.names = T, sep = "\t", quote = F)
   Splicing.res[[i]] <- fGseaRES
}

# Heatmap plotting Purples

Splicing.res <- lapply(Splicing.res, function(x){
   rownames(x) <- x$pathway
   x
})

Splicing.bin.df <- data.frame(Factor = Splicing.res$PanAS$pathway,
                              PanAS = NA,
                              TissuesAS = NA,
                              StressAS = NA)

Splicing.res$PanAS
Splicing.res$TissuesAS
Splicing.res$StressAS

for(i in 1:nrow(Splicing.bin.df)){
   for(j in 2:ncol(Splicing.bin.df)){
      hit <- which(Splicing.res[[j-1]]$pathway == Splicing.bin.df$Factor[i])
      if(length(hit) == 0){
         Splicing.bin.df[i,j] <- 1
      }else if(length(hit) > 1){
         stop("check check")
      }else if(length(hit) == 1){
         Splicing.bin.df[i,j]  <- Splicing.res[[j-1]]$padj[hit]
         if(Splicing.bin.df[i,j]  > 0.1){
            Splicing.bin.df[i,j]  <- 1
         }
      }
   }
}

rownames(Splicing.bin.df) <- Splicing.bin.df$Factor
Splicing.bin.df <- Splicing.bin.df[,-1]
Splicing.bin.df <- -log10(Splicing.bin.df)

for(i in 1:nrow(Splicing.bin.df)){
   for(j in 1:ncol(Splicing.bin.df)){
      if(Splicing.bin.df[i,j] > 5 ){
         Splicing.bin.df[i,j] <- 5
      }
   }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Purples"))(500)

ht.Bin.Splicing <- Heatmap(Splicing.bin.df, name = "Bins.-log10(FDR).Splicing", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

## Expression

# Compute enrichments

PanGE<- read.table("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Batch_Corrected/PanGE_GeneID.txt")
TissuesGE<- read.table("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Batch_Corrected/TissuesGE_GenesID.txt")
StressGE<- read.table("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Batch_Corrected/StressGE_GenesID.txt")

Expression <- list(PanGE = unique(PanGE$V1[PanGE$V1 %in% unlist(Bin_db)]), 
                 TissuesGE = unique(TissuesGE$V1[TissuesGE$V1 %in% unlist(Bin_db)]), 
                 StressGE = unique(StressGE$V1[StressGE$V1 %in% unlist(Bin_db)]))
Universe <- unique(unlist(Bin_db, use.names = F))
Expression.res <- list()
for(i in names(Expression)){
   data <- Expression[[i]]
   #data$ID <- rownames(data)
   #data <- data[,c("ID", "logFC")]
   #data <- deframe(data)
   fGseaRES <- fora(pathways = Bin_db, genes = data, universe = Universe, minSize = 5, maxSize = 1000)
   write.table(fGseaRES[,-ncol(fGseaRES)], file = paste0(output.path, i, "_fora.txt"), col.names = T, sep = "\t", quote = F)
   Expression.res[[i]] <- fGseaRES
}

# Heatmap plotting Greens

Expression.bin.df <- data.frame(Factor = Splicing.res$PanAS$pathway,
                              PanGE = NA,
                              TissuesGE = NA,
                              StressGE = NA)

for(i in 1:nrow(Expression.bin.df)){
   for(j in 2:ncol(Expression.bin.df)){
      hit <- which(Expression.res[[j-1]]$pathway == Expression.bin.df$Factor[i])
      if(length(hit) == 0){
         Expression.bin.df[i,j] <- 1
      }else if(length(hit) > 1){
         stop("check check")
      }else if(length(hit) == 1){
         Expression.bin.df[i,j]  <- Expression.res[[j-1]]$padj[hit]
         if(Expression.bin.df[i,j]  > 0.1){
            Expression.bin.df[i,j]  <- 1
         }
      }
   }
}

rownames(Expression.bin.df) <- Expression.bin.df$Factor
Expression.bin.df <- Expression.bin.df[,-1]
Expression.bin.df <- -log10(Expression.bin.df)

for(i in 1:nrow(Expression.bin.df)){
   for(j in 1:ncol(Expression.bin.df)){
      if(Expression.bin.df[i,j] > 5 ){
         Expression.bin.df[i,j] <- 5
      }
   }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Greens"))(500)

ht.Bin.Expression <- Heatmap(Expression.bin.df, name = "Bins.-log10(FDR).Expression", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

ht.Bin.Splicing + ht.Bin.Expression

## Re do modules heatmap same order for plotting and with ORA for comparisons

Modules <- list()
for(i in levels(as.factor(ModuleInfo.WGCNA$Module))){
   Modules[[i]] <- ModuleInfo.WGCNA$GeneID[ModuleInfo.WGCNA$Module == i]
}

Universe <- unique(unlist(Bin_db, use.names = F))
Modules.res <- list()

for(i in names(Modules)){
   data <- Modules[[i]]
   #data$ID <- rownames(data)
   #data <- data[,c("ID", "logFC")]
   #data <- deframe(data)
   fGseaRES <- fora(pathways = Bin_db, genes = data, universe = Universe, minSize = 5, maxSize = 1000)
   write.table(fGseaRES[,-ncol(fGseaRES)], file = paste0(output.path, i, "_fora.txt"), col.names = T, sep = "\t", quote = F)
   Modules.res[[i]] <- fGseaRES
}


Heatmap.Modules <- data.frame(Factor = Splicing.res$PanAS$pathway, 
                                "M06" = NA, "M04" = NA, "M19" = NA, "M05" = NA, "M16" = NA,
                                "M15" = NA, "M07" = NA, "M18" = NA, "M09" = NA, "M17" = NA,
                                "M12" = NA, "M03" = NA, "M13" = NA, "M08" = NA, "M10" = NA,
                                "M20" = NA, "M14" = NA, "M11" = NA, "M02" = NA, "M01" = NA)

for(i in 1:nrow(Heatmap.Modules)){
   for(j in 2:ncol(Heatmap.Modules)){
      
      colon <- which(names(Modules.res) == colnames(Heatmap.Modules)[j])
      hit <- which(Modules.res[[colon]]$pathway == Heatmap.Modules$Factor[i])
      
      if(length(hit) == 0){
         Heatmap.Modules[i,j] <- 1
      }else if(length(hit) > 1){
         stop("check check")
      }else if(length(hit) == 1){
         Heatmap.Modules[i,j]  <- Modules.res[[colon]]$padj[hit]
         if(Heatmap.Modules[i,j]  > 0.1){
            Heatmap.Modules[i,j]  <- 1
         }
      }
   }
}

rownames(Heatmap.Modules) <- Heatmap.Modules$Factor
Heatmap.Modules <- Heatmap.Modules[,-1]
Heatmap.Modules <- -log10(Heatmap.Modules)

for(i in 1:nrow(Heatmap.Modules)){
   for(j in 1:ncol(Heatmap.Modules)){
      if(Heatmap.Modules[i,j] > 5 ){
         Heatmap.Modules[i,j] <- 5
      }
   }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Greens"))(500)

ht.Bin.Modules <- Heatmap(Heatmap.Modules, name = "Bins.-log10(FDR).Modules", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)


ht.Bin.Splicing + ht.Bin.Expression + ht.Bin.Modules

#### 4. Integration ####

#### 4.0 Tissues ####

#### Evolutionary Transcriptomics ####

setwd("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Descriptive/EvolutionaryTranscriptomics/")
Tissues.TAI <- as.data.frame(t(data0)[,TissuesSamples])
Tissues.TAI$GeneID <- rownames(Tissues.TAI) 

Tissues.TAI.df <- data.frame(Phylostratum = NA,
                             GeneID = Tissues.TAI$GeneID,
                             Bud = rowMeans(Tissues.TAI[,grep("Bud", colnames(Tissues.TAI))]),
                             Needle = rowMeans(Tissues.TAI[,grep("Needle", colnames(Tissues.TAI))]),
                             Xylem = rowMeans(Tissues.TAI[,grep("Xylem", colnames(Tissues.TAI))]),
                             Phloem = rowMeans(Tissues.TAI[,grep("Phloem", colnames(Tissues.TAI))])
                             )

View(Tissues.TAI.df)
for(i in 1:nrow(Tissues.TAI.df)){
   Tissues.TAI.df$Phylostratum[i] <- Gene.annotation$PhyloStratum[which(Gene.annotation$GeneID == Tissues.TAI.df$GeneID[i])]
}

Tissues.TAI.df <- Tissues.TAI.df[!is.na(Tissues.TAI.df$Phylostratum),]
PlotSignature(Tissues.TAI.df)

#### MOFA ####

Expression <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Batch_Corrected/Final_VST_lowFilt_Presence_batchCorrected.txt")
Splicing <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/Final_PSI_lowFilter_PEandTE_batchCorrectedTechStudy_Normalized_PresenceRange.txt")
rownames(Splicing) <- Splicing$events.names
Splicing <- Splicing[,-1]

### 4.0.1 Split Splicing Events by AS type and general AS

Splicing.IR <- Splicing[which(rownames(Splicing) %in% Event.annotation.split$EVENTID[Event.annotation.split$Type == "IR"]),]
dim(Splicing.IR)
Splicing.ES <- Splicing[which(rownames(Splicing) %in% Event.annotation.split$EVENTID[Event.annotation.split$Type == "ES"]),]
dim(Splicing.ES)
Splicing.AltAD <- Splicing[which(rownames(Splicing) %in% Event.annotation.split$EVENTID[Event.annotation.split$Type == "AltAD"]),]
dim(Splicing.AltAD)
Splicing.AS <- Splicing[which(rownames(Splicing) %in% Event.annotation.split$EVENTID[Event.annotation.split$Type == "AS"]),]
dim(Splicing.AS)

AllTranscriptionalList.split <- list(Expression = as.matrix(Expression), IR = as.matrix(Splicing.IR),
                                     ES = as.matrix(Splicing.ES), AltAD = as.matrix(Splicing.AltAD),
                                     AS = as.matrix(Splicing.AS))

AllTranscriptionalList <- list(Expression = as.matrix(Expression), Splicing = as.matrix(Splicing))

AllTranscriptionalList.split.tissues <- lapply(AllTranscriptionalList.split, function(x){
   x[,TissuesSamples]
})

AllTranscriptionalList.tissues <- lapply(AllTranscriptionalList, function(x){
   x[,TissuesSamples]
})

### define groups

groups <- c(rep("Bud",16), rep("Needle",27), rep("vascular",8))

### filter 0 variance features

## Splitted

AllTranscriptionalList.split.tissues.Expression.vars <- apply(AllTranscriptionalList.split.tissues$Expression, 1, var)
AllTranscriptionalList.split.tissues.IR.vars <- apply(AllTranscriptionalList.split.tissues$IR, 1, function(x){
   var(x, na.rm = TRUE)
})
AllTranscriptionalList.split.tissues.ES.vars <- apply(AllTranscriptionalList.split.tissues$ES, 1, function(x){
   var(x, na.rm = TRUE)
})
AllTranscriptionalList.split.tissues.AltAD.vars <- apply(AllTranscriptionalList.split.tissues$AltAD, 1, function(x){
   var(x, na.rm = TRUE)
})
AllTranscriptionalList.split.tissues.AS.vars <- apply(AllTranscriptionalList.split.tissues$AS, 1, function(x){
   var(x, na.rm = TRUE)
})

names(AllTranscriptionalList.split.tissues.Expression.vars)[AllTranscriptionalList.split.Expression.vars > 0]
names(AllTranscriptionalList.split.tissues.IR.vars)[AllTranscriptionalList.split.IR.vars > 0]
names(AllTranscriptionalList.split.tissues.ES.vars)[AllTranscriptionalList.split.ES.vars > 0]
names(AllTranscriptionalList.split.tissues.AltAD.vars)[AllTranscriptionalList.split.AltAD.vars > 0]
names(AllTranscriptionalList.split.tissues.AS.vars)[AllTranscriptionalList.split.AS.vars > 0]

## Top 5,000 AS 10,000 Genes

Splitted.Genes <- names(AllTranscriptionalList.split.tissues.Expression.vars)[order(AllTranscriptionalList.split.tissues.Expression.vars, decreasing = T)[1:10000]]
Splitted.IR <- names(AllTranscriptionalList.split.tissues.IR.vars)[order(AllTranscriptionalList.split.tissues.IR.vars, decreasing = T)[1:5000]]
Splitted.ES <- names(AllTranscriptionalList.split.tissues.ES.vars)
Splitted.AltAD <- names(AllTranscriptionalList.split.tissues.AltAD.vars)[order(AllTranscriptionalList.split.tissues.AltAD.vars, decreasing = T)[1:5000]]
Splitted.AS <- names(AllTranscriptionalList.split.tissues.AS.vars)[order(AllTranscriptionalList.split.tissues.AS.vars, decreasing = T)[1:5000]]

AllTranscriptionalList.split.tissues.top <- AllTranscriptionalList.split.tissues

AllTranscriptionalList.split.tissues.top$Expression <- AllTranscriptionalList.split.tissues$Expression[rownames(AllTranscriptionalList.split.tissues$Expression) %in% Splitted.Genes, ]
AllTranscriptionalList.split.tissues.top$IR <- AllTranscriptionalList.split.tissues$IR[rownames(AllTranscriptionalList.split.tissues$IR) %in% Splitted.IR, ]
AllTranscriptionalList.split.tissues.top$ES <- AllTranscriptionalList.split.tissues$ES[rownames(AllTranscriptionalList.split.tissues$ES) %in% Splitted.ES, ]
AllTranscriptionalList.split.tissues.top$AltAD <- AllTranscriptionalList.split.tissues$AltAD[rownames(AllTranscriptionalList.split.tissues$AltAD) %in% Splitted.AltAD, ]
AllTranscriptionalList.split.tissues.top$AS <- AllTranscriptionalList.split.tissues$AS[rownames(AllTranscriptionalList.split.tissues$AS) %in% Splitted.AS, ]

lapply(AllTranscriptionalList.split.tissues, dim)
lapply(AllTranscriptionalList.split.tissues.top, dim)

rownames(AllTranscriptionalList.split.tissues$IR) <- paste0(rownames(AllTranscriptionalList.split.tissues$IR), "_IR")
rownames(AllTranscriptionalList.split.tissues$ES) <- paste0(rownames(AllTranscriptionalList.split.tissues$ES), "_ES")
rownames(AllTranscriptionalList.split.tissues$AltAD) <- paste0(rownames(AllTranscriptionalList.split.tissues$AltAD), "_AltAD")
rownames(AllTranscriptionalList.split.tissues$AS) <- paste0(rownames(AllTranscriptionalList.split.tissues$AS), "_AS")

write.table(AllTranscriptionalList.split.tissues$Expression, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/MOFA_tissues_split_all/Expression.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.tissues$IR, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/MOFA_tissues_split_all/IR.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.tissues$ES, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/MOFA_tissues_split_all/ES.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.tissues$AltAD, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/MOFA_tissues_split_all/AltAD.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.tissues$AS, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/MOFA_tissues_split_all/AS.txt", sep = "\t", col.names = T, row.names = T, quote = F)

rownames(AllTranscriptionalList.split.tissues.top$IR) <- paste0(rownames(AllTranscriptionalList.split.tissues.top$IR), "_IR")
rownames(AllTranscriptionalList.split.tissues.top$ES) <- paste0(rownames(AllTranscriptionalList.split.tissues.top$ES), "_ES")
rownames(AllTranscriptionalList.split.tissues.top$AltAD) <- paste0(rownames(AllTranscriptionalList.split.tissues.top$AltAD), "_AltAD")
rownames(AllTranscriptionalList.split.tissues.top$AS) <- paste0(rownames(AllTranscriptionalList.split.tissues.top$AS), "_AS")

rownames(AllTranscriptionalList.split.tissues.top$IR) <- gsub("_IR", "", rownames(AllTranscriptionalList.split.tissues.top$IR))
rownames(AllTranscriptionalList.split.tissues.top$ES) <- gsub("_ES", "", rownames(AllTranscriptionalList.split.tissues.top$ES))
rownames(AllTranscriptionalList.split.tissues.top$AltAD) <- gsub("_AltAD", "", rownames(AllTranscriptionalList.split.tissues.top$AltAD))
rownames(AllTranscriptionalList.split.tissues.top$AS) <- gsub("_AS", "", rownames(AllTranscriptionalList.split.tissues.top$AS))

write.table(AllTranscriptionalList.split.tissues.top$Expression, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/MOFA_tissues_split_top/Expression.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.tissues.top$IR, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/MOFA_tissues_split_top/IR.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.tissues.top$ES, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/MOFA_tissues_split_top/ES.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.tissues.top$AltAD, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/MOFA_tissues_split_top/AltAD.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.tissues.top$AS, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/MOFA_tissues_split_top/AS.txt", sep = "\t", col.names = T, row.names = T, quote = F)

## future think about HVF for each group and regulation level

## General

AllTranscriptionalList.Expression.vars <- apply(AllTranscriptionalList.tissues$Expression, 1, var)
AllTranscriptionalList.Splicing.vars <- apply(AllTranscriptionalList.tissues$Splicing, 1, function(x){
   var(x, na.rm = TRUE)
})

## Top 10,000 both

General.Genes <- names(AllTranscriptionalList.Expression.vars)[order(AllTranscriptionalList.Expression.vars, decreasing = T)[1:10000]]
General.Splicing <- names(AllTranscriptionalList.Splicing.vars)[order(AllTranscriptionalList.Splicing.vars, decreasing = T)[1:10000]]

AllTranscriptionalList.tissues.top <- AllTranscriptionalList.tissues

AllTranscriptionalList.tissues.top$Expression <- AllTranscriptionalList.tissues$Expression[rownames(AllTranscriptionalList.tissues$Expression) %in% General.Genes,]
AllTranscriptionalList.tissues.top$Splicing <- AllTranscriptionalList.tissues$Splicing[rownames(AllTranscriptionalList.tissues$Splicing) %in% General.Splicing,]

lapply(AllTranscriptionalList.tissues, dim)
lapply(AllTranscriptionalList.tissues.top, dim)

write.table(AllTranscriptionalList.tissues$Expression, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/MOFA_tissues_all_all/Expression.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.tissues$Splicing, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/MOFA_tissues_all_all/Splicing.txt", sep = "\t", col.names = T, row.names = T, quote = F)

write.table(AllTranscriptionalList.tissues.top$Expression, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/MOFA_tissues_all_top/Expression.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.tissues.top$Splicing, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/MOFA_tissues_all_top/Splicing.txt", sep = "\t", col.names = T, row.names = T, quote = F)


## future think about HVF for each group and regulation level

AllTranscriptionalList.split.tissues <- list(Expression = as.matrix(read.table("/data/Transcriptional_MOFA/MOFA_tissues_split_all/Expression.txt")),
                                             IR = as.matrix(read.table("/data/Transcriptional_MOFA/MOFA_tissues_split_all/IR.txt")),
                                             ES = as.matrix(read.table("/data/Transcriptional_MOFA/MOFA_tissues_split_all/ES.txt")),
                                             AltAD = as.matrix(read.table("/data/Transcriptional_MOFA/MOFA_tissues_split_all/AltAD.txt")),
                                             AS = as.matrix(read.table("/data/Transcriptional_MOFA/MOFA_tissues_split_all/AS.txt"))
                                             )

AllTranscriptionalList.split.tissues.top <- list(Expression = as.matrix(read.table("/data/Transcriptional_MOFA/MOFA_tissues_split_top/Expression.txt")),
                                             IR = as.matrix(read.table("/data/Transcriptional_MOFA/MOFA_tissues_split_top/IR.txt")),
                                             ES = as.matrix(read.table("/data/Transcriptional_MOFA/MOFA_tissues_split_top/ES.txt")),
                                             AltAD = as.matrix(read.table("/data/Transcriptional_MOFA/MOFA_tissues_split_top/AltAD.txt")),
                                             AS = as.matrix(read.table("/data/Transcriptional_MOFA/MOFA_tissues_split_top/AS.txt"))
)

AllTranscriptionalList.tissues <- list(Expression = as.matrix(read.table("/data/Transcriptional_MOFA/MOFA_tissues_all_all/Expression.txt")),
                                             Splicing = as.matrix(read.table("/data/Transcriptional_MOFA/MOFA_tissues_all_all/Splicing.txt"))
)

AllTranscriptionalList.tissues.top <- list(Expression = as.matrix(read.table("/data/Transcriptional_MOFA/MOFA_tissues_all_top/Expression.txt")),
                                                 Splicing = as.matrix(read.table("/data/Transcriptional_MOFA/MOFA_tissues_all_top/Splicing.txt"))
)

MOFAgrouped <- create_mofa(data = AllTranscriptionalList.split.tissues, save_metadata = TRUE)

## 4.0.2 Prepare MOFA and run

plot_data_overview(MOFAgrouped)
sample_metadata <- data.frame(
   sample = colnames(MOFAgrouped@data$Expression$group1),
   tissues = groups <- c(rep("Bud",16), rep("Needle",27), rep("Phloem",2), rep("Xylem",6)),
   group = groups,
   sample_number = c(1:51)
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
outfileGrouped <- "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/"

MOFAgrouped <- prepare_mofa(MOFAgrouped, data_options = data_opts, model_options = model_opts, training_options = train_opts)
MOFAgrouped.trained <- run_mofa(MOFAgrouped, outfile = outfileGrouped)

MOFAgrouped <- prepare_mofa(MOFAgrouped, data_options = data_opts, model_options = model_opts, training_options = train_opts)
MOFAgrouped.trained <- run_mofa(MOFAgrouped, outfile = outfileGrouped)

MOFAgrouped <- load_model("F:/ATLAS_Pra/Protein_module/0.Process/Transcriptional_MOFA/MOFA_tissues_all_top/MOFA_tissues_all_top_400_false.hdf5")
MOFAgrouped.trained <- MOFAgrouped

## 4.0.3 Plot General

plot_factor_cor(MOFAgrouped.trained)
plot_variance_explained(MOFAgrouped.trained, x = "group", y = "factor")
plot_variance_explained(MOFAgrouped.trained, x = "group", y = "factor", plot_total = T)

MOFAgrouped.trained@cache$variance_explained
MOFAgrouped.trained@cache$variance_explained$r2_per_factor$Bud
MOFAgrouped.trained@cache$variance_explained$r2_per_factor$Needle
MOFAgrouped.trained@cache$variance_explained$r2_per_factor$vascular

a <- plot_variance_explained(MOFAgrouped.trained, x="view", y="factor", factors = c(1:8))
b <- plot_variance_explained(MOFAgrouped.trained, x="group", y="factor", plot_total = T, factors = c(1:8))[[2]]
PlotVariance <- grid.arrange(b, a, ncol = 1, nrow = 2)

## 4.0.4 LFs Inspection

# Single

MOFAgrouped.trained@samples_metadata <- sample_metadata

plot_factor(MOFAgrouped.trained, 
            factor = c(1:8),
            color_by = "tissues",
            shape_by = "tissues", dot_size = 4
)

# Multiple

plot_factors(MOFAgrouped.trained, 
             factors = c(2,3),
             color_by = "tissues", dot_size = 3, shape_by = "tissues"
)

plot_factors(MOFAgrouped.trained, 
             factors = c(1,2,3),
             color_by = "stress", dot_size = 3, shape_by = "sex"
)

# Weights of biologically relevant LFs

LF1 <- plot_top_weights(MOFAgrouped.trained,
                        view = c("Expression", "Splicing"),
                        factor = 1,
                        nfeatures = 20,     # Number of features to highlight
                        scale = T,          # Scale weights from -1 to 1
                        abs = F  ) 

LF3 <- plot_top_weights(MOFAgrouped.trained,
                        view = c("Expression", "Splicing"),
                        factor = 3,
                        nfeatures = 20,     # Number of features to highlight
                        scale = T,          # Scale weights from -1 to 1
                        abs = F  ) 

gridExtra::grid.arrange(LF1, LF3, ncol = 2, nrow = 1)

# Enrichment Analyses #

# Mercator

Bin_db

Bin_db_splicing <- list()
bin.names.splicing <- levels(as.factor(Event.annotation.bin$Mercator4.Bin))
for(i in bin.names.splicing){
   Bin_db_splicing[[i]] <- Event.annotation.bin$EVENTID[Event.annotation.bin$Mercator4.Bin == i]
}
names(Bin_db_splicing) <- c("PS", "Redox hom", "Phytohormone act", "Chromatin org", "Cell division", "DNA dmg response", "RNA biosynthesis", "RNA processing", 
                            "Protein biosynthesis", "Protein modification", "Protein hom",
                            "Cellular respiration", "Cytoskeleton org", "Cell wall org", "Vessicle trafficking", "Protein translocation", "Solute transport",
                            "Nutrient uptake", "External stimuli resp", "Multi-process reg", "Plant reproduction", "Plant organogen",
                            "Carbohydrate met", "Clade-specific met", "Unknown", 
                            "Aminoacid met", "Lipid met", "Enzyme classification", "Nucleotide met", "Coenzyme met", "Polyamine met", "Secondary met")


# Gene Age

Gene.annotation <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Gene_annotation.txt", header = T)
Gene.annotation.age <- Gene.annotation[!is.na(Gene.annotation$PhyloStratum),]

table(Gene.annotation.age$PhyloStratum)

PS_db <- list()
PS_names <- levels(as.factor(Gene.annotation.age$PhyloStratum))

for(i in PS_names){
   PS_db[[i]] <- Gene.annotation.age$GeneID[Gene.annotation.age$PhyloStratum == i]
}

names(PS_db) <- paste0("PS", names(PS_db))
PS_db

PS_db_splicing <- list()
PS_names_splicing <- levels(as.factor(Event.annotation.bin$Phylostratum))
for(i in PS_names_splicing){
   PS_db_splicing[[i]] <- unique(Event.annotation.bin$EVENTID[Event.annotation.bin$Phylostratum == i])
}

names(PS_db_splicing) <- paste0("PS", names(PS_db_splicing))
PS_db_splicing

# Gene family founder

table(Gene.annotation$GeneFamilyFounder.Age)

Gene.annotation.family <- Gene.annotation
Gene.annotation.family <- Gene.annotation.family[!is.na(Gene.annotation.family$GeneFamilyFounder.Age),]

PSF_db <- list()
PSF_names <- levels(as.factor(Gene.annotation.family$GeneFamilyFounder.Age))

for(i in PSF_names){
   PSF_db[[i]] <- Gene.annotation.family$GeneID[Gene.annotation.family$GeneFamilyFounder.Age == i]
}

names(PSF_db) <- paste0("PS", names(PSF_db))
PSF_db

PSF_db_splicing <- list()
PSF_names_splicing <- levels(as.factor(Event.annotation.bin$GeneFamilyFounder.Age))

for(i in PSF_names_splicing){
   PSF_db_splicing[[i]] <- unique(Event.annotation.bin$EVENTID[Event.annotation.bin$GeneFamilyFounder.Age == i])
}

names(PSF_db_splicing) <- paste0("PS", names(PSF_db_splicing))
PSF_db_splicing

# Compute enrichments; LF1-4; Pos and neg; Mercator, Gene Age and Family foundation

Bin.paths <- list_to_matrix(Bin_db)
PS.paths <- list_to_matrix(PS_db)
PSF.paths <- list_to_matrix(PSF_db)

Bin_splicing.paths <- list_to_matrix(Bin_db_splicing)
PS_splicing.paths <- list_to_matrix(PS_db_splicing)
PSF_splicing.paths <- list_to_matrix(PSF_db_splicing)

# Expression

enrichment.parametric.Bin.all.Expression <- run_enrichment(MOFAgrouped.trained,
                                                view = "Expression", factors = c(1,3),
                                                feature.sets = t(Bin.paths),
                                                sign = "all", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, p.adj.method = "fdr", alpha = 0.05
                                                #nperm = 5000
)


enrichment.parametric.PS.all.Expression <- run_enrichment(MOFAgrouped.trained,
                                               view = "Expression", factors = c(1,3),
                                               feature.sets = t(PS.paths),
                                               sign = "all", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                               # statistical.test = "permutation",
                                               min.size = 5, alpha = 0.05
                                               #nperm = 5000
)

enrichment.parametric.PSF.all.Expression <- run_enrichment(MOFAgrouped.trained,
                                                view = "Expression", factors = c(1,3),
                                                feature.sets = t(PSF.paths),
                                                sign = "all", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                # statistical.test = "permutation",
                                                min.size = 5, alpha = 0.05
                                                #nperm = 5000
)

# IR

enrichment.parametric.Bin.all.IR <- run_enrichment(MOFAgrouped.trained,
                                                           view = "IR", factors = c(1,5),
                                                           feature.sets = t(Bin_splicing.paths),
                                                           sign = "all", set.statistic = "rank.sum",
                                                           # statistical.test = "permutation",
                                                           min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                           #nperm = 5000
)


enrichment.parametric.PS.all.IR <- run_enrichment(MOFAgrouped.trained,
                                                          view = "IR", factors = c(1,5),
                                                          feature.sets = t(PS_splicing.paths),
                                                          sign = "all", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                                          # statistical.test = "permutation",
                                                          min.size = 5, alpha = 0.1
                                                          #nperm = 5000
)

enrichment.parametric.PSF.all.IR <- run_enrichment(MOFAgrouped.trained,
                                                           view = "IR", factors = c(1,5),
                                                           feature.sets = t(PSF_splicing.paths),
                                                           sign = "all", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                           # statistical.test = "permutation",
                                                           min.size = 5, alpha = 0.1
                                                           #nperm = 5000
)

# ES

enrichment.parametric.Bin.all.ES <- run_enrichment(MOFAgrouped.trained,
                                                   view = "ES", factors = c(1,5),
                                                   feature.sets = t(Bin_splicing.paths),
                                                   sign = "all", set.statistic = "rank.sum",
                                                   # statistical.test = "permutation",
                                                   min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                   #nperm = 5000
)


enrichment.parametric.PS.all.ES <- run_enrichment(MOFAgrouped.trained,
                                                  view = "ES", factors = c(1,5),
                                                  feature.sets = t(PS_splicing.paths),
                                                  sign = "all", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                                  # statistical.test = "permutation",
                                                  min.size = 5, alpha = 0.1
                                                  #nperm = 5000
)

enrichment.parametric.PSF.all.ES <- run_enrichment(MOFAgrouped.trained,
                                                   view = "ES", factors = c(1,5),
                                                   feature.sets = t(PSF_splicing.paths),
                                                   sign = "all", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                   # statistical.test = "permutation",
                                                   min.size = 5, alpha = 0.1
                                                   #nperm = 5000
)

# AltAD

enrichment.parametric.Bin.all.AltAD <- run_enrichment(MOFAgrouped.trained,
                                                   view = "AltAD", factors = c(1,5),
                                                   feature.sets = t(Bin_splicing.paths),
                                                   sign = "all", set.statistic = "rank.sum",
                                                   # statistical.test = "permutation",
                                                   min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                   #nperm = 5000
)


enrichment.parametric.PS.all.AltAD <- run_enrichment(MOFAgrouped.trained,
                                                  view = "AltAD", factors = c(1,5),
                                                  feature.sets = t(PS_splicing.paths),
                                                  sign = "all", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                                  # statistical.test = "permutation",
                                                  min.size = 5, alpha = 0.1
                                                  #nperm = 5000
)

enrichment.parametric.PSF.all.AltAD <- run_enrichment(MOFAgrouped.trained,
                                                   view = "AltAD", factors = c(1,5),
                                                   feature.sets = t(PSF_splicing.paths),
                                                   sign = "all", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                   # statistical.test = "permutation",
                                                   min.size = 5, alpha = 0.1
                                                   #nperm = 5000
)

# AS

enrichment.parametric.Bin.all.AS <- run_enrichment(MOFAgrouped.trained,
                                                      view = "AS", factors = c(1,5),
                                                      feature.sets = t(Bin_splicing.paths),
                                                      sign = "all", set.statistic = "rank.sum",
                                                      # statistical.test = "permutation",
                                                      min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                      #nperm = 5000
)


enrichment.parametric.PS.all.AS <- run_enrichment(MOFAgrouped.trained,
                                                     view = "AS", factors = c(1,5),
                                                     feature.sets = t(PS_splicing.paths),
                                                     sign = "all", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                                     # statistical.test = "permutation",
                                                     min.size = 5, alpha = 0.1
                                                     #nperm = 5000
)

enrichment.parametric.PSF.all.AS <- run_enrichment(MOFAgrouped.trained,
                                                      view = "AS", factors = c(1,5),
                                                      feature.sets = t(PSF_splicing.paths),
                                                      sign = "all", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                      # statistical.test = "permutation",
                                                      min.size = 5, alpha = 0.1
                                                      #nperm = 5000
)

# Splicing

enrichment.parametric.Bin.all.Splicing <- run_enrichment(MOFAgrouped.trained,
                                                   view = "Splicing", factors = c(1,3),
                                                   feature.sets = t(Bin_splicing.paths),
                                                   sign = "all", set.statistic = "rank.sum",
                                                   # statistical.test = "permutation",
                                                   min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                   #nperm = 5000
)


enrichment.parametric.PS.all.Splicing <- run_enrichment(MOFAgrouped.trained,
                                                  view = "Splicing", factors = c(1,3),
                                                  feature.sets = t(PS_splicing.paths),
                                                  sign = "all", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                                  # statistical.test = "permutation",
                                                  min.size = 5, alpha = 0.1
                                                  #nperm = 5000
)

enrichment.parametric.PSF.all.Splicing <- run_enrichment(MOFAgrouped.trained,
                                                   view = "Splicing", factors = c(1,3),
                                                   feature.sets = t(PSF_splicing.paths),
                                                   sign = "all", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                   # statistical.test = "permutation",
                                                   min.size = 5, alpha = 0.1
                                                   #nperm = 5000
)

## plot enrich

# expression

order.bin.exp <- sort(rownames(enrichment.parametric.Bin.all.Expression$pval.adj))

Bin.df.exp <- t(enrichment.parametric.Bin.all.Expression$pval.adj)[,order.bin.exp]

for(i in 1:nrow(Bin.df.exp)){
   for(j in 1:ncol(Bin.df.exp)){
      if(Bin.df.exp[i,j] > 0.1 ){
         Bin.df.exp[i,j] <- 1
      }
   }
}

Bin.df.exp <- -log10(Bin.df.exp)

for(i in 1:nrow(Bin.df.exp)){
   for(j in 1:ncol(Bin.df.exp)){
      if(Bin.df.exp[i,j] > 10 ){
         Bin.df.exp[i,j] <- 10
      }
   }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Greens"))(500)

ht.Bin.exp <- Heatmap(Bin.df.exp, name = "Bins.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

order.ps.exp <- sort(rownames(enrichment.parametric.PS.all.Expression$pval.adj))

PS.df.exp <- t(enrichment.parametric.PS.all.Expression$pval.adj)[,order.ps.exp]

for(i in 1:nrow(PS.df.exp)){
   for(j in 1:ncol(PS.df.exp)){
      if(PS.df.exp[i,j] > 0.1 ){
         PS.df.exp[i,j] <- 1
      }
   }
}

PS.df.exp <- -log10(PS.df.exp)

for(i in 1:nrow(PS.df.exp)){
   for(j in 1:ncol(PS.df.exp)){
      if(PS.df.exp[i,j] > 10 ){
         PS.df.exp[i,j] <- 10
      }
   }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Purples"))(100)

ht.PS.exp <- Heatmap(PS.df.exp, name = "PS.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "white", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

order.psf.exp <- sort(rownames(enrichment.parametric.PSF.all.Expression$pval.adj))

PSF.df.exp <- t(enrichment.parametric.PSF.all.Expression$pval.adj)[,order.psf.exp]

for(i in 1:nrow(PSF.df.exp)){
   for(j in 1:ncol(PSF.df.exp)){
      if(PSF.df.exp[i,j] > 0.1 ){
         PSF.df.exp[i,j] <- 1
      }
   }
}

PSF.df.exp <- -log10(PSF.df.exp)

for(i in 1:nrow(PSF.df.exp)){
   for(j in 1:ncol(PSF.df.exp)){
      if(PSF.df.exp[i,j] > 10 ){
         PSF.df.exp[i,j] <- 10
      }
   }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100)

ht.PSF.exp <- Heatmap(PSF.df.exp, name = "PSF.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

draw(ht.Bin.exp + ht.PS.exp + ht.PSF.exp, ht_gap = unit(1, "cm"))

# splicing

order.bin.exp <- sort(rownames(enrichment.parametric.Bin.all.Expression$pval.adj))

Bin.df.splicing <- data.frame(Bin = order.bin.exp, Factor1 = NA, Factor3 = NA)

for(i in 1:nrow(Bin.df.splicing)){
   hit <- which(rownames(enrichment.parametric.Bin.all.Splicing$pval.adj) == Bin.df.splicing$Bin[i])
   if(length(hit) == 0){
      Bin.df.splicing$Factor1[i] <- 1
      Bin.df.splicing$Factor3[i] <- 1
   }else if(length(hit) > 0){
      Bin.df.splicing$Factor1[i] <- enrichment.parametric.Bin.all.Splicing$pval.adj[hit,1]
      Bin.df.splicing$Factor3[i] <- enrichment.parametric.Bin.all.Splicing$pval.adj[hit,2]
   }
}

rownames(Bin.df.splicing) <- Bin.df.splicing$Bin
Bin.df.splicing <- Bin.df.splicing[,-1]
Bin.df.splicing <- t(Bin.df.splicing)[,order.bin.exp]

for(i in 1:nrow(Bin.df.splicing)){
   for(j in 1:ncol(Bin.df.splicing)){
      if(Bin.df.splicing[i,j] > 0.1 ){
         Bin.df.splicing[i,j] <- 1
      }
   }
}

Bin.df.splicing <- -log10(Bin.df.splicing)

for(i in 1:nrow(Bin.df.splicing)){
   for(j in 1:ncol(Bin.df.splicing)){
      if(Bin.df.splicing[i,j] > 10 ){
         Bin.df.splicing[i,j] <- 10
      }
   }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Greens"))(500)

ht.Bin.splicing <- Heatmap(Bin.df.splicing, name = "Bins.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "white", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

order.ps.exp <- sort(rownames(enrichment.parametric.PS.all.Expression$pval.adj))

PS.df.splicing <- data.frame(PS = order.ps.exp, Factor1 = NA, Factor3 = NA)

for(i in 1:nrow(PS.df.splicing)){
   hit <- which(rownames(enrichment.parametric.PS.all.Splicing$pval.adj) == PS.df.splicing$PS[i])
   if(length(hit) == 0){
      PS.df.splicing$Factor1[i] <- 1
      PS.df.splicing$Factor3[i] <- 1
   }else if(length(hit) > 0){
      PS.df.splicing$Factor1[i] <- enrichment.parametric.PS.all.Splicing$pval.adj[hit,1]
      PS.df.splicing$Factor3[i] <- enrichment.parametric.PS.all.Splicing$pval.adj[hit,2]
   }
}

rownames(PS.df.splicing) <- PS.df.splicing$PS
PS.df.splicing <- PS.df.splicing[,-1]
PS.df.splicing <- t(PS.df.splicing)[,order.ps.exp]

for(i in 1:nrow(PS.df.splicing)){
   for(j in 1:ncol(PS.df.splicing)){
      if(PS.df.splicing[i,j] > 0.1 ){
         PS.df.splicing[i,j] <- 1
      }
   }
}

PS.df.splicing <- -log10(PS.df.splicing)

for(i in 1:nrow(PS.df.splicing)){
   for(j in 1:ncol(PS.df.splicing)){
      if(PS.df.splicing[i,j] > 10 ){
         PS.df.splicing[i,j] <- 10
      }
   }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Purples"))(100)

ht.PS.splicing <- Heatmap(PS.df.splicing, name = "PS.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

order.psf.exp <- sort(rownames(enrichment.parametric.PSF.all.Expression$pval.adj))

PSF.df.splicing <- data.frame(PSF = order.psf.exp, Factor1 = NA, Factor3 = NA)

for(i in 1:nrow(PSF.df.splicing)){
   hit <- which(rownames(enrichment.parametric.PSF.all.Splicing$pval.adj) == PSF.df.splicing$PSF[i])
   if(length(hit) == 0){
      PSF.df.splicing$Factor1[i] <- 1
      PSF.df.splicing$Factor3[i] <- 1
   }else if(length(hit) > 0){
      PSF.df.splicing$Factor1[i] <- enrichment.parametric.PSF.all.Splicing$pval.adj[hit,1]
      PSF.df.splicing$Factor3[i] <- enrichment.parametric.PSF.all.Splicing$pval.adj[hit,2]
   }
}

rownames(PSF.df.splicing) <- PSF.df.splicing$PSF
PSF.df.splicing <- PSF.df.splicing[,-1]
PSF.df.splicing <- t(PSF.df.splicing)[,order.psf.exp]

for(i in 1:nrow(PSF.df.splicing)){
   for(j in 1:ncol(PSF.df.splicing)){
      if(PSF.df.splicing[i,j] > 0.1 ){
         PSF.df.splicing[i,j] <- 1
      }
   }
}

PSF.df.splicing <- -log10(PSF.df.splicing)

for(i in 1:nrow(PSF.df.splicing)){
   for(j in 1:ncol(PSF.df.splicing)){
      if(PSF.df.splicing[i,j] > 10 ){
         PSF.df.splicing[i,j] <- 10
      }
   }
}
col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100)

ht.PSF.splicing <- Heatmap(PSF.df.splicing, name = "PSF.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "white", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

draw(ht.Bin.splicing + ht.PS.exp + ht.PSF.splicing, ht_gap = unit(1, "cm"))

#### 4.1 Biotic/all stresses ####

#### Evolutionary Transcriptomics ####

StressSamples <- c(13:16, 18,19,21,22,24:35,49:137)
Stress.TAI <- as.data.frame(t(data0)[,StressSamples])
Stress.TAI$GeneID <- rownames(Stress.TAI) 

Stress.TAI.df <- data.frame(Phylostratum = NA,
                             GeneID = Stress.TAI$GeneID,
                             Control_Heat = rowMeans(Stress.TAI[,grep("HS_Control", colnames(Stress.TAI))]),
                             Stress_Heat_1 = rowMeans(Stress.TAI[,grep("HS_Stress_1", colnames(Stress.TAI))]),
                             Stress_Heat_2 = rowMeans(Stress.TAI[,grep("HS_Stress_2", colnames(Stress.TAI))]),
                             Control_DO = rowMeans(Stress.TAI[,grep("DO_Control", colnames(Stress.TAI))]),
                             Stress_DO_1 = rowMeans(Stress.TAI[,grep("DO_Stress_DO_1", colnames(Stress.TAI))]),
                             Stress_DO_2 = rowMeans(Stress.TAI[,grep("DO_Stress_DO_2", colnames(Stress.TAI))]),
                             Control_PH = rowMeans(Stress.TAI[,grep("PH_Control", colnames(Stress.TAI))]),
                             Stress_PH = rowMeans(Stress.TAI[,grep("PH_Stress_PH", colnames(Stress.TAI))]),
                             Control_FU_1 = rowMeans(Stress.TAI[,grep("FU_Control_Controldmg_1", colnames(Stress.TAI))]),
                             Control_FU_2 = rowMeans(Stress.TAI[,grep("FU_Control_Controldmg_2", colnames(Stress.TAI))]),
                             Control_FU_3 = rowMeans(Stress.TAI[,grep("FU_Control_Controldmg_3", colnames(Stress.TAI))]),
                             Control_FU_4 = rowMeans(Stress.TAI[,grep("FU_Control_Controldmg_4", colnames(Stress.TAI))]),
                             Stress_FU_1 = rowMeans(Stress.TAI[,grep("FU_Stress_FU_1", colnames(Stress.TAI))]),
                             Stress_FU_2 = rowMeans(Stress.TAI[,grep("FU_Stress_FU_2", colnames(Stress.TAI))]),
                             Stress_FU_3 = rowMeans(Stress.TAI[,grep("FU_Stress_FU_3", colnames(Stress.TAI))]),
                             Stress_FU_4 = rowMeans(Stress.TAI[,grep("FU_Stress_FU_4", colnames(Stress.TAI))])
                             )


View(Stress.TAI.df)
for(i in 1:nrow(Stress.TAI.df)){
   Stress.TAI.df$Phylostratum[i] <- Gene.annotation$PhyloStratum[which(Gene.annotation$GeneID == Stress.TAI.df$GeneID[i])]
}

Stress.TAI.df <- Stress.TAI.df[!is.na(Stress.TAI.df$Phylostratum),]
dim(Stress.TAI.df)
PlotSignature(Stress.TAI.df)
PlotSignature(Stress.TAI.df[,c(1:2,3:5)]) # HS: 0.00164
PlotSignature(Stress.TAI.df[,c(1:2,6:8)]) # DO: 0.000187
PlotSignature(Stress.TAI.df[,c(1:2,9:10)]) # PH: 0.3
PlotSignature(Stress.TAI.df[,c(1:2,11:18)]) # PH: 0.232

Stress.TAI.jic <- Stress.TAI

#### MOFA ####

### 4.1.1 Split Splicing Events by AS type and general AS

Splicing.IR <- Splicing[which(rownames(Splicing) %in% Event.annotation.split$EVENTID[Event.annotation.split$Type == "IR"]),]
dim(Splicing.IR)
Splicing.ES <- Splicing[which(rownames(Splicing) %in% Event.annotation.split$EVENTID[Event.annotation.split$Type == "ES"]),]
dim(Splicing.ES)
Splicing.AltAD <- Splicing[which(rownames(Splicing) %in% Event.annotation.split$EVENTID[Event.annotation.split$Type == "AltAD"]),]
dim(Splicing.AltAD)
Splicing.AS <- Splicing[which(rownames(Splicing) %in% Event.annotation.split$EVENTID[Event.annotation.split$Type == "AS"]),]
dim(Splicing.AS)

AllTranscriptionalList.split <- list(Expression = as.matrix(Expression), IR = as.matrix(Splicing.IR),
                                     ES = as.matrix(Splicing.ES), AltAD = as.matrix(Splicing.AltAD),
                                     AS = as.matrix(Splicing.AS))

AllTranscriptionalList <- list(Expression = as.matrix(Expression), Splicing = as.matrix(Splicing))

StressSamplesBiotic <- c(11:16, 18:19, 
                   #21:22, 
                   24:35, 
                   #53:56, 
                   57:66, 67:72, 73:93, 94:141) 

StressSamplesAll <- c(11:16, 18:19, 
                         21:22, 
                         24:35, 
                         53:56, 
                         57:66, 67:72, 73:93, 94:141) 

AllTranscriptionalList.split.stress.biotic <- lapply(AllTranscriptionalList.split, function(x){
   x[,StressSamplesBiotic]
})

AllTranscriptionalList.split.stress.all <- lapply(AllTranscriptionalList.split, function(x){
   x[,StressSamplesAll]
})

AllTranscriptionalList.stress.biotic <- lapply(AllTranscriptionalList, function(x){
   x[,StressSamplesBiotic]
})

AllTranscriptionalList.stress.all <- lapply(AllTranscriptionalList, function(x){
   x[,StressSamplesAll]
})

### define groups

groupsBiotic <- c(rep("FU", 6),  rep("DO", 2), rep("PH", 12), rep("DO", 10), rep("PH", 6), rep("FU", 69))
groupsAll <- c(rep("FU", 6),  rep("DO", 2), rep("HS",2), rep("PH", 12), rep("HS", 4), rep("DO", 10), rep("PH", 6), rep("FU", 69))

### filter by variance features

## Splitted biotic

AllTranscriptionalList.split.stress.biotic.Expression.vars <- apply(AllTranscriptionalList.split.stress.biotic$Expression, 1, var)
AllTranscriptionalList.split.stress.biotic.IR.vars <- apply(AllTranscriptionalList.split.stress.biotic$IR, 1, function(x){
   var(x, na.rm = TRUE)
})
AllTranscriptionalList.split.stress.biotic.ES.vars <- apply(AllTranscriptionalList.split.stress.biotic$ES, 1, function(x){
   var(x, na.rm = TRUE)
})
AllTranscriptionalList.split.stress.biotic.AltAD.vars <- apply(AllTranscriptionalList.split.stress.biotic$AltAD, 1, function(x){
   var(x, na.rm = TRUE)
})
AllTranscriptionalList.split.stress.biotic.AS.vars <- apply(AllTranscriptionalList.split.stress.biotic$AS, 1, function(x){
   var(x, na.rm = TRUE)
})

AllTranscriptionalList.split.stress.biotic.all <-  AllTranscriptionalList.split.stress.biotic

AllTranscriptionalList.split.stress.biotic.all$Expression <- AllTranscriptionalList.split.stress.biotic.all$Expression[rownames(AllTranscriptionalList.split.stress.biotic.all$Expression) %in% names(AllTranscriptionalList.split.stress.biotic.Expression.vars)[AllTranscriptionalList.split.stress.biotic.Expression.vars > 0],]
AllTranscriptionalList.split.stress.biotic.all$IR <- AllTranscriptionalList.split.stress.biotic.all$IR[rownames(AllTranscriptionalList.split.stress.biotic.all$IR) %in% names(AllTranscriptionalList.split.stress.biotic.IR.vars)[AllTranscriptionalList.split.stress.biotic.IR.vars > 0],]
AllTranscriptionalList.split.stress.biotic.all$ES <- AllTranscriptionalList.split.stress.biotic.all$ES[rownames(AllTranscriptionalList.split.stress.biotic.all$ES) %in% names(AllTranscriptionalList.split.stress.biotic.ES.vars)[AllTranscriptionalList.split.stress.biotic.ES.vars > 0],]
AllTranscriptionalList.split.stress.biotic.all$AltAD <- AllTranscriptionalList.split.stress.biotic.all$AltAD[rownames(AllTranscriptionalList.split.stress.biotic.all$AltAD) %in% names(AllTranscriptionalList.split.stress.biotic.AltAD.vars)[AllTranscriptionalList.split.stress.biotic.AltAD.vars > 0],]
AllTranscriptionalList.split.stress.biotic.all$AS <- AllTranscriptionalList.split.stress.biotic.all$AS[rownames(AllTranscriptionalList.split.stress.biotic.all$AS) %in% names(AllTranscriptionalList.split.stress.biotic.AS.vars)[AllTranscriptionalList.split.stress.biotic.AS.vars > 0],]

## Splitted Stressall

AllTranscriptionalList.split.stress.all.Expression.vars <- apply(AllTranscriptionalList.split.stress.all$Expression, 1, var)
AllTranscriptionalList.split.stress.all.IR.vars <- apply(AllTranscriptionalList.split.stress.all$IR, 1, function(x){
   var(x, na.rm = TRUE)
})
AllTranscriptionalList.split.stress.all.ES.vars <- apply(AllTranscriptionalList.split.stress.all$ES, 1, function(x){
   var(x, na.rm = TRUE)
})
AllTranscriptionalList.split.stress.all.AltAD.vars <- apply(AllTranscriptionalList.split.stress.all$AltAD, 1, function(x){
   var(x, na.rm = TRUE)
})
AllTranscriptionalList.split.stress.all.AS.vars <- apply(AllTranscriptionalList.split.stress.all$AS, 1, function(x){
   var(x, na.rm = TRUE)
})

AllTranscriptionalList.split.stress.all.all <-  AllTranscriptionalList.split.stress.all

AllTranscriptionalList.split.stress.all.all$Expression <- AllTranscriptionalList.split.stress.all.all$Expression[rownames(AllTranscriptionalList.split.stress.all.all$Expression) %in% names(AllTranscriptionalList.split.stress.all.Expression.vars)[AllTranscriptionalList.split.stress.all.Expression.vars > 0],]
AllTranscriptionalList.split.stress.all.all$IR <- AllTranscriptionalList.split.stress.all.all$IR[rownames(AllTranscriptionalList.split.stress.all.all$IR) %in% names(AllTranscriptionalList.split.stress.all.IR.vars)[AllTranscriptionalList.split.stress.all.IR.vars > 0],]
AllTranscriptionalList.split.stress.all.all$ES <- AllTranscriptionalList.split.stress.all.all$ES[rownames(AllTranscriptionalList.split.stress.all.all$ES) %in% names(AllTranscriptionalList.split.stress.all.ES.vars)[AllTranscriptionalList.split.stress.all.ES.vars > 0],]
AllTranscriptionalList.split.stress.all.all$AltAD <- AllTranscriptionalList.split.stress.all.all$AltAD[rownames(AllTranscriptionalList.split.stress.all.all$AltAD) %in% names(AllTranscriptionalList.split.stress.all.AltAD.vars)[AllTranscriptionalList.split.stress.all.AltAD.vars > 0],]
AllTranscriptionalList.split.stress.all.all$AS <- AllTranscriptionalList.split.stress.all.all$AS[rownames(AllTranscriptionalList.split.stress.all.all$AS) %in% names(AllTranscriptionalList.split.stress.all.AS.vars)[AllTranscriptionalList.split.stress.all.AS.vars > 0],]

lapply(AllTranscriptionalList.split.stress.all, dim)
lapply(AllTranscriptionalList.split.stress.all.all, dim)

## Export StressBiotic splitted all

write.table(AllTranscriptionalList.split.stress.biotic.all$Expression, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressBiotic_split_all/Expression.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.stress.biotic.all$IR, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressBiotic_split_all/IR.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.stress.biotic.all$ES, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressBiotic_split_all/ES.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.stress.biotic.all$AltAD, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressBiotic_split_all/AltAD.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.stress.biotic.all$AS, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressBiotic_split_all/AS.txt", sep = "\t", col.names = T, row.names = T, quote = F)

## Export StressAll splitted all

write.table(AllTranscriptionalList.split.stress.all.all$Expression, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressAll_split_all/Expression.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.stress.all.all$IR, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressAll_split_all/IR.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.stress.all.all$ES, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressAll_split_all/ES.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.stress.all.all$AltAD, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressAll_split_all/AltAD.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.stress.all.all$AS, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressAll_split_all/AS.txt", sep = "\t", col.names = T, row.names = T, quote = F)

## All StressBiotic all

## Export StressAll all all

## All StressAll all

## Export StressAll splitted all

## Split StressBiotic Top 5,000 AS 10,000 Genes: per view and per group

AllTranscriptionalList.split.stress.biotic.top <- AllTranscriptionalList.split.stress.biotic

FU.Genes <- apply(AllTranscriptionalList.split.stress.biotic$Expression[,grep(pattern = "FU_", x = colnames(AllTranscriptionalList.split.stress.biotic$Expression), fixed = T)], 1, var)
FU.Genes <- names(FU.Genes[order(FU.Genes, decreasing = T)])[1:10000]
PH.Genes <- apply(AllTranscriptionalList.split.stress.biotic$Expression[,grep(pattern = "PH_", x = colnames(AllTranscriptionalList.split.stress.biotic$Expression), fixed = T)], 1, var)
PH.Genes <- names(PH.Genes[order(PH.Genes, decreasing = T)])[1:10000]
DO.Genes <- apply(AllTranscriptionalList.split.stress.biotic$Expression[,grep(pattern = "DO_", x = colnames(AllTranscriptionalList.split.stress.biotic$Expression), fixed = T)], 1, var)
DO.Genes <- names(DO.Genes[order(DO.Genes, decreasing = T)])[1:10000]
Expression.Genes <- unique(c(FU.Genes, PH.Genes, DO.Genes)) 
AllTranscriptionalList.split.stress.biotic.top$Expression <- AllTranscriptionalList.split.stress.biotic.top$Expression[rownames(AllTranscriptionalList.split.stress.biotic.top$Expression) %in% Expression.Genes,]


FU.IR <- apply(AllTranscriptionalList.split.stress.biotic$IR[,grep(pattern = "FU_", x = colnames(AllTranscriptionalList.split.stress.biotic$IR), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
FU.IR <- names(FU.IR[order(FU.IR, decreasing = T)])[1:5000]
PH.IR <- apply(AllTranscriptionalList.split.stress.biotic$IR[,grep(pattern = "PH_", x = colnames(AllTranscriptionalList.split.stress.biotic$IR), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
PH.IR <- names(PH.IR[order(PH.IR, decreasing = T)])[1:5000]
DO.IR <- apply(AllTranscriptionalList.split.stress.biotic$IR[,grep(pattern = "DO_", x = colnames(AllTranscriptionalList.split.stress.biotic$IR), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
DO.IR <- names(DO.IR[order(DO.IR, decreasing = T)])[1:5000]
IR.Events <- unique(c(FU.IR, PH.IR, DO.IR)) 
AllTranscriptionalList.split.stress.biotic.top$IR <- AllTranscriptionalList.split.stress.biotic.top$IR[rownames(AllTranscriptionalList.split.stress.biotic.top$IR) %in% IR.Events,]


FU.ES <- apply(AllTranscriptionalList.split.stress.biotic$ES[,grep(pattern = "FU_", x = colnames(AllTranscriptionalList.split.stress.biotic$ES), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
FU.ES <- names(FU.ES[order(FU.ES, decreasing = T)])[1:5000]
PH.ES <- apply(AllTranscriptionalList.split.stress.biotic$ES[,grep(pattern = "PH_", x = colnames(AllTranscriptionalList.split.stress.biotic$ES), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
PH.ES <- names(PH.ES[order(PH.ES, decreasing = T)])[1:5000]
DO.ES <- apply(AllTranscriptionalList.split.stress.biotic$ES[,grep(pattern = "DO_", x = colnames(AllTranscriptionalList.split.stress.biotic$ES), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
DO.ES <- names(DO.ES[order(DO.ES, decreasing = T)])[1:5000]
ES.Events <- unique(c(FU.ES, PH.ES, DO.ES)) 
AllTranscriptionalList.split.stress.biotic.top$ES <- AllTranscriptionalList.split.stress.biotic.top$ES[rownames(AllTranscriptionalList.split.stress.biotic.top$ES) %in% ES.Events,]


FU.AltAD <- apply(AllTranscriptionalList.split.stress.biotic$AltAD[,grep(pattern = "FU_", x = colnames(AllTranscriptionalList.split.stress.biotic$AltAD), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
FU.AltAD <- names(FU.AltAD[order(FU.AltAD, decreasing = T)])[1:5000]
PH.AltAD <- apply(AllTranscriptionalList.split.stress.biotic$AltAD[,grep(pattern = "PH_", x = colnames(AllTranscriptionalList.split.stress.biotic$AltAD), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
PH.AltAD <- names(PH.AltAD[order(PH.AltAD, decreasing = T)])[1:5000]
DO.AltAD <- apply(AllTranscriptionalList.split.stress.biotic$AltAD[,grep(pattern = "DO_", x = colnames(AllTranscriptionalList.split.stress.biotic$AltAD), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
DO.AltAD <- names(DO.AltAD[order(DO.AltAD, decreasing = T)])[1:5000]
AltAD.Events <- unique(c(FU.AltAD, PH.AltAD, DO.AltAD)) 
AllTranscriptionalList.split.stress.biotic.top$AltAD <- AllTranscriptionalList.split.stress.biotic.top$AltAD[rownames(AllTranscriptionalList.split.stress.biotic.top$AltAD) %in% AltAD.Events,]


FU.AS <- apply(AllTranscriptionalList.split.stress.biotic$AS[,grep(pattern = "FU_", x = colnames(AllTranscriptionalList.split.stress.biotic$AS), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
FU.AS <- names(FU.AS[order(FU.AS, decreasing = T)])[1:5000]
PH.AS <- apply(AllTranscriptionalList.split.stress.biotic$AS[,grep(pattern = "PH_", x = colnames(AllTranscriptionalList.split.stress.biotic$AS), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
PH.AS <- names(PH.AS[order(PH.AS, decreasing = T)])[1:5000]
DO.AS <- apply(AllTranscriptionalList.split.stress.biotic$AS[,grep(pattern = "DO_", x = colnames(AllTranscriptionalList.split.stress.biotic$AS), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
DO.AS <- names(DO.AS[order(DO.AS, decreasing = T)])[1:5000]
AS.Events <- unique(c(FU.AS, PH.AS, DO.AS)) 
AllTranscriptionalList.split.stress.biotic.top$AS <- AllTranscriptionalList.split.stress.biotic.top$AS[rownames(AllTranscriptionalList.split.stress.biotic.top$AS) %in% AS.Events,]

lapply(AllTranscriptionalList.split.stress.biotic, dim)
lapply(AllTranscriptionalList.split.stress.biotic.top, dim)

write.table(AllTranscriptionalList.split.stress.biotic.top$Expression, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressBiotic_split_top/Expression.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.stress.biotic.top$IR, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressBiotic_split_top/IR.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.stress.biotic.top$ES, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressBiotic_split_top/ES.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.stress.biotic.top$AltAD, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressBiotic_split_top/AltAD.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.stress.biotic.top$AS, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressBiotic_split_top/AS.txt", sep = "\t", col.names = T, row.names = T, quote = F)

## Split StressBiotic Top 5,000 AS 10,000 Genes: per view and per group

AllTranscriptionalList.split.stress.all.top <- AllTranscriptionalList.split.stress.all

FU.Genes <- apply(AllTranscriptionalList.split.stress.all$Expression[,grep(pattern = "FU_", x = colnames(AllTranscriptionalList.split.stress.all$Expression), fixed = T)], 1, var)
FU.Genes <- names(FU.Genes[order(FU.Genes, decreasing = T)])[1:10000]
PH.Genes <- apply(AllTranscriptionalList.split.stress.all$Expression[,grep(pattern = "PH_", x = colnames(AllTranscriptionalList.split.stress.all$Expression), fixed = T)], 1, var)
PH.Genes <- names(PH.Genes[order(PH.Genes, decreasing = T)])[1:10000]
DO.Genes <- apply(AllTranscriptionalList.split.stress.all$Expression[,grep(pattern = "DO_", x = colnames(AllTranscriptionalList.split.stress.all$Expression), fixed = T)], 1, var)
DO.Genes <- names(DO.Genes[order(DO.Genes, decreasing = T)])[1:10000]
Expression.Genes <- unique(c(FU.Genes, PH.Genes, DO.Genes)) 
AllTranscriptionalList.split.stress.all.top$Expression <- AllTranscriptionalList.split.stress.all.top$Expression[rownames(AllTranscriptionalList.split.stress.all.top$Expression) %in% Expression.Genes,]


FU.IR <- apply(AllTranscriptionalList.split.stress.all$IR[,grep(pattern = "FU_", x = colnames(AllTranscriptionalList.split.stress.all$IR), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
FU.IR <- names(FU.IR[order(FU.IR, decreasing = T)])[1:5000]
PH.IR <- apply(AllTranscriptionalList.split.stress.all$IR[,grep(pattern = "PH_", x = colnames(AllTranscriptionalList.split.stress.all$IR), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
PH.IR <- names(PH.IR[order(PH.IR, decreasing = T)])[1:5000]
DO.IR <- apply(AllTranscriptionalList.split.stress.all$IR[,grep(pattern = "DO_", x = colnames(AllTranscriptionalList.split.stress.all$IR), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
DO.IR <- names(DO.IR[order(DO.IR, decreasing = T)])[1:5000]
IR.Events <- unique(c(FU.IR, PH.IR, DO.IR)) 
AllTranscriptionalList.split.stress.all.top$IR <- AllTranscriptionalList.split.stress.all.top$IR[rownames(AllTranscriptionalList.split.stress.all.top$IR) %in% IR.Events,]


FU.ES <- apply(AllTranscriptionalList.split.stress.all$ES[,grep(pattern = "FU_", x = colnames(AllTranscriptionalList.split.stress.all$ES), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
FU.ES <- names(FU.ES[order(FU.ES, decreasing = T)])[1:5000]
PH.ES <- apply(AllTranscriptionalList.split.stress.all$ES[,grep(pattern = "PH_", x = colnames(AllTranscriptionalList.split.stress.all$ES), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
PH.ES <- names(PH.ES[order(PH.ES, decreasing = T)])[1:5000]
DO.ES <- apply(AllTranscriptionalList.split.stress.all$ES[,grep(pattern = "DO_", x = colnames(AllTranscriptionalList.split.stress.all$ES), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
DO.ES <- names(DO.ES[order(DO.ES, decreasing = T)])[1:5000]
ES.Events <- unique(c(FU.ES, PH.ES, DO.ES)) 
AllTranscriptionalList.split.stress.all.top$ES <- AllTranscriptionalList.split.stress.all.top$ES[rownames(AllTranscriptionalList.split.stress.all.top$ES) %in% ES.Events,]


FU.AltAD <- apply(AllTranscriptionalList.split.stress.all$AltAD[,grep(pattern = "FU_", x = colnames(AllTranscriptionalList.split.stress.all$AltAD), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
FU.AltAD <- names(FU.AltAD[order(FU.AltAD, decreasing = T)])[1:5000]
PH.AltAD <- apply(AllTranscriptionalList.split.stress.all$AltAD[,grep(pattern = "PH_", x = colnames(AllTranscriptionalList.split.stress.all$AltAD), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
PH.AltAD <- names(PH.AltAD[order(PH.AltAD, decreasing = T)])[1:5000]
DO.AltAD <- apply(AllTranscriptionalList.split.stress.all$AltAD[,grep(pattern = "DO_", x = colnames(AllTranscriptionalList.split.stress.all$AltAD), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
DO.AltAD <- names(DO.AltAD[order(DO.AltAD, decreasing = T)])[1:5000]
AltAD.Events <- unique(c(FU.AltAD, PH.AltAD, DO.AltAD)) 
AllTranscriptionalList.split.stress.all.top$AltAD <- AllTranscriptionalList.split.stress.all.top$AltAD[rownames(AllTranscriptionalList.split.stress.all.top$AltAD) %in% AltAD.Events,]


FU.AS <- apply(AllTranscriptionalList.split.stress.all$AS[,grep(pattern = "FU_", x = colnames(AllTranscriptionalList.split.stress.all$AS), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
FU.AS <- names(FU.AS[order(FU.AS, decreasing = T)])[1:5000]
PH.AS <- apply(AllTranscriptionalList.split.stress.all$AS[,grep(pattern = "PH_", x = colnames(AllTranscriptionalList.split.stress.all$AS), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
PH.AS <- names(PH.AS[order(PH.AS, decreasing = T)])[1:5000]
DO.AS <- apply(AllTranscriptionalList.split.stress.all$AS[,grep(pattern = "DO_", x = colnames(AllTranscriptionalList.split.stress.all$AS), fixed = T)], 1, function(x){
   var(x, na.rm = TRUE)
})
DO.AS <- names(DO.AS[order(DO.AS, decreasing = T)])[1:5000]
AS.Events <- unique(c(FU.AS, PH.AS, DO.AS)) 
AllTranscriptionalList.split.stress.all.top$AS <- AllTranscriptionalList.split.stress.all.top$AS[rownames(AllTranscriptionalList.split.stress.all.top$AS) %in% AS.Events,]

lapply(AllTranscriptionalList.split.stress.all, dim)
lapply(AllTranscriptionalList.split.stress.all.top, dim)

write.table(AllTranscriptionalList.split.stress.all.top$Expression, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressAll_split_top/Expression.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.stress.all.top$IR, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressAll_split_top/IR.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.stress.all.top$ES, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressAll_split_top/ES.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.stress.all.top$AltAD, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressAll_split_top/AltAD.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(AllTranscriptionalList.split.stress.all.top$AS, file = "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/Stress/MOFA_StressAll_split_top/AS.txt", sep = "\t", col.names = T, row.names = T, quote = F)

## Prepare list for MOFA model tunning

AllTranscriptionalList.split.stress.all.all <- list(Expression = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressAll_split_all/Expression.txt")),
                                             IR = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressAll_split_all/IR.txt")),
                                             ES = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressAll_split_all/ES.txt")),
                                             AltAD = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressAll_split_all/AltAD.txt")),
                                             AS = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressAll_split_all/AS.txt"))
)

AllTranscriptionalList.split.stress.all.top <- list(Expression = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressAll_split_top/Expression.txt")),
                                                 IR = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressAll_split_top/IR.txt")),
                                                 ES = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressAll_split_top/ES.txt")),
                                                 AltAD = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressAll_split_top/AltAD.txt")),
                                                 AS = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressAll_split_top/AS.txt"))
)

AllTranscriptionalList.split.stress.biotic.all <- list(Expression = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressBiotic_split_all/Expression.txt")),
                                                       IR = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressBiotic_split_all/IR.txt")),
                                                       ES = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressBiotic_split_all/ES.txt")),
                                                       AltAD = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressBiotic_split_all/AltAD.txt")),
                                                       AS = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressBiotic_split_all/AS.txt"))
)

AllTranscriptionalList.split.stress.biotic.top <- list(Expression = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressBiotic_split_top/Expression.txt")),
                                                       IR = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressBiotic_split_top/IR.txt")),
                                                       ES = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressBiotic_split_top/ES.txt")),
                                                       AltAD = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressBiotic_split_top/AltAD.txt")),
                                                       AS = as.matrix(read.table("/data/Transcriptional_MOFA/Stress/MOFA_StressBiotic_split_top/AS.txt"))
)

Check <- list(Expression = as.matrix(read.table("/data/Transcriptional_MOFA/ExpressionCheckStress.txt")))
Check <- list(Expression = as.matrix(read.table("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/ExpressionCheckStress.txt")))
Check <- Check$Expression[,StressSamplesAll]
write.table(Check, "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/ExpressionCheckStress.txt", sep = "\t", col.names = T, row.names = T, quote = F)   

groupsBiotic <- c(rep("FU", 6),  rep("DO", 2), rep("PH", 12), rep("DO", 10), rep("PH", 6), rep("FU", 69))
groupsAll <- c(rep("FU", 6),  rep("DO", 2), rep("HS",2), rep("PH", 12), rep("HS", 4), rep("DO", 10), rep("PH", 6), rep("FU", 69))

MOFAgrouped <- create_mofa(data = AllTranscriptionalList.split.tissues, groups = groupsAll, save_metadata = TRUE)

## 4.1.2 Prepare MOFA and run

plot_data_overview(MOFAgrouped)
sample_metadata <- data.frame(
   sample = colnames(AllTranscriptionalList.split.stress.all.all$Expression),
   stress = labelitass$Stress[StressSamplesAll],
   time = labelitass$Intensity[StressSamplesAll],
   genotype = labelitass$Genotype[StressSamplesAll],
   treatment = labelitass$Treatment[StressSamplesAll],
   treatment.time = paste0(labelitass$Treatment[StressSamplesAll], "_", labelitass$Intensity[StressSamplesAll]),
   group = groupsAll,
   sample_number = c(1:length(groupsAll))
)

sample_metadata <- data.frame(
   sample = colnames(AllTranscriptionalList.split.stress.biotic.all$Expression),
   stress = labelitass$Stress[StressSamplesBiotic],
   time = labelitass$Intensity[StressSamplesBiotic],
   genotype = labelitass$Genotype[StressSamplesBiotic],
   treatment = labelitass$Treatment[StressSamplesBiotic],
   group = groupsBiotic,
   sample_number = c(1:length(groupsBiotic))
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
outfileGrouped <- "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/MOFA/"

MOFAgrouped <- prepare_mofa(MOFAgrouped, data_options = data_opts, model_options = model_opts, training_options = train_opts)
MOFAgrouped.trained <- run_mofa(MOFAgrouped, outfile = outfileGrouped)

MOFAgrouped <- prepare_mofa(MOFAgrouped, data_options = data_opts, model_options = model_opts, training_options = train_opts)
MOFAgrouped.trained <- run_mofa(MOFAgrouped, outfile = outfileGrouped)

MOFAgrouped <- load_model("F:/ATLAS_Pra/Protein_module/0.Process/Transcriptional_MOFA/Stress/MOFA_StressAll_split_top/MOFA_stressAll_split_top_400_false.hdf5")
MOFAgrouped <- load_model("F:/ATLAS_Pra/Protein_module/0.Process/Transcriptional_MOFA/Check_400_false.hdf5")
MOFAgrouped.trained <- MOFAgrouped

## 4.1.3 Plot General

plot_factor_cor(MOFAgrouped.trained)
plot_variance_explained(MOFAgrouped.trained, x = "group", y = "factor")
plot_variance_explained(MOFAgrouped.trained, x = "group", y = "factor", plot_total = T)

MOFAgrouped.trained@cache$variance_explained
MOFAgrouped.trained@cache$variance_explained$r2_per_factor$DO
MOFAgrouped.trained@cache$variance_explained$r2_per_factor$FU
MOFAgrouped.trained@cache$variance_explained$r2_per_factor$HS
MOFAgrouped.trained@cache$variance_explained$r2_per_factor$PH

a <- plot_variance_explained(MOFAgrouped.trained, x="view", y="factor", factors = c(1:8))
b <- plot_variance_explained(MOFAgrouped.trained, x="group", y="factor", plot_total = T, factors = c(1:8))[[2]]
PlotVariance <- grid.arrange(b, a, ncol = 1, nrow = 2)

## 4.1.4 LFs Inspection

# Single

MOFAgrouped.trained@samples_metadata <- sample_metadata

plot_factor(MOFAgrouped.trained, 
            factor = c(3),
            color_by = "treatment.time",
            shape_by = "genotype", dot_size = 4
)

# Multiple

plot_factors(MOFAgrouped.trained, 
             factors = c(1,2,3),
             color_by = "time", dot_size = 4, shape_by = "genotype"
)

plot_factors(MOFAgrouped.trained, 
             factors = c(1,2,3),
             color_by = "stress", dot_size = 3, shape_by = "sex"
)

# Weights of biologically relevant LFs

LF1 <- plot_top_weights(MOFAgrouped.trained,
                        view = c("Expression", "IR", "ES", "AltAD", "AS"),
                        factor = 1,
                        nfeatures = 10,     # Number of features to highlight
                        scale = T,          # Scale weights from -1 to 1
                        abs = F  ) 

LF2 <- plot_top_weights(MOFAgrouped.trained,
                        view = c("Expression", "IR", "ES", "AltAD", "AS"),
                        factor = 2,
                        nfeatures = 10,     # Number of features to highlight
                        scale = T,          # Scale weights from -1 to 1
                        abs = F  ) 

LF3 <- plot_top_weights(MOFAgrouped.trained,
                        view = c("Expression", "IR", "ES", "AltAD", "AS"),
                        factor = 3,
                        nfeatures = 10,     # Number of features to highlight
                        scale = T,          # Scale weights from -1 to 1
                        abs = F  )

gridExtra::grid.arrange(LF1, LF2, LF3, ncol = 3, nrow = 1)

# Enrichment Analyses #

# Mercator

Bin_db

Bin_db_splicing <- list()
bin.names.splicing <- levels(as.factor(Event.annotation.bin$Mercator4.Bin))
for(i in bin.names.splicing){
   Bin_db_splicing[[i]] <- Event.annotation.bin$EVENTID[Event.annotation.bin$Mercator4.Bin == i]
}
names(Bin_db_splicing) <- c("PS", "Redox hom", "Phytohormone act", "Chromatin org", "Cell division", "DNA dmg response", "RNA biosynthesis", "RNA processing", 
                            "Protein biosynthesis", "Protein modification", "Protein hom",
                            "Cellular respiration", "Cytoskeleton org", "Cell wall org", "Vessicle trafficking", "Protein translocation", "Solute transport",
                            "Nutrient uptake", "External stimuli resp", "Multi-process reg", "Plant reproduction", "Plant organogen",
                            "Carbohydrate met", "Clade-specific met", "Unknown", 
                            "Aminoacid met", "Lipid met", "Enzyme classification", "Nucleotide met", "Coenzyme met", "Polyamine met", "Secondary met")


# Gene Age

Gene.annotation <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Expression/Gene_annotation.txt", header = T)
Gene.annotation.age <- Gene.annotation[!is.na(Gene.annotation$PhyloStratum),]

table(Gene.annotation.age$PhyloStratum)

PS_db <- list()
PS_names <- levels(as.factor(Gene.annotation.age$PhyloStratum))

for(i in PS_names){
   PS_db[[i]] <- Gene.annotation.age$GeneID[Gene.annotation.age$PhyloStratum == i]
}

names(PS_db) <- paste0("PS", names(PS_db))
PS_db

PS_db_splicing <- list()
PS_names_splicing <- levels(as.factor(Event.annotation.bin$Phylostratum))
for(i in PS_names_splicing){
   PS_db_splicing[[i]] <- unique(Event.annotation.bin$EVENTID[Event.annotation.bin$Phylostratum == i])
}

names(PS_db_splicing) <- paste0("PS", names(PS_db_splicing))
PS_db_splicing

# Gene family founder

table(Gene.annotation$GeneFamilyFounder.Age)

Gene.annotation.family <- Gene.annotation
Gene.annotation.family <- Gene.annotation.family[!is.na(Gene.annotation.family$GeneFamilyFounder.Age),]

PSF_db <- list()
PSF_names <- levels(as.factor(Gene.annotation.family$GeneFamilyFounder.Age))

for(i in PSF_names){
   PSF_db[[i]] <- Gene.annotation.family$GeneID[Gene.annotation.family$GeneFamilyFounder.Age == i]
}

names(PSF_db) <- paste0("PS", names(PSF_db))
PSF_db

PSF_db_splicing <- list()
PSF_names_splicing <- levels(as.factor(Event.annotation.bin$GeneFamilyFounder.Age))

for(i in PSF_names_splicing){
   PSF_db_splicing[[i]] <- unique(Event.annotation.bin$EVENTID[Event.annotation.bin$GeneFamilyFounder.Age == i])
}

names(PSF_db_splicing) <- paste0("PS", names(PSF_db_splicing))
PSF_db_splicing

# Compute enrichments; LF1-4; Pos and neg; Mercator, Gene Age and Family foundation

Bin.paths <- list_to_matrix(Bin_db)
PS.paths <- list_to_matrix(PS_db)
PSF.paths <- list_to_matrix(PSF_db)

Bin_splicing.paths <- list_to_matrix(Bin_db_splicing)
PS_splicing.paths <- list_to_matrix(PS_db_splicing)
PSF_splicing.paths <- list_to_matrix(PSF_db_splicing)

# Expression

enrichment.parametric.Bin.all.Expression <- run_enrichment(MOFAgrouped.trained,
                                                           view = "Expression", factors = c(1,2,3),
                                                           feature.sets = t(Bin.paths),
                                                           sign = "all", set.statistic = "rank.sum",
                                                           # statistical.test = "permutation",
                                                           min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                           #nperm = 5000
)


enrichment.parametric.PS.all.Expression <- run_enrichment(MOFAgrouped.trained,
                                                          view = "Expression", factors = c(1,2,3),
                                                          feature.sets = t(PS.paths),
                                                          sign = "all", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                                          # statistical.test = "permutation",
                                                          min.size = 5, alpha = 0.1
                                                          #nperm = 5000
)

enrichment.parametric.PSF.all.Expression <- run_enrichment(MOFAgrouped.trained,
                                                           view = "Expression", factors = c(1,2,3),
                                                           feature.sets = t(PSF.paths),
                                                           sign = "all", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                           # statistical.test = "permutation",
                                                           min.size = 5, alpha = 0.1
                                                           #nperm = 5000
)

# IR

Bin_splicing.paths.IR <- Bin_splicing.paths[rownames(Bin_splicing.paths) %in% Event.annotation.split$EVENTID[Event.annotation.split$Type == "IR"],]
PS_splicing.paths.IR <- PS_splicing.paths[rownames(PS_splicing.paths) %in% Event.annotation.split$EVENTID[Event.annotation.split$Type == "IR"],]
PSF_splicing.paths.IR <- PSF_splicing.paths[rownames(PSF_splicing.paths) %in% Event.annotation.split$EVENTID[Event.annotation.split$Type == "IR"],]

enrichment.parametric.Bin.all.IR <- run_enrichment(MOFAgrouped.trained,
                                                   view = "IR", factors = c(1,2,3),
                                                   feature.sets = t(Bin_splicing.paths.IR),
                                                   sign = "all", set.statistic = "rank.sum",
                                                   # statistical.test = "permutation",
                                                   min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                   #nperm = 5000
)


enrichment.parametric.PS.all.IR <- run_enrichment(MOFAgrouped.trained,
                                                  view = "IR", factors = c(1,2,3),
                                                  feature.sets = t(PS_splicing.paths.IR),
                                                  sign = "all", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                                  # statistical.test = "permutation",
                                                  min.size = 5, alpha = 0.1
                                                  #nperm = 5000
)

enrichment.parametric.PSF.all.IR <- run_enrichment(MOFAgrouped.trained,
                                                   view = "IR", factors = c(1,2,3),
                                                   feature.sets = t(PSF_splicing.paths.IR),
                                                   sign = "all", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                   # statistical.test = "permutation",
                                                   min.size = 5, alpha = 0.1
                                                   #nperm = 5000
)

# ES

Bin_splicing.paths.ES <- Bin_splicing.paths[rownames(Bin_splicing.paths) %in% Event.annotation.split$EVENTID[Event.annotation.split$Type == "ES"],]
PS_splicing.paths.ES <- PS_splicing.paths[rownames(PS_splicing.paths) %in% Event.annotation.split$EVENTID[Event.annotation.split$Type == "ES"],]
PSF_splicing.paths.ES <- PSF_splicing.paths[rownames(PSF_splicing.paths) %in% Event.annotation.split$EVENTID[Event.annotation.split$Type == "ES"],]

enrichment.parametric.Bin.all.ES <- run_enrichment(MOFAgrouped.trained,
                                                   view = "ES", factors = c(1,2,3),
                                                   feature.sets = t(Bin_splicing.paths.ES),
                                                   sign = "all", set.statistic = "rank.sum",
                                                   # statistical.test = "permutation",
                                                   min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                   #nperm = 5000
)


enrichment.parametric.PS.all.ES <- run_enrichment(MOFAgrouped.trained,
                                                  view = "ES", factors = c(1,2,3),
                                                  feature.sets = t(PS_splicing.paths.ES),
                                                  sign = "all", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                                  # statistical.test = "permutation",
                                                  min.size = 5, alpha = 0.1
                                                  #nperm = 5000
)

enrichment.parametric.PSF.all.ES <- run_enrichment(MOFAgrouped.trained,
                                                   view = "ES", factors = c(1,2,3),
                                                   feature.sets = t(PSF_splicing.paths.ES),
                                                   sign = "all", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                   # statistical.test = "permutation",
                                                   min.size = 5, alpha = 0.1
                                                   #nperm = 5000
)

# AltAD

Bin_splicing.paths.AltAD <- Bin_splicing.paths[rownames(Bin_splicing.paths) %in% Event.annotation.split$EVENTID[Event.annotation.split$Type == "AltAD"],]
PS_splicing.paths.AltAD <- PS_splicing.paths[rownames(PS_splicing.paths) %in% Event.annotation.split$EVENTID[Event.annotation.split$Type == "AltAD"],]
PSF_splicing.paths.AltAD <- PSF_splicing.paths[rownames(PSF_splicing.paths) %in% Event.annotation.split$EVENTID[Event.annotation.split$Type == "AltAD"],]

enrichment.parametric.Bin.all.AltAD <- run_enrichment(MOFAgrouped.trained,
                                                      view = "AltAD", factors = c(1,2,3),
                                                      feature.sets = t(Bin_splicing.paths.AltAD),
                                                      sign = "all", set.statistic = "rank.sum",
                                                      # statistical.test = "permutation",
                                                      min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                      #nperm = 5000
)


enrichment.parametric.PS.all.AltAD <- run_enrichment(MOFAgrouped.trained,
                                                     view = "AltAD", factors = c(1,2,3),
                                                     feature.sets = t(PS_splicing.paths.AltAD),
                                                     sign = "all", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                                     # statistical.test = "permutation",
                                                     min.size = 5, alpha = 0.1
                                                     #nperm = 5000
)

enrichment.parametric.PSF.all.AltAD <- run_enrichment(MOFAgrouped.trained,
                                                      view = "AltAD", factors = c(1,2,3),
                                                      feature.sets = t(PSF_splicing.paths.AltAD),
                                                      sign = "all", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                      # statistical.test = "permutation",
                                                      min.size = 5, alpha = 0.1
                                                      #nperm = 5000
)

# AS

enrichment.parametric.Bin.all.AS <- run_enrichment(MOFAgrouped.trained,
                                                   view = "AS", factors = c(1,2,3),
                                                   feature.sets = t(Bin_splicing.paths),
                                                   sign = "all", set.statistic = "rank.sum",
                                                   # statistical.test = "permutation",
                                                   min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                   #nperm = 5000
)


enrichment.parametric.PS.all.AS <- run_enrichment(MOFAgrouped.trained,
                                                  view = "AS", factors = c(1,5),
                                                  feature.sets = t(PS_splicing.paths),
                                                  sign = "all", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                                  # statistical.test = "permutation",
                                                  min.size = 5, alpha = 0.1
                                                  #nperm = 5000
)

enrichment.parametric.PSF.all.AS <- run_enrichment(MOFAgrouped.trained,
                                                   view = "AS", factors = c(1,5),
                                                   feature.sets = t(PSF_splicing.paths),
                                                   sign = "all", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                   # statistical.test = "permutation",
                                                   min.size = 5, alpha = 0.1
                                                   #nperm = 5000
)

# Splicing

enrichment.parametric.Bin.all.Splicing <- run_enrichment(MOFAgrouped.trained,
                                                         view = "Splicing", factors = c(1,3),
                                                         feature.sets = t(Bin_splicing.paths),
                                                         sign = "all", set.statistic = "rank.sum",
                                                         # statistical.test = "permutation",
                                                         min.size = 5, p.adj.method = "fdr", alpha = 0.1
                                                         #nperm = 5000
)


enrichment.parametric.PS.all.Splicing <- run_enrichment(MOFAgrouped.trained,
                                                        view = "Splicing", factors = c(1,3),
                                                        feature.sets = t(PS_splicing.paths),
                                                        sign = "all", p.adj.method = "fdr",  set.statistic = "rank.sum",
                                                        # statistical.test = "permutation",
                                                        min.size = 5, alpha = 0.1
                                                        #nperm = 5000
)

enrichment.parametric.PSF.all.Splicing <- run_enrichment(MOFAgrouped.trained,
                                                         view = "Splicing", factors = c(1,3),
                                                         feature.sets = t(PSF_splicing.paths),
                                                         sign = "all", p.adj.method = "fdr", set.statistic = "rank.sum",
                                                         # statistical.test = "permutation",
                                                         min.size = 5, alpha = 0.1
                                                         #nperm = 5000
)

## plot enrich

# expression

order.bin.exp <- sort(rownames(enrichment.parametric.Bin.all.Expression$pval.adj))

Bin.df.exp <- t(enrichment.parametric.Bin.all.Expression$pval.adj)[,order.bin.exp]

for(i in 1:nrow(Bin.df.exp)){
   for(j in 1:ncol(Bin.df.exp)){
      if(Bin.df.exp[i,j] > 0.1 ){
         Bin.df.exp[i,j] <- 1
      }
   }
}

Bin.df.exp <- -log10(Bin.df.exp)

for(i in 1:nrow(Bin.df.exp)){
   for(j in 1:ncol(Bin.df.exp)){
      if(Bin.df.exp[i,j] > 10 ){
         Bin.df.exp[i,j] <- 10
      }
   }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Greens"))(500)

ht.Bin.exp <- Heatmap(Bin.df.exp, name = "Bins.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "white", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

order.ps.exp <- sort(rownames(enrichment.parametric.PS.all.Expression$pval.adj))

PS.df.exp <- t(enrichment.parametric.PS.all.Expression$pval.adj)[,order.ps.exp]

for(i in 1:nrow(PS.df.exp)){
   for(j in 1:ncol(PS.df.exp)){
      if(PS.df.exp[i,j] > 0.1 ){
         PS.df.exp[i,j] <- 1
      }
   }
}

PS.df.exp <- -log10(PS.df.exp)

for(i in 1:nrow(PS.df.exp)){
   for(j in 1:ncol(PS.df.exp)){
      if(PS.df.exp[i,j] > 10 ){
         PS.df.exp[i,j] <- 10
      }
   }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Purples"))(100)

ht.PS.exp <- Heatmap(t(PS.df.exp), name = "PS.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

order.psf.exp <- sort(rownames(enrichment.parametric.PSF.all.Expression$pval.adj))

PSF.df.exp <- t(enrichment.parametric.PSF.all.Expression$pval.adj)[,order.psf.exp]

for(i in 1:nrow(PSF.df.exp)){
   for(j in 1:ncol(PSF.df.exp)){
      if(PSF.df.exp[i,j] > 0.1 ){
         PSF.df.exp[i,j] <- 1
      }
   }
}

PSF.df.exp <- -log10(PSF.df.exp)

for(i in 1:nrow(PSF.df.exp)){
   for(j in 1:ncol(PSF.df.exp)){
      if(PSF.df.exp[i,j] > 10 ){
         PSF.df.exp[i,j] <- 10
      }
   }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100)

ht.PSF.exp <- Heatmap(t(PSF.df.exp), name = "PSF.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

draw(ht.Bin.exp + ht.PS.exp + ht.PSF.exp, ht_gap = unit(1, "cm"))
draw(ht.Bin.exp %v% ht.PS.exp %v% ht.PSF.exp, ht_gap = unit(1, "cm"))

# splicing

order.bin.exp <- sort(rownames(enrichment.parametric.Bin.all.Expression$pval.adj))

Bin.df.IR <- data.frame(Bin = order.bin.exp, Factor1 = NA, Factor2 = NA, Factor3 = NA)

for(i in 1:nrow(Bin.df.IR)){
   hit <- which(rownames(enrichment.parametric.Bin.all.IR$pval.adj) == Bin.df.IR$Bin[i])
   if(length(hit) == 0){
      Bin.df.IR$Factor1[i] <- 1
      Bin.df.IR$Factor2[i] <- 1
      Bin.df.IR$Factor3[i] <- 1
   }else if(length(hit) > 0){
      Bin.df.IR$Factor1[i] <- enrichment.parametric.Bin.all.IR$pval.adj[hit,1]
      Bin.df.IR$Factor2[i] <- enrichment.parametric.Bin.all.IR$pval.adj[hit,2]
      Bin.df.IR$Factor3[i] <- enrichment.parametric.Bin.all.IR$pval.adj[hit,3]
   }
}

rownames(Bin.df.IR) <- Bin.df.IR$Bin
Bin.df.IR <- Bin.df.IR[,-1]
Bin.df.IR <- t(Bin.df.IR)[,order.bin.exp]

for(i in 1:nrow(Bin.df.IR)){
   for(j in 1:ncol(Bin.df.IR)){
      if(Bin.df.IR[i,j] > 0.1 ){
         Bin.df.IR[i,j] <- 1
      }
   }
}

Bin.df.IR <- -log10(Bin.df.IR)

for(i in 1:nrow(Bin.df.IR)){
   for(j in 1:ncol(Bin.df.IR)){
      if(Bin.df.IR[i,j] > 10 ){
         Bin.df.IR[i,j] <- 10
      }
   }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Greens"))(500)

ht.Bin.IR <- Heatmap(t(Bin.df.IR), name = "Bins.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

order.ps.exp <- sort(rownames(enrichment.parametric.PS.all.Expression$pval.adj))

PS.df.IR <- data.frame(PS = order.ps.exp, Factor1 = NA, Factor2 = NA, Factor3 = NA)

for(i in 1:nrow(PS.df.IR)){
   hit <- which(rownames(enrichment.parametric.PS.all.IR$pval.adj) == PS.df.IR$PS[i])
   if(length(hit) == 0){
      PS.df.IR$Factor1[i] <- 1
      PS.df.IR$Factor2[i] <- 1
      PS.df.IR$Factor3[i] <- 1
   }else if(length(hit) > 0){
      PS.df.IR$Factor1[i] <- enrichment.parametric.PS.all.IR$pval.adj[hit,1]
      PS.df.IR$Factor2[i] <- enrichment.parametric.PS.all.IR$pval.adj[hit,2]
      PS.df.IR$Factor3[i] <- enrichment.parametric.PS.all.IR$pval.adj[hit,3]
   }
}

rownames(PS.df.IR) <- PS.df.IR$PS
PS.df.IR <- PS.df.IR[,-1]
PS.df.IR <- t(PS.df.IR)[,order.ps.exp]

for(i in 1:nrow(PS.df.IR)){
   for(j in 1:ncol(PS.df.IR)){
      if(PS.df.IR[i,j] > 0.1 ){
         PS.df.IR[i,j] <- 1
      }
   }
}

PS.df.IR <- -log10(PS.df.IR)

for(i in 1:nrow(PS.df.IR)){
   for(j in 1:ncol(PS.df.IR)){
      if(PS.df.IR[i,j] > 10 ){
         PS.df.IR[i,j] <- 10
      }
   }
}

col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Purples"))(100)

ht.PS.IR <- Heatmap(PS.df.IR, name = "PS.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "black", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

order.psf.exp <- sort(rownames(enrichment.parametric.PSF.all.Expression$pval.adj))

PSF.df.IR <- data.frame(PSF = order.psf.exp, Factor1 = NA, Factor2 = NA, Factor3 = NA)

for(i in 1:nrow(PSF.df.IR)){
   hit <- which(rownames(enrichment.parametric.PSF.all.IR$pval.adj) == PSF.df.IR$PSF[i])
   if(length(hit) == 0){
      PSF.df.IR$Factor1[i] <- 1
      PSF.df.IR$Factor2[i] <- 1
      PSF.df.IR$Factor3[i] <- 1
   }else if(length(hit) > 0){
      PSF.df.IR$Factor1[i] <- enrichment.parametric.PSF.all.IR$pval.adj[hit,1]
      PSF.df.IR$Factor2[i] <- enrichment.parametric.PSF.all.IR$pval.adj[hit,2]
      PSF.df.IR$Factor3[i] <- enrichment.parametric.PSF.all.IR$pval.adj[hit,3]
   }
}

rownames(PSF.df.IR) <- PSF.df.IR$PSF
PSF.df.IR <- PSF.df.IR[,-1]
PSF.df.IR <- t(PSF.df.IR)[,order.psf.exp]

for(i in 1:nrow(PSF.df.IR)){
   for(j in 1:ncol(PSF.df.IR)){
      if(PSF.df.IR[i,j] > 0.1 ){
         PSF.df.IR[i,j] <- 1
      }
   }
}

PSF.df.IR <- -log10(PSF.df.IR)

for(i in 1:nrow(PSF.df.IR)){
   for(j in 1:ncol(PSF.df.IR)){
      if(PSF.df.IR[i,j] > 10 ){
         PSF.df.IR[i,j] <- 10
      }
   }
}
col_quantity <- colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100)

ht.PSF.IR <- Heatmap(PSF.df.IR, name = "PSF.-log10(FDR)", col = col_quantity, rect_gp = gpar(col = "white", lwd = 2), show_column_dend = F, show_column_names = T, cluster_columns = F, cluster_rows = F, show_row_names = T,
)

draw(ht.Bin.IR + ht.PS.exp + ht.PSF.exp, ht_gap = unit(1, "cm"))
draw(ht.Bin.IR %v% ht.PS.exp %v% ht.PSF.exp, ht_gap = unit(1, "cm"))

#### 5. Plot barplots of isoform expression balance from validated events ####

EventNormCounts <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/6.Transcriptional_Analyses/Tables/Splicing/Batch_Corrected/WithGroup/0.countsEvents_lowFilter_PEandTE_batchCorrectedTechStudy_Normalized.txt")
View(EventNormCounts)

# ERD6

ERD6.count <- EventNormCounts[EventNormCounts$events.names == "bcc_456625|Cycle_78",]
ERD6.count <- ERD6.count[,-1]
ERD6.count <- ERD6.count[!is.na(ERD6.count),]
ERD6.count <- t(ERD6.count[1:2,])
ERD6.count <- ERD6.count[-1,]
ERD6.count <- ERD6.count[c(1:16, 18:42), ]


ERD6 <- data.frame(Tissues = c(rep("Bud", 2),rep("Needle", 2)),
                   Isoform = c("Long", "Short", "Long", "Short"),
                   Expression = c(mean(ERD6.count[1:16,1])/sum(mean(ERD6.count[1:16,1]), mean(ERD6.count[1:16,2])), 
                                  mean(ERD6.count[1:16,2])/sum(mean(ERD6.count[1:16,1]), mean(ERD6.count[1:16,2])),
                                  mean(ERD6.count[17:41,1])/sum(mean(ERD6.count[17:41,1]), mean(ERD6.count[17:41,2])), 
                                  mean(ERD6.count[17:41,2])/sum(mean(ERD6.count[17:41,1]), mean(ERD6.count[17:41,2]))
                                  )
                   )

pal_annotation <- c("lightsteelblue", "lightgoldenrod1")
plot1 <- ggplot(ERD6,aes(x=Tissues,y=Expression,fill=Isoform,color=Events,width=0.75),size=0.01)+
   geom_bar(stat="identity",color="black",size=0.1)+theme_classic()+labs(x=NULL,y="Expression")+
   scale_fill_manual(values=pal_annotation)+scale_colour_manual(values=rep("#000000",11))+
   theme(legend.key = element_rect(color="gray",size=0.1),legend.key.size=unit(0.2,"cm"),
         legend.position="right",legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.3, unit="cm"),
         axis.text.x = element_text(color="black",size=5,face="bold"),axis.text.y = element_text(color="black",size=5.5),
         axis.title.x = element_text(size=6,face="bold"),axis.title.y = element_text(size=6.5,face="bold"),
         legend.title = element_text(size=6.5,face="bold"),legend.text = element_text(size=5.5),
         panel.border = element_rect(colour = "black", fill=NA, size=0.2,linetype=3))

# ELF4

ELF4.count <- EventNormCounts[EventNormCounts$events.names == "bcc_668500|Cycle_4",]
ELF4.count <- ELF4.count[,-1]
ELF4.count <- ELF4.count[!is.na(ELF4.count),]
ELF4.count <- t(ELF4.count[1:2,])
ELF4.count <- ELF4.count[-1,]
ELF4.count <- ELF4.count[c(1:16, 18:42), ]

ELF4 <- data.frame(Tissues = c(rep("Bud", 2),rep("Needle", 2)),
                   Isoform = c("Long", "Short", "Long", "Short"),
                   Expression = c(mean(ELF4.count[1:16,1])/sum(mean(ELF4.count[1:16,1]), mean(ELF4.count[1:16,2])), 
                                  mean(ELF4.count[1:16,2])/sum(mean(ELF4.count[1:16,1]), mean(ELF4.count[1:16,2])),
                                  mean(ELF4.count[17:41,1])/sum(mean(ELF4.count[17:41,1]), mean(ELF4.count[17:41,2])), 
                                  mean(ELF4.count[17:41,2])/sum(mean(ELF4.count[17:41,1]), mean(ELF4.count[17:41,2]))
                   )
)

pal_annotation <- c("lightsteelblue", "lightgoldenrod1")
plot2 <- ggplot(ELF4,aes(x=Tissues,y=Expression,fill=Isoform,color=Events,width=0.75),size=0.01)+
   geom_bar(stat="identity",color="black",size=0.1)+theme_classic()+labs(x=NULL,y="Expression")+
   scale_fill_manual(values=pal_annotation)+scale_colour_manual(values=rep("#000000",11))+
   theme(legend.key = element_rect(color="gray",size=0.1),legend.key.size=unit(0.2,"cm"),
         legend.position="right",legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.3, unit="cm"),
         axis.text.x = element_text(color="black",size=5,face="bold"),axis.text.y = element_text(color="black",size=5.5),
         axis.title.x = element_text(size=6,face="bold"),axis.title.y = element_text(size=6.5,face="bold"),
         legend.title = element_text(size=6.5,face="bold"),legend.text = element_text(size=5.5),
         panel.border = element_rect(colour = "black", fill=NA, size=0.2,linetype=3))

# cpSRP43

cpSRP43.count <- EventNormCounts[EventNormCounts$events.names == "bcc_629945|Cycle_2594791",]
cpSRP43.count <- cpSRP43.count[,-1]
cpSRP43.count <- cpSRP43.count[!is.na(cpSRP43.count),]
cpSRP43.count <- t(cpSRP43.count[1:2,])
cpSRP43.count <- cpSRP43.count[-1,]
cpSRP43.count <- cpSRP43.count[c(1:16, 18:42), ]

cpSRP43 <- data.frame(Tissues = c(rep("Bud", 2),rep("Needle", 2)),
                   Isoform = c("Long", "Short", "Long", "Short"),
                   Expression = c(mean(cpSRP43.count[1:16,1])/sum(mean(cpSRP43.count[1:16,1]), mean(cpSRP43.count[1:16,2])), 
                                  mean(cpSRP43.count[1:16,2])/sum(mean(cpSRP43.count[1:16,1]), mean(cpSRP43.count[1:16,2])),
                                  mean(cpSRP43.count[17:41,1])/sum(mean(cpSRP43.count[17:41,1]), mean(cpSRP43.count[17:41,2])), 
                                  mean(cpSRP43.count[17:41,2])/sum(mean(cpSRP43.count[17:41,1]), mean(cpSRP43.count[17:41,2]))
                   )
)

pal_annotation <- c("lightsteelblue", "lightgoldenrod1")
plot3 <- ggplot(cpSRP43,aes(x=Tissues,y=Expression,fill=Isoform,color=Events,width=0.75),size=0.01)+
   geom_bar(stat="identity",color="black",size=0.1)+theme_classic()+labs(x=NULL,y="Expression")+
   scale_fill_manual(values=pal_annotation)+scale_colour_manual(values=rep("#000000",11))+
   theme(legend.key = element_rect(color="gray",size=0.1),legend.key.size=unit(0.2,"cm"),
         legend.position="right",legend.margin=margin(t=-0.1, r=-0.1, b=-0.1, l=-0.3, unit="cm"),
         axis.text.x = element_text(color="black",size=5,face="bold"),axis.text.y = element_text(color="black",size=5.5),
         axis.title.x = element_text(size=6,face="bold"),axis.title.y = element_text(size=6.5,face="bold"),
         legend.title = element_text(size=6.5,face="bold"),legend.text = element_text(size=5.5),
         panel.border = element_rect(colour = "black", fill=NA, size=0.2,linetype=3))
