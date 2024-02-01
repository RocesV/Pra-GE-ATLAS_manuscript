#####################################
# Pra-GE-ATLAS - Annotation Seeker  #
#####################################

###################################################################################################
# Annotation.Seeker == Parse all the annotations computed by Mercator, IP5, dammit and eggnog-mapper
###################################################################################################

#### 0.Load libraries and pkgs ####

suppressPackageStartupMessages(library(cowsay)) 
suppressPackageStartupMessages(library(multicolor)) 
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(progress))

#### 1.Dammit-Parser ####

cat("\n Importing GFF3 and NameMap... \n")

dammit.output <- "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/2.dammit/EviGene_okayall/EviGene_allokay_send_2/EviGene_allokay_send/Pra-GE-ATLAS_EviGene_okayall.cdna.dammit.gff3"
dammit.namemap <- "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/2.dammit/EviGene_okayall/EviGene_allokay_send_2/EviGene_allokay_send/Pra-GE-ATLAS_EviGene_okayall.cdna.dammit.namemap.csv"

dammit.GFF3 <- read.delim(dammit.output, header = F, sep = "\t")
dammit.GFF3 <- dammit.GFF3[-1,]
dammit.names <- read.csv(dammit.namemap, header = T, sep = ",")
dammit.GFF3.noTransdecoder <- dammit.GFF3[which(dammit.GFF3$V2 != "transdecoder"),]

identified <- length(dammit.names$renamed[which(dammit.names$renamed %in% dammit.GFF3.noTransdecoder$V1)])
total <- nrow(dammit.names)
percentage.identification <- round(identified / total * 100,2)
dammit.GFF3.noTransdecoder$V9 <- as.character(dammit.GFF3.noTransdecoder$V9)

cat(paste0("\n ", percentage.identification, " % sequences have been identified (", identified, ") \n"))

final.df <- as.data.frame(matrix(NA, nrow(dammit.names), ncol(dammit.names) + 7))
colnames(final.df) <- c("TranscriptID", "dammitID", "All","uniref", "=sprot","OrthoDB", "Pfam:", "Rfam:", "Epiomides")
final.df[,c(1,2)] <- dammit.names
rownames(final.df) <- final.df$dammitID
rm(dammit.names)
rm(dammit.GFF3)

## After awk scripts for collapsing same transcript info
dammit.table <- read.delim("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/2.dammit/EviGene_okayall/EviGene_allokay_send_2/EviGene_allokay_send/Pra-GE-ATLAS_EviGene_okayall.cdna.dammit.gff3.table", sep = " ")
final.df <- final.df[order(final.df$dammitID),]
dammit.table <- dammit.table[order(dammit.table$X..gff.version),]
final.df$All[final.df$dammitID %in% dammit.table$X..gff.version] <- dammit.table$X 

write.table(final.df[,c(1,2,3)],  "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/4.FinalAnnotation/Pra-GE-ATLAS_dammit.txt", col.names = T, sep = "\t", quote = F, row.names = F)

#### 2.Mercator-Parser ####

setwd("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/1.Mercator4/EviGene_okayall/Pra-GE-ATLAS_Pra_alloksReduced0905/")
path <- "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/1.Mercator4/EviGene_okayall/Pra-GE-ATLAS_Pra_alloksReduced0905/"

final.df <- final.df[,-c(4:9)]
colnames(final.df) <- c("TranscriptID", "dammitID", "dammit_All")
final.df$MercatorBinCo <- NA
final.df$MercatorBin <- NA
final.df$MercatorDescription <- NA
final.df$MercatorName <- NA

MercatorRes <- read.delim("Pra-GE-ATLAS0905Redo.results.txt", header = T, sep = "\t")
MercatorRes$IDENTIFIER <- gsub("pra", "Pra", x = MercatorRes$IDENTIFIER)
MercatorRes$IDENTIFIER <- gsub("ge_atlas_radiata", "GE_ATLAS_Radiata", x = MercatorRes$IDENTIFIER)
MercatorRes$IDENTIFIER <- gsub("'", "", x = MercatorRes$IDENTIFIER, fixed = T)
MercatorRes <- MercatorRes[order(MercatorRes$IDENTIFIER),]
MercatorRes <- MercatorRes[-which(MercatorRes$IDENTIFIER == ""),]
MercatorRes$BINCODE <- gsub("'", "", x = MercatorRes$BINCODE, fixed = T)
MercatorRes$NAME <- gsub(" ", "", x = MercatorRes$NAME, fixed = T)
MercatorRes$DESCRIPTION <- gsub(" ", "", x = MercatorRes$DESCRIPTION, fixed = T)
MercatorRes$BIN <- do.call(rbind,strsplit(x = MercatorRes$BINCODE, split = ".", fixed = T))[,1]
write.table(MercatorRes, "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/1.Mercator4/EviGene_okayall/Pra-GE-ATLAS_Pra_alloksReduced0905/Mercator_formated.txt", sep = "\t", col.names = T, quote = F, row.names = F)
## After awk script for collapsing same transcript info
MercatorRes_BinCo <- read.delim("Mercator_formated_Bins.txt", header = T, sep = " ")
MercatorRes_Bin <- read.delim("Mercator_formated_Bin.txt", header = T, sep = " ")
MercatorRes_Names <- read.delim("Mercator_formated_Name.txt", header = T, sep = " ")
MercatorRes_Description <- read.delim("Mercator_formated_Description.txt", header = T, sep = " ")
length(which((final.df$TranscriptID == MercatorRes_BinCo$IDENTIFIER)==TRUE)) # check
length(which((final.df$TranscriptID == MercatorRes_BinCo$IDENTIFIER)==TRUE)) # check
length(which((final.df$TranscriptID == MercatorRes_BinCo$IDENTIFIER)==TRUE)) # check
length(which((final.df$TranscriptID == MercatorRes_Bin$IDENTIFIER)==TRUE))
final.df$MercatorBinCo <- MercatorRes_BinCo$BINCODE
final.df$MercatorDescription <- MercatorRes_Description$DESCRIPTION
final.df$MercatorName <- MercatorRes_Names$NAME
final.df$MercatorBin <- MercatorRes_Bin$BIN
write.table(final.df,  "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/4.FinalAnnotation/Pra-GE-ATLAS_dammitMercator.txt", col.names = T, sep = "\t", quote = F, row.names = F)

#### 3.IP5-Parser ####

setwd("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/1.InterProScan5/EviGene_allokay/cat/")
interpro.IP5_SignatureAccesion <- read.delim("Pra-GE-ATLAS-Pra_alloksReduced0905_IP5_5.names.sorted.tsv", header = F, sep = "\t")
interpro.IP5_Accesion <- read.delim("Pra-GE-ATLAS-Pra_alloksReduced0905_IP5_12.names.sorted.tsv", header = F, sep = "\t")
interpro.IP5_Description <- read.delim("Pra-GE-ATLAS-Pra_alloksReduced0905_IP5_13.names.sorted.tsv", header = F, sep = "\t")
interpro.IP5_GO <- read.delim("Pra-GE-ATLAS-Pra_alloksReduced0905_IP5_14.names.sorted.tsv", header = F, sep = "\t")
final.df$IP5_SignatureAccesion <- NA
final.df$IP5_Accesion <- NA
final.df$IP5_Description <- NA
final.df$IP5_GO <- NA

length(which((interpro.IP5_Accesion$V1 %in% final.df$TranscriptID)== TRUE)) == nrow(interpro.IP5_Accesion)
final.df <- final.df[order(final.df$TranscriptID),]
interpro.IP5_Accesion <- interpro.IP5_Accesion[order(interpro.IP5_Accesion$V1),] 
final.df$IP5_Accesion[final.df$TranscriptID %in% interpro.IP5_Accesion$V1] <- interpro.IP5_Accesion$V2

length(which((interpro.IP5_Description$V1 %in% final.df$TranscriptID)== TRUE)) == nrow(interpro.IP5_Description)
final.df <- final.df[order(final.df$TranscriptID),]
interpro.IP5_Description <- interpro.IP5_Description[order(interpro.IP5_Description$V1),] 
final.df$IP5_Description[final.df$TranscriptID %in% interpro.IP5_Description$V1] <- interpro.IP5_Description$V2

length(which((interpro.IP5_GO$V1 %in% final.df$TranscriptID)== TRUE)) == nrow(interpro.IP5_GO)
final.df <- final.df[order(final.df$TranscriptID),]
interpro.IP5_GO <- interpro.IP5_GO[order(interpro.IP5_GO$V1),] 
final.df$IP5_GO[final.df$TranscriptID %in% interpro.IP5_GO$V1] <- interpro.IP5_GO$V2

length(which((interpro.IP5_SignatureAccesion$V1 %in% final.df$TranscriptID)== TRUE)) == nrow(interpro.IP5_SignatureAccesion)
final.df <- final.df[order(final.df$TranscriptID),]
interpro.IP5_SignatureAccesion <- interpro.IP5_SignatureAccesion[order(interpro.IP5_SignatureAccesion$V1),] 
final.df$IP5_SignatureAccesion[final.df$TranscriptID %in% interpro.IP5_SignatureAccesion$V1] <- interpro.IP5_SignatureAccesion$V2

write.table(final.df,  "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/4.FinalAnnotation/Pra-GE-ATLAS_dammitMercatorIP5.txt", col.names = T, sep = "\t", quote = F, row.names = F)

#### 4.Eggnog-mapper-Parser ####

setwd("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/3.Eggnog-mapper/Pra-GE-ATLAS_pepsAll_done/cat/")
EGGNOGMAPPER.NAME9 <- read.delim("Pra-GE-ATLAS-Pra_alloksReduced0905_eggnogmapper_NAME9.names.sorted.annotations", header = F, sep = "\t")
EGGNOGMAPPER.SYMBOL10 <- read.delim("Pra-GE-ATLAS-Pra_alloksReduced0905_eggnogmapper_SYMBOL10.names.sorted.annotations", header = F, sep = "\t")
EGGNOGMAPPER.GO11 <- read.delim("Pra-GE-ATLAS-Pra_alloksReduced0905_eggnogmapper_GO11.names.sorted.annotations", header = F, sep = "\t")
EGGNOGMAPPER.DESC22 <- read.delim("Pra-GE-ATLAS-Pra_alloksReduced0905_eggnogmapper_DESC22.names.sorted.annotations", header = F, sep = "\t")
final.df$EGGNOGMAPPER_NAME <- NA
final.df$EGGNOGMAPPER_SYMBOL <- NA
final.df$EGGNOGMAPPER_GO <- NA
final.df$EGGNOGMAPPER_DESC <- NA

length(which((EGGNOGMAPPER.NAME9$V1 %in% final.df$TranscriptID)== TRUE)) == nrow(EGGNOGMAPPER.NAME9)
final.df <- final.df[order(final.df$TranscriptID),]
EGGNOGMAPPER.NAME9 <- EGGNOGMAPPER.NAME9[order(EGGNOGMAPPER.NAME9$V1),] 
final.df$EGGNOGMAPPER_NAME[final.df$TranscriptID %in% EGGNOGMAPPER.NAME9$V1] <- EGGNOGMAPPER.NAME9$V2

length(which((EGGNOGMAPPER.SYMBOL10$V1 %in% final.df$TranscriptID)== TRUE)) == nrow(EGGNOGMAPPER.SYMBOL10)
final.df <- final.df[order(final.df$TranscriptID),]
EGGNOGMAPPER.SYMBOL10 <- EGGNOGMAPPER.SYMBOL10[order(EGGNOGMAPPER.SYMBOL10$V1),] 
final.df$EGGNOGMAPPER_SYMBOL[final.df$TranscriptID %in% EGGNOGMAPPER.SYMBOL10$V1] <- EGGNOGMAPPER.SYMBOL10$V2

length(which((EGGNOGMAPPER.GO11$V1 %in% final.df$TranscriptID)== TRUE)) == nrow(EGGNOGMAPPER.GO11)
final.df <- final.df[order(final.df$TranscriptID),]
EGGNOGMAPPER.GO11 <- EGGNOGMAPPER.GO11[order(EGGNOGMAPPER.GO11$V1),] 
final.df$EGGNOGMAPPER_GO[final.df$TranscriptID %in% EGGNOGMAPPER.GO11$V1] <- EGGNOGMAPPER.GO11$V2

length(which((EGGNOGMAPPER.DESC22$V1 %in% final.df$TranscriptID)== TRUE)) == nrow(EGGNOGMAPPER.DESC22)
final.df <- final.df[order(final.df$TranscriptID),]
EGGNOGMAPPER.DESC22 <- EGGNOGMAPPER.DESC22[order(EGGNOGMAPPER.DESC22$V1),] 
final.df$EGGNOGMAPPER_DESC[final.df$TranscriptID %in% EGGNOGMAPPER.DESC22$V1] <- EGGNOGMAPPER.DESC22$V2

write.table(final.df,  "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/4.FinalAnnotation/Pra-GE-ATLAS_dammitMercatorIP5Eggnog.txt", col.names = T, sep = "\t", quote = F, row.names = F)

#### 5.ProteinHit-Parser ####

setwd("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/4.FinalAnnotation/blast_output/")
ProtHit <- read.delim("Prots-ATLAS_correspondence.sorted.reduced.tsv", header = F, sep = "\t")
final.df$ProteinID <- NA

length(which((ProtHit$V1 %in% final.df$TranscriptID)== TRUE)) == nrow(ProtHit)
final.df <- final.df[order(final.df$TranscriptID),]
ProtHit <- ProtHit[order(ProtHit$V1),] 
final.df$ProteinID[final.df$TranscriptID %in% ProtHit$V1] <- ProtHit$V2

write.table(final.df,  "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/4.FinalAnnotation/Pra-GE-ATLAS_dammitMercatorIP5EggnogProtID.txt", col.names = T, sep = "\t", quote = F, row.names = F)

#### 6.Reduce Final Table Annotation for downstream analyses ####

final.df.Reduced <- final.df[,c(1,2,16,5,14,15,13,12)]
colnames(final.df.Reduced)
write.table(final.df.Reduced,  "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/3.Annotation/4.FinalAnnotation/Pra-GE-ATLAS_ReducedAnaltyses.txt", col.names = T, sep = "\t", quote = F, row.names = F)
