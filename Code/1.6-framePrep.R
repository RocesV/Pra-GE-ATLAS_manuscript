#### 6-frame proteomics prepare Pra-GE-ATLAS_Pra

### Load libraries

library(Biostrings)
library(stringr)
library(foreach)
library(parallel)
library(doParallel)
library(seqRFLP)
library(doSNOW)
library(progress)

### Program

setwd("G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/1.Assemblies/2.EvidentialGene/final-repeat/0.send/")
WorkingDirectory <- "G:/Pra-GE-ATLAS/ATLAS/Pinus_radiata/1.Assemblies/2.EvidentialGene/final-repeat/0.send/"

for(i in list.files(WorkingDirectory)){
    all <- readAAStringSet(filepath = paste0(WorkingDirectory,i))
    all_splitted <- list()
    
    ## 0. splitting peptides by X  
    # parallel fashion
    n.cores <- parallel::detectCores()-2
    my.cluster <- parallel::makeCluster(n.cores)
    doSNOW::registerDoSNOW(my.cluster)
    parallel::clusterExport(cl = my.cluster, varlist = c("all"))
    all_splitted <- foreach::foreach(k = 1:length(all), 
                               .combine = rbind, .verbose = T) %dopar%
      {
        pep <- strsplit(as.character(all[[k]]), split="X", fixed = T)
        df <- cbind(paste0(strsplit(names(all[k]), split = " ", fixed= T)[[1]][1],"|",1:length(unlist(pep))), do.call(cbind, pep))
        return(df)
      }
    parallel::stopCluster(cl = my.cluster)
    # sequential
    pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                           total = length(all),
                           complete = "=",   # Completion bar character
                           incomplete = "-", # Incomplete bar character
                           current = ">",    # Current bar character
                           clear = FALSE,    # If TRUE, clears the bar when finish
                           width = 100)      # Width of the progress bar
    for(k in 1:length(all)){
      pb$tick()
      all_splitted[[k]] <- strsplit(as.character(all[[k]]), split="X", fixed = T)
      all_splitted[[k]] <- cbind(paste0(strsplit(names(all[k]), split = " ", fixed= T)[[1]][1],"|",1:length(unlist(all_splitted[[k]]))), do.call(cbind, all_splitted[[k]]))
    }
    
    ## 1.format and new columns for easy use
    if(!is.data.frame(all_splitted)){ all_splitted_f <- do.call(rbind, all_splitted)}
    colnames(all_splitted_f) <- c("ID","Sequence")
    frame <- str_sub(string = do.call(rbind,strsplit(x = all_splitted_f[,1], split="|", fixed =T))[,1], start = -2)
    frame <- gsub(pattern = "_1", replacement = "F", x = frame)
    frame <- gsub(pattern = "_2", replacement = "F", x = frame)
    frame <- gsub(pattern = "_3", replacement = "F", x = frame)
    frame <- gsub(pattern = "_4", replacement = "R", x = frame)
    frame <- gsub(pattern = "_5", replacement = "R", x = frame)
    frame <- gsub(pattern = "_6", replacement = "R", x = frame)
    all_splitted_f <- as.data.frame(all_splitted_f)
    all_splitted_f$frame <- frame
    sequence <- paste0(do.call(rbind,strsplit(x = all_splitted_f[,1], split="_", fixed =T))[,1],"_",do.call(rbind,strsplit(x = all_splitted_f[,1], split="_", fixed =T))[,2],"_",
                       do.call(rbind,strsplit(x = all_splitted_f[,1], split="_", fixed =T))[,3],"_",do.call(rbind,strsplit(x = all_splitted_f[,1], split="_", fixed =T))[,4],"_",
                       do.call(rbind,strsplit(x = all_splitted_f[,1], split="_", fixed =T))[,5])
    sequence <- paste0(do.call(rbind,strsplit(x = aa50[,1], split="_", fixed =T))[,1],"_",do.call(rbind,strsplit(x = aa50[,1], split="_", fixed =T))[,2],"_",
                       do.call(rbind,strsplit(x = aa50[,1], split="_", fixed =T))[,3],"_",do.call(rbind,strsplit(x = aa50[,1], split="_", fixed =T))[,4],"_",
                       do.call(rbind,strsplit(x = aa50[,1], split="_", fixed =T))[,5])
    aa50$IDsequence <- sequence
    sequence <- paste0(do.call(rbind,strsplit(all_splitted_f[,1], split="_"))[,1])
    all_splitted_f$IDsequence <- sequence
    length <- nchar(all_splitted_f$Sequence)
    all_splitted_f$length <- length
    rm(frame)
    rm(length)
    rm(sequence)
    gc()
    
    ## 2.Filter by aa length (25,50,75)
    aa25 <- all_splitted_f[-which(all_splitted_f$length < 25),1:2]
    aa50 <- all_splitted_f[-which(all_splitted_f$length < 50),1:2]
    aa75 <- all_splitted_f[-which(all_splitted_f$length < 75),1:2]
    seqinr::write.fasta(sequences = as.list(aa25$Sequence), names = aa25$ID, 
                        file.out = paste0(WorkingDirectory,i,"/3.25aa/",i,"_aa25.aa"),
                        as.string = T, open = "w")
    seqinr::write.fasta(sequences = as.list(aa50$Sequence), names = aa50$ID, 
                        file.out = paste0(WorkingDirectory,i,"/4.50aa/",i,"_aa50.aa"),
                        as.string = T, open = "w")
    seqinr::write.fasta(sequences = as.list(aa75$Sequence), names = aa75$ID, 
                        file.out = paste0(WorkingDirectory,i,"/5.75aa/",i,"_aa75.aa"),
                        as.string = T, open = "w")
    rm(aa25)
    rm(aa75)
    rm(aa50)
    gc()
    
   ## 3.Filter by aa length and select longest for each sense &| sequence
   aa50 <- all_splitted_f[-which(all_splitted_f$length < 50),]
   transcriptsaa50 <- levels(as.factor(aa50$IDsequence))
   # sequential
   a <- Sys.time()
   Longestbyframe <- list()
   pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                          total = length(transcriptsaa50),
                          complete = "=",   # Completion bar character
                          incomplete = "-", # Incomplete bar character
                          current = ">",    # Current bar character
                          clear = FALSE,    # If TRUE, clears the bar when finish
                          width = 100)      # Width of the progress bar
   for(k in 1:length(transcriptsaa50)){
     pb$tick()
     rows.F <- which(aa50$IDsequence == transcriptsaa50[k] & aa50$frame == "F")
     rows.R <- which(aa50$IDsequence == transcriptsaa50[k] & aa50$frame == "R")
     max.rows.F <- aa50[rows.F[which.max(aa50$length[rows.F])],]
     max.rows.R <- aa50[rows.R[which.max(aa50$length[rows.R])],]
     Longestbyframe[[k]] <- rbind(max.rows.F, max.rows.R)
   } 
   longestbyframee <- do.call(rbind, Longestbyframe)
   onlyoneframe.rows <- longestbyframee$IDsequence %in% names(table(longestbyframee$IDsequence)[which(table(longestbyframee$IDsequence) == 1)])
   onlyoneframe <- longestbyframee[onlyoneframe.rows,]
   twoframe <- longestbyframee[!onlyoneframe.rows,]
   twoframe.F <- twoframe[which(twoframe$frame == "F"),]
   twoframe.R <- twoframe[which(twoframe$frame == "R"),]
   if(!length(which((twoframe.F$IDsequence == twoframe.R$IDsequence) ==TRUE)) == nrow(twoframe.F)){print("\n Stop something went bad boy... \n")}
   longestt <- rbind(onlyoneframe,
                     twoframe.F[(twoframe.F$length > twoframe.R$length),],
                     twoframe.R[(twoframe.R$length > twoframe.F$length),],
                     twoframe.F[(twoframe.F$length == twoframe.R$length),],
                     twoframe.R[(twoframe.F$length == twoframe.R$length),])
   longestt <- longestt[order(longestt$ID),]
   b <- Sys.time()
   b-a
   seqinr::write.fasta(sequences = as.list(longestbyframee$Sequence), names = longestbyframee$ID, 
                       file.out = paste0(WorkingDirectory,i,"/6.Longestbysnseand50aa/",i,"_longestbysenseaa50.aa"),
                       as.string = T, open = "w")
   seqinr::write.fasta(sequences = as.list(longestt$Sequence), names = longestt$ID, 
                       file.out = paste0(WorkingDirectory,i,"/7.Longestand50aa/",i,"_longest50aa.aa"),
                       as.string = T, open = "w")
   rm(longestt)
   rm(longestbyframee)
   rm(Longestbyframe)
   rm(twoframe.F)
   rm(twoframe.R)
   rm(twoframe)
   rm(onlyoneframe)
   rm(onlyoneframe.rows)
   gc()
   
   # parallel 
   c <- Sys.time()
   n.cores <- parallel::detectCores()-2
   my.cluster <- parallel::makeCluster(n.cores)
   doSNOW::registerDoSNOW(my.cluster)
   parallel::clusterExport(cl = my.cluster, varlist = c("aa50", "transcriptsaa50"))
   longestbyframee <- foreach::foreach(k = 1:length(transcriptsaa50) , 
                               .combine = rbind, .verbose = T) %dopar%
     {
       rows.F <- which(aa50$IDsequence == transcriptsaa50[k]  & aa50$frame == "F")
       rows.R <- which(aa50$IDsequence == transcriptsaa50[k]  & aa50$frame == "R")
       max.rows.F <- aa50[rows.F[which.max(aa50$length[rows.F])],]
       max.rows.R <- aa50[rows.R[which.max(aa50$length[rows.R])],]
       max.rows <- rbind(max.rows.F, max.rows.R)
       return(max.rows)
     }
   parallel::stopCluster(cl = my.cluster)
   onlyoneframe.rows <- longestbyframee$IDsequence %in% names(table(longestbyframee$IDsequence)[which(table(longestbyframee$IDsequence) == 1)])
   onlyoneframe <- longestbyframee[onlyoneframe.rows,]
   twoframe <- longestbyframee[!onlyoneframe.rows,]
   twoframe.F <- twoframe[which(twoframe$frame == "F"),]
   twoframe.R <- twoframe[which(twoframe$frame == "R"),]
   if(!length(which((twoframe.F$IDsequence == twoframe.R$IDsequence) ==TRUE)) == nrow(twoframe.F)){print("\n Stop something went bad boy... \n")}
   longestt <- rbind(onlyoneframe,
                     twoframe.F[(twoframe.F$length > twoframe.R$length),],
                     twoframe.R[(twoframe.R$length > twoframe.F$length),],
                     twoframe.F[(twoframe.F$length == twoframe.R$length),],
                     twoframe.R[(twoframe.F$length == twoframe.R$length),])
   longestt <- longestt[order(longestt$ID),]
   d <- Sys.time()
   d-c
   seqinr::write.fasta(sequences = as.list(longestbyframee$Sequence), names = longestbyframee$ID, 
                       file.out = paste0(WorkingDirectory,i,"/6.Longestbysenseand50aa/",i,"_longestbysenseaa50.aa"),
                       as.string = T, open = "w")
   seqinr::write.fasta(sequences = as.list(longestt$Sequence), names = longestt$ID, 
                       file.out = paste0(WorkingDirectory,i,"/7.Longestand50aa/",i,"_longest50aa.aa"),
                       as.string = T, open = "w")
   rm(longestt)
   rm(longestbyframee)
   rm(twoframe.F)
   rm(twoframe.R)
   rm(twoframe)
   rm(onlyoneframe)
   rm(onlyoneframe.rows)
   gc()
  }
