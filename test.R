setwd("~/Libs/projects/drugScreen/")
dir()

### 0. get the drug's positions on the Well & Plate
pw <- read.table("Prestwick_384.tab",header = F, sep = "\t",quote = ""
                 ,colClasses = "character")
colnames(pw) <- c("drugName")
pw$Well <- paste(rep(LETTERS[1:16],each=24),rep(sprintf("%02d",1:24),16),sep="")
pw$Plate <- rep(sprintf("%02d",1:4),each=384)
pw$Library <- "PW"
pw$code <- paste(pw$Library,pw$Plate,pw$Well,sep="")
rownames(pw) <- as.character(pw$code)
# write.table(pw[,c(1:4)],"drugList.standard_format.txt",quote = F, sep = "\t", col.names = T, row.names = F)

### 1. get the drug's name in raw drug screen data
htc.be007 <- read.table("MMC007_BE_CRD_PW_1uM_RAW_DATA.tab"
                        ,header = T
                        ,sep = "\t"
                        ,quote = ""
                        ,colClasses = c(rep("character",7),rep("numeric",2)))
head(htc.be007)
htc.be007$Plate.new <- as.character(sprintf("%02d",as.numeric(sub("[AB]","",htc.be007$Plate,perl=T))))
htc.be007$code <- paste(htc.be007$Library,htc.be007$Plate.new,htc.be007$Well,sep="")
htc.be007$drugName <- pw[htc.be007$code,"drugName"]
htc.be007$Sample <- paste(htc.be007$Disease,htc.be007$Patient
                          ,htc.be007$Library,htc.be007$Dose
                          ,sub("[0-9]+","rep",htc.be007$Plate,perl=T)
                          ,sep="_")
# write.table(htc.be007[,c(1:9)], "rawData.standard_format.txt", quote = F, sep = "\t",col.names = T,row.names = F)

### 2. output raw data for each drug
pw.data <- pw
sample.list <- paste(rep(unique(htc.be007$Sample),2)
                     ,rep(c('Area','Count'),each=4)
                     ,sep=".")
for ( i in sample.list){
    pw.data[[i]] <- 0
}
for ( i in 1:nrow(pw.data)){
    for ( j in sample.list){
        tmp.code <- pw.data[i,"code"]
        tmp.sample <- unlist(strsplit(j,"[.]",perl=T))[1]
        tmp.subject <- unlist(strsplit(j,"[.]",perl = T))[2]
        pw.data[i,j] <- subset(htc.be007, code==tmp.code
                               & Sample==tmp.sample)[[tmp.subject]]
    }
}

### 3. Average replications
pw.data$BE.Area.mean <- rowMeans(pw.data[,grep("BE.*Area",colnames(pw.data),perl = T)])
pw.data$BE.Count.mean <- rowMeans(pw.data[,grep("BE.*Count",colnames(pw.data),perl = T)])
pw.data$CRD.Area.mean <- rowMeans(pw.data[,grep("CRD.*Area",colnames(pw.data),perl = T)])
pw.data$CRD.Count.mean <- rowMeans(pw.data[,grep("CRD.*Count",colnames(pw.data),perl = T)])

### 4. Normalization
# use all the values without quality control.

