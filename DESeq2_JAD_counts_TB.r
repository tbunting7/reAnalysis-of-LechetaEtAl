BiocManager::install("tximportData", force = T)
BiocManager::install("tximport", force = T)
BiocManager::install("DESeq2", force = T)
BiocManager::install('PCAtools')
install.packages('emmeans')
install.packages('sjPlot')
install.packages('VennDiagram')
install.packages("gghighlight")

library(tximportData)
library(tximport)
library(DESeq2)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(car)
library(ggpubr)
library(PCAtools)
library(emmeans)
library(sjPlot)
library(VennDiagram)
library(gghighlight)
library(stringr)
library(edgeR)

getwd()
wd <- "/media/raglandlab/ExtraDrive1/infoGen24/tristanb/starMapped/"
setwd(wd)

# create a file list containing all the files with gene results
file_list <- list.files(path = wd, pattern = '*e.out.tab', recursive = T)
file_list


# read in your sample information document
sampleInfo<-read.csv('/media/raglandlab/ExtraDrive1/infoGen24/tristanb/sampleInfoTemperature.csv',stringsAsFactors = F)
sampleInfo$temperature <- factor(sampleInfo$temperature, levels = c("25","4","37"))

# empty table to populate with sample ID's with the first parts of the filenames
table_all<-c()

## Loops through file list and pull out ID's and pulls out expected counts for each sample ID and for each gene
# Example
# gene id     sample 1      sample 2
# gene 1         5             2
# gene 2        10            400

## Use when analysis method was RNAStar
#file=1
for (file in file_list){
  sampleID<-str_match(file,"([A-Za-z]*[0-9])")[,1]
  if (is.na(sampleID)){
    sampleID<-str_match(file,"([0-9])*_D[A-Za-z]*__[A-Za-z]*_[A-Za-z]*[0-9]")[,1]
  }
  print(sampleID)
  df<-read.table(paste0(wd,file),header=F,sep="\t",row.names=NULL)
  df2<-df[c(-1:-4), #remove the first 4 rows
          c(1,3) #keep only columns 1 (geneID) and 3 (expected counts)
  ]
  colnames(df2)[1] <- paste("gene")
  colnames(df2)[2] <- paste(sampleID,"_exp_count", sep=""
  )
  if (length(table_all) ==0){table_all<-df2
  } else {
    table_all<-merge(table_all,df2,by=c(1)) }
}

# see how many genes are in the count table
length(table_all$gene)
write.table(table_all, file = "/media/raglandlab/ExtraDrive1/infoGen24/tristanb/STAR_table.txt", sep = "\t")

## This function removes genes that do not have a least 1 count in at least 50% of the samples
filterMinCount<- function(x) {
  pres<-x >=1
  out=F
  if ((sum(pres)/length(pres)) >= 0.5) {out=T}
  return(out)
}


filterInd<-apply(table_all[,(-c(1))],1,filterMinCount) 
table_all<-table_all[filterInd,]
nrow(table_all)
# 9880 genes 
countsMatrix <- table_all[,(-c(1))]

# create a data frame to use in the DESeq2 models
InfoAll <- sampleInfo

## create a DESeq2 model
dds <- DESeqDataSetFromMatrix(countData = round(countsMatrix), # count data needs to be integers, i.e. no decimals. So round the matrix
                              colData = InfoAll, 
                              design = ~ temperature)
dds <- DESeq(dds)

# see if there are any samples that cluster together in a strange way or represent outliers
plotMDS(dds)

# compare contrasting conditions for results
res <- results(dds, contrast = c("temperature", "37", "25"), alpha = 0.05)
res <- res[order(res$padj),]
res.df <- as.data.frame(res)
res.df$gene <- rownames(res.df)
res.df.sig <- res.df %>% 
  filter(padj < 0.05)
length(res.df.sig$gene)
