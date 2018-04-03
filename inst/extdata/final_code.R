# Data processing
library(tibble)
library("phyloseq")

# Loading data
sal_A <- read.table("david14_otus_metadata_saliva_A.txt", header = TRUE, row.names = 1)
stool_B <- read.table("david14_otus_metadata_stool_B.txt",header = TRUE, row.names = 1)
stool_A <- read.table("david14_otus_metadata_stool_A.txt", header = TRUE, row.names = 1)

# Merging data
sal_A<-t(sal_A)
sal_A <-data.frame(sal_A)
sal_A<-add_column(sal_A, sample_type = "Saliva", .after = "collection.day")
stool_A<-t(stool_A)
stool_A <-data.frame(stool_A)
stool_A<-add_column(stool_A, sample_type = "Fecal", .after = "collection.day")
stool_B<-t(stool_B)
stool_B <-data.frame(stool_B)
stool_B<-add_column(stool_B, sample_type = "Fecal", .after = "collection.day")
sample_whole<-rbind(sal_A,stool_A,stool_B)
rm(stool_A)
rm(stool_B)
rm(stool.A)
rm(sal_A)


# Meta Data
meta.W<-sample_whole[,5432:5444 ]
colnames(meta.W)[1] <- "time"
meta.W<-add_column(meta.W, gender = "male", .after = "time")
meta.W<-add_column(meta.W, subject = ifelse(meta.W$AGE == 26, "A", "B"), .after = "gender")
meta.W<-add_column(meta.W, sample = rownames(meta.W), .after = "sample_type")
meta.W<-add_column(meta.W[-15],age = as.integer(26) , .after = "sample")



#filter data (as in David et al. 2014)

toDrop <- c("SAMPLE_Stool76", "SAMPLE_Stool77","SAMPLE_Stool339", "SAMPLE_Stool340", "SAMPLE_Stool341", "SAMPLE_Stool342", "SAMPLE_Stool343", "SAMPLE_Stool344", "SAMPLE_Stool345", "SAMPLE_Stool346", "SAMPLE_Stool349", "SAMPLE_Stool348", "SAMPLE_Stool347","SAMPLE_Stool540","SAMPLE_Stool541")
meta.W<- meta.W[!(rownames(meta.W) %in% toDrop), ]
sample_whole<-sample_whole[!(rownames(sample_whole) %in% toDrop), ]

sample.whole <-sample_whole[,1:5431]
rm(sample_whole)
sample.whole <- as.matrix(sample.whole)

## Phyloseq object

# 1) otu_table

seq.data<- sample.whole+1
seq.data<- as.matrix(t(sample.whole))


#2) sample_data
metadata.W <- meta.W


#change otu names to character dataframe
coln.A <- names(metadata.W[,c(2,3,6:15)])
metadata.W[, coln.A] <- sapply(metadata.W[, coln.A], as.character)


#3) #otu
otu <- read.table("otu_table-metadata.txt", sep = "\t", fill= TRUE)
colnames(otu)[1] <- "OTU"
colnames(otu)[2] <- "Domain"
colnames(otu)[3] <- "Phylum"
colnames(otu)[4] <- "Class"
colnames(otu)[5] <- "Order"
colnames(otu)[6] <- "Family"
colnames(otu)[7] <- "Genus"
colnames(otu)[8] <- "Species"
rownames(otu) = otu[,"OTU"]
otu <- otu[,c(1:9)]



seqa_otu = otu_table(seq.data, taxa_are_rows = TRUE)
metaa_sample = sample_data(metadata.W)
taxaa_tax = tax_table(as.matrix(otu))
david_2014<- phyloseq(seqa_otu, metaa_sample, taxaa_tax)

save(david_2014, file = "david_2014.rda")
rm(david_2014,metaa_sample, meta.W,coln.A,metadata.W,otu,sample.whole,seq.data,seqa_otu,taxaa_tax,toDrop)
