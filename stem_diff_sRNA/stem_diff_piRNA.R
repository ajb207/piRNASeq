

#library(edgeR)
#library(biomaRt)
require(ggplot2)
#library(gplots)
require(knitr)
require(plyr)
require(Biostrings)
#require(ReactomePA)

EnsembltoHGNC <- function(countsData){
  
  
  require(biomaRt)
  
  ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org") # v82
  listMarts(mart=ensembl)
  # list all the Datasets in Ensembl biomart I connected
  listDatasets(ensembl) # 69 of them
  
  # I can also save Datasets available to me to a dataframe
  myDatasets<-listDatasets(ensembl)
  #View(myDatasets)
  # Ensembl currently contains >50 datasets~species
  
  # create a handle to connect Homo sapiens gene data
  ###n/a### ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                    host = "www.ensembl.org", dataset="hsapiens_gene_ensembl")
  
  # find gene symbols for affy ids from chip HG-U133+v2
  #myAttributes<-listAttributes(ensembl)
  
  # or you can just use this to get all rows that have U133
  #myAttributes[grep("Ensembl",myAttributes$description),]
  
  
  # To build a query for biomart, you need 3 things:
  # provide attriutes (col names), filters (to filter rows), and values to query.
  annotation = getBM(attributes=c("ensembl_gene_id" , "hgnc_symbol" ), 
                     filters="ensembl_gene_id",  values=row.names(countsData),
                     mart = ensembl)
  
  countsData$ensembl_gene_id <- rownames(countsData)
  countsData <- merge(countsData, annotation, by = "ensembl_gene_id")
  #countsData <- countsData[countsData$hgnc_symbol != "",]
  #rownames(countsData) <- countsData$hgnc_symbol
  #countsData <- countsData[-c(1,22)]
  return(countsData)
  
}

#load in the piRNA -> mRNA alignments
alignments <- read.csv("../alignments/stem_diff/transcript_genename_alignments/sa_og_25_33_no_sRNA_list_3UTR_transcriptID_alignments.txt", sep="\t", stringsAsFactors=FALSE, header = FALSE)
alignments <- alignments[,c(1,3,5)]
colnames(alignments) <- c("piRNA", "Target", "Alignment")
alignments <- unique(alignments)
#alignments$Gene <- alignments$Target



test <- as.data.frame(do.call(rbind, strsplit(alignments$Target[1:10000000], "\\|")))
test2 <- as.data.frame(do.call(rbind, strsplit(alignments$Target[10000001:20000000], "\\|")))
test3 <- as.data.frame(do.call(rbind, strsplit(alignments$Target[20000001:33234619], "\\|")))

test.all <- rbind(test,test2,test3)
test.all$V1 <- as.character(test.all$V1)
test.all$V2 <- as.character(test.all$V2)
rm(test); rm(test2);rm(test3)
alignments$Transcript_ID <- test.all$V1
alignments$Target <- test.all$V2
rm(test.all)


#list of protein coding genes
#only protein coding genes are of interest
Homo_sapiens <- read.delim("../../Homo_sapiens.gene_info", header=FALSE, comment.char="#", stringsAsFactors = FALSE)
pc_genes <- Homo_sapiens[Homo_sapiens$V10 == "protein-coding",]
pc_genes <- unlist(pc_genes$V3)
rm(Homo_sapiens)



#only keep alignments to protein coding genes
alignments <- alignments[alignments$Target %in% pc_genes,]
rm(pc_genes)


alignments <-read.table("final_alignments.txt")

alignments <- alignments[,1:3]
alignments <- unique(alignments)


#load in the sequence counts for each cell line
all_counts <- read.delim("../counts/sa_og_25_33_no_sRNA_counts.txt", stringsAsFactors=FALSE, header = TRUE, sep = "")
rownames(all_counts) <- all_counts$piRNA


S7_counts <- subset(all_counts, select = -piRNA)
colnames(S7_counts) <- c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P", "W", "X", "Y", "Z")
piRNA_totals <- as.data.frame(colSums(S7_counts))
colnames(piRNA_totals) <- "Counts"
barplot(piRNA_totals$Counts, names.arg = rownames(piRNA_totals), main = "piRNA Total Counts")
S7_counts <- S7_counts[rowSums(S7_counts) > 0,]

norm_counts <- S7_counts
norm_counts$I <- norm_counts$I*(sum(norm_counts$A/norm_counts$I)) 
norm_counts$J <- norm_counts$J*(sum(norm_counts$B/norm_counts$J)) 
norm_counts$K <- norm_counts$K*(sum(norm_counts$C/norm_counts$K)) 
norm_counts$L <- norm_counts$L*(sum(norm_counts$D/norm_counts$L)) 
norm_counts$M <- norm_counts$M*(sum(norm_counts$E/norm_counts$M)) 
norm_counts$N <- norm_counts$N*(sum(norm_counts$F/norm_counts$N)) 
norm_counts$O <- norm_counts$O*(sum(norm_counts$G/norm_counts$O)) 
norm_counts$P <- norm_counts$P*(sum(norm_counts$H/norm_counts$P)) 

S7_counts <- norm_counts
S7_counts$piRNA <- rownames(S7_counts)


#calculate the number of targets a gene had
piRNA_number_of_targets <- as.data.frame(table(alignments$piRNA))
colnames(piRNA_number_of_targets) <- c("piRNA", "NumbTar")

gene_number_targeted <- as.data.frame(table(alignments$Target))
colnames(gene_number_targeted) <- c("Gene", "NumbTar")

S7_counts <- merge(S7_counts, piRNA_number_of_targets, by = "piRNA")


rownames(S7_counts) <- S7_counts$piRNA
S7_counts <- subset(S7_counts, select = -piRNA)


for(row in 1:nrow(S7_counts)){
  
  S7_counts[row,1:20] <- S7_counts[row,1:20] / S7_counts[row,21]
}

S7_counts$piRNA <- rownames(S7_counts)
S7_counts<- merge(S7_counts, alignments, by = "piRNA")


total_targeted <- ddply(S7_counts, "Target", numcolwise(sum))
rm(alignments)
rm(all_counts)





total_targeted[,2:21] <- round(total_targeted[,2:21],5)

colnames(total_targeted) <- c("Target", "A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P", "W", "X", "Y", "Z", "NumTar")
total_targeted[,2:21] <- round(total_targeted[,2:21],5)
total_targeted[,2:21] <- signif(total_targeted[,2:21],5)

write.table(total_targeted, "total_targeted.txt", quote = FALSE)

total_targeted <- read.table("total_targeted.txt")



AI <- as.data.frame(total_targeted[,c(1,2,10)])
AI$difference <- AI$A-AI$I

BJ <- as.data.frame(total_targeted[,c(1,3,11)])
BJ$difference <- BJ$B-BJ$J

CK <- as.data.frame(total_targeted[,c(1,4,12)])
CK$difference <- CK$C-CK$K

DL <- as.data.frame(total_targeted[,c(1,5,13)])
DL$difference <- DL$D-DL$L

EM <- as.data.frame(total_targeted[,c(1,6,14)])
EM$difference <- EM$E-EM$M

FN <- as.data.frame(total_targeted[,c(1,7,15)])
FN$difference <- FN$F-FN$N

GO <- as.data.frame(total_targeted[,c(1,8,16)])
GO$difference <- GO$G-GO$O

HP <- as.data.frame(total_targeted[,c(1,9,17)])
HP$difference <- HP$H-HP$P


