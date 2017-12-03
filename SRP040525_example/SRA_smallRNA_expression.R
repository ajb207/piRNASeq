


require(limma)
require(knitr)
require(edgeR)
require(biomaRt)
require(DESeq2)
require(edgeR)
require(ggplot2)
require(reshape)
require(reshape2)
require(stringdist)
require(apcluster)
require(qqman)
require(Biostrings)
require(motifRG)
require(stringr)
require(scales)
require(NMF)


#force R to not use scientific notation
options(scipen=999)

#compute various information about the alignments.
alignmentInfo <- function(alignments){
  
  
  #find the distance alignment occured from 3 prime end of transcript
  #alignments$bp_end <- alignments$length - alignments$bp_start
  #count the number of mismatches
  #alignments$number_mismatches <- str_count(alignments$mismatches, ">")
  
  #finding total targets for each piRNA
  target_number <- as.data.frame(table(unique(subset(alignments, select = c(piRNA,hgnc_symbol)))[,1]))
  colnames(target_number) <- c("piRNA", "targets")
  alignments <- merge(alignments, target_number, by = "piRNA")
  
  #find total time gene is targeted by unique piRNA
  target_number <- as.data.frame(table(as.data.frame(unique(subset(alignments, select = c(piRNA,hgnc_symbol)))[,2])))
  colnames(target_number) <- c("hgnc_symbol", "gene_frequency")
  alignments <- merge(alignments, target_number, by = "hgnc_symbol")
  
  #denote weather the piRNA and/or the gene is differentially expressed
  alignments$DE_piRNA <- FALSE
  alignments$DE_piRNA[alignments$piRNA %in% rownames(voom_top)] <- TRUE
  
  alignments$DE_gene <- FALSE
  alignments$DE_gene[alignments$hgnc_symbol %in% top_genes$hgnc_symbol] <- TRUE
  
  #denote whether the 1T 10A bias is seen
  #alignments$`1T` <-  ifelse(sapply(alignments$piRNA, substr,1,1) == "T",TRUE,FALSE) 
  #alignments$`10A` <-  ifelse(sapply(alignments$piRNA, substr,10,10) == "A",TRUE,FALSE)
  #alignments$`1T10A` <- ifelse((alignments$`1T` & alignments$`10A`), TRUE, FALSE)
  
  
  target_number <- as.data.frame(table(unique(subset(alignments, select = c(piRNA,transcript_id)))[,2]))
  colnames(target_number) <- c("transcript", "targets")
  #rpkm <- merge(subset(alignments, select = c(transcript, length)) , target_number, by = "transcript")
  #rpkm$rpkm <- rpkm$targets/((rpkm$length)/1000)
  
  #alignments <- merge(alignments, subset(rpkm, select = c(transcript, rpkm)), by = "transcript")
  
  return(alignments)
}





#import the desired read counts

read_counts <- read.delim("../sRNA_counts/SRP040525_counts.txt", sep=",", stringsAsFactors=FALSE, row.names = 1)[,]
#read_counts <- read_counts[,4:19]
colnames(read_counts) <- gsub(colnames(read_counts), pattern ="_sa_collapsed.fq.", replacement =  "")
og <- read.table("../alignments/SRP040525_sa_alignments_genome_list.txt", header = FALSE, stringsAsFactors = FALSE)
rRNA <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_rRNA_list.txt", header = FALSE, stringsAsFactors = FALSE)
YRNA <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_YRNA_list.txt", header = FALSE, stringsAsFactors = FALSE)
miRNA <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_miRNA_list.txt", header = FALSE, stringsAsFactors = FALSE)
read_counts <- read_counts[rownames(read_counts) %in% og$V1,]
read_counts <- read_counts[!rownames(read_counts) %in% rRNA$V1,]
read_counts <- read_counts[!rownames(read_counts) %in% YRNA$V1,]
read_counts <- read_counts[!rownames(read_counts) %in% miRNA$V1,]
read_counts <- read_counts[nchar(rownames(read_counts)) >= 25 & nchar(rownames(read_counts)) <= 33,]




read_counts <- read_counts[,c(5,11,2,8)]


#factor for the cell tyoe
celltype <- factor(c(rep("D",2),  #D is differentiated of non-GSC
                     rep("S",2)), #S is stem or GSC
                   levels = c("D", "S"))

#factor for the paired samples
kd <- factor(c("WT", "KD", "WT", "KD"), levels = c("WT","KD"))
#design matix of the experiment
#cellline <- factor(c(1,1:5,1,1:5), levels = c(2:5,1))

design <-model.matrix(~ kd + celltype )

#converting to DGEList and keeping only expressions significantly above 0
dge_counts <- DGEList(counts=read_counts[rowSums(cpm(read_counts) >= 1) >= 2,])
#dge_counts <- dge_counts[rowSums(cpm(dge_counts) >=1) >=7,]
dge_counts <- calcNormFactors(dge_counts, method = c("TMM"))

#voom normalization
voom_counts <- cpm(dge_counts,log=TRUE, prior.count = 3)
voom_counts <- voom(dge_counts, design = design)
#MDSplot
plotMDS(dge_counts,  cex = 1, main = "MDS Plot for piRNA Expressions", dim.plot = c(1,2))#, labels = c("WTD3", "WTD4", "WTD5", "WTD6", "KDD3", "KDD4", "KDD5", "KDD6", "WTS1", "WTS2", "KDS1", "KDS2" ))

#create the linear fit for the counts
fit <- lmFit(voom_counts,design)
fit <- eBayes(fit, trend = TRUE)
fit <- eBayes(fit, trend = F)


#voom fold changesmiRNA_alignments
voom_sRNA <- topTable(fit, coef = "celltypeS", sort.by="logFC", number =nrow(voom_counts) )
voom_top <- topTable(fit,coef = "celltypeS", sort.by="logFC", number =nrow(voom_counts),lfc = 2,p.value = 0.05)
#voomNormExp <- voomCounts$E
#rm(celltype);rm(design.array); rm(design); rm(dge_counts)
#my name function

#adding the reads names into its own column for merging
voom_sRNA$sRNA <- rownames(voom_sRNA)
voom_top_sRNA <- voom_top
voom_top_sRNA$sRNA <- rownames(voom_top_sRNA)







#tRNA <- read.table("../alignments/biogenesis_small_sa_og_norRNA/small_sa_og_norRNA_alignments_tRNA.txt", header = FALSE, stringsAsFactors = FALSE)
lncRNA <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_lncRNA_list.txt", header = FALSE, stringsAsFactors = FALSE)
unique <- read.table("../alignments/SRP040525_sa_alignments_uniquemapping_list.txt", header = FALSE, stringsAsFactors = FALSE)
piRBase <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_piRNA_list.txt", header = FALSE, stringsAsFactors = FALSE)
piRBaseplus15 <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_piRBaseplus15_list.txt", header = FALSE, stringsAsFactors = FALSE)
YRNA <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_YRNA_list.txt", header = FALSE, stringsAsFactors = FALSE)
snoRNA <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_snoRNA_list.txt", header = FALSE, stringsAsFactors = FALSE)
miRNA <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_miRNA_list.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
mRNA <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_cDNA_list.txt", header = FALSE, stringsAsFactors = FALSE)
rmsk <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_repeatmasker.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
colnames(rmsk) <- c("sRNA", "Repeat Name")
rmsk$`Repeat Name` <- gsub(pattern = "_", replacement = "", rmsk$`Repeat Name`)
rmsk_plus40 <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_repeatmaskerplus40.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(rmsk_plus40) <- c("sRNA", "Extended Repeat Name")
rmsk_plus40$`Extended Repeat Name` <- gsub(pattern = "_", replacement = "", rmsk_plus40$`Extended Repeat Name`)
tRNA <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_repeatmasker_tRNA.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(tRNA) <- c("sRNA", "tRNA Name")
tRNA_plus40 <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_repeatmaskerplus40_tRNA.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(tRNA_plus40) <- c("sRNA", "Extended tRNA Name")



voom_top_sRNA$'Targeting Sequence' <- strtrim(voom_top_sRNA$sRNA,25)
voom_top_sRNA$'Targeting Sequence Minus 1st 5`' <- substring(voom_top_sRNA$'Targeting Sequence',2)
voom_top_sRNA$piRBase <- ifelse(voom_top_sRNA$sRNA %in% piRBase$V1, T, F)
voom_top_sRNA$'Extended piRBase' <- ifelse(voom_top_sRNA$sRNA %in% piRBaseplus15$V1, T, F)
voom_top_sRNA$'Unique Mapping' <- ifelse(voom_top_sRNA$sRNA %in% unique$V1, T, F)
voom_top_sRNA$miRNA <- ifelse(voom_top_sRNA$sRNA %in% miRNA$V1, T, F)
voom_top_sRNA$mRNA <- ifelse(voom_top_sRNA$sRNA %in% mRNA$V1, T, F)
voom_top_sRNA$lncRNA <- ifelse(voom_top_sRNA$sRNA %in% lncRNA$V1, T, F)
voom_top_sRNA$snoRNA <- ifelse(voom_top_sRNA$sRNA %in% snoRNA$V1, T, F)
voom_top_sRNA$YRNA <- ifelse(voom_top_sRNA$sRNA %in% YRNA$V1, T, F)
voom_top_sRNA$tRNA <- ifelse(voom_top_sRNA$sRNA %in% tRNA$V1, T, F)
voom_top_sRNA <- merge(voom_top_sRNA, tRNA[,c(1:2)], by = 'sRNA', all.x = T)
voom_top_sRNA$'Extended tRNA' <- ifelse(voom_top_sRNA$sRNA %in% tRNA_plus40$'sRNA', T, F)
voom_top_sRNA <- merge(voom_top_sRNA, tRNA_plus40[,c(1:2)], by = 'sRNA', all.x = T)
voom_top_sRNA$'Repeat Masker' <- ifelse(voom_top_sRNA$sRNA %in% rmsk$sRNA, T, F)
voom_top_sRNA <- merge(voom_top_sRNA, rmsk[,c(1:2)], by = 'sRNA', all.x = T)
voom_top_sRNA$'Extended Repeat Masker' <- ifelse(voom_top_sRNA$sRNA %in% rmsk_plus40$V1, T, F)




write.table(voom_top_sRNA, "../data/SRP040525_DE_putative_piRNA_25_33.csv", sep = ",", row.names = F, quote = F)

#annotating all reads



voom_sRNA$'Targeting Sequence' <- strtrim(voom_sRNA$sRNA,25)
voom_sRNA$'Targeting Sequence Minus 1st 5`' <- substring(voom_sRNA$'Targeting Sequence',2)
voom_sRNA$piRBase <- ifelse(voom_sRNA$sRNA %in% piRBase$V1, T, F)
voom_sRNA$'Extended piRBase' <- ifelse(voom_sRNA$sRNA %in% piRBaseplus15$V1, T, F)
voom_sRNA$'Unique Mapping' <- ifelse(voom_sRNA$sRNA %in% unique$V1, T, F)
voom_sRNA$miRNA <- ifelse(voom_sRNA$sRNA %in% miRNA$V1, T, F)
voom_sRNA$mRNA <- ifelse(voom_sRNA$sRNA %in% mRNA$V1, T, F)
voom_sRNA$lncRNA <- ifelse(voom_sRNA$sRNA %in% lncRNA$V1, T, F)
voom_sRNA$snoRNA <- ifelse(voom_sRNA$sRNA %in% snoRNA$V1, T, F)
voom_sRNA$YRNA <- ifelse(voom_sRNA$sRNA %in% YRNA$V1, T, F)
voom_sRNA$tRNA <- ifelse(voom_sRNA$sRNA %in% tRNA$V1, T, F)
voom_sRNA <- merge(voom_sRNA, tRNA[,c(1:2)], by = 'sRNA', all.x = T)
voom_sRNA$'Extended tRNA' <- ifelse(voom_sRNA$sRNA %in% tRNA_plus40$'sRNA', T, F)
voom_sRNA <- merge(voom_sRNA, tRNA_plus40[,c(1:2)], by = 'sRNA', all.x = T)
voom_sRNA$'Repeat Masker' <- ifelse(voom_sRNA$sRNA %in% rmsk$sRNA, T, F)
voom_sRNA <- merge(voom_sRNA, rmsk[,c(1:2)], by = 'sRNA', all.x = T)
voom_sRNA$'Extended Repeat Masker' <- ifelse(voom_sRNA$sRNA %in% rmsk_plus40$V1, T, F)


write.table(voom_sRNA, "../data/SRP040525_putative_piRNA_25_33.csv", sep = ",", row.names = F, quote = F)
