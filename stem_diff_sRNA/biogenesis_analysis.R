

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


options(scipen=999)




#import the desired read counts
read_counts <- read.delim("../sRNA_counts/small_sa_og_norRNA_25_33_counts_table.txt", sep=",", stringsAsFactors=FALSE, row.names = 1)[,]
read_counts <- read_counts[,4:19]
read_counts <- read_counts[,c(9:16, 1:8)]
read_counts <- read_counts[,-c(3,11)]



#filter reads
rRNA <- read.table("../data/small_sa_og_norRNA_25_33_rRNA_list.txt", header = FALSE, stringsAsFactors = FALSE)
read_counts <- read_counts[!(rownames(read_counts) %in% rRNA$V1) ,]
#piRNAdb <- read.table("../data/small_sa_og_norRNA_25_33_piRNAdb_list.txt", header = FALSE, stringsAsFactors = FALSE)
ncRNA <- read.table("../data/small_sa_og_norRNA_25_33_ncRNA_list.txt", header = FALSE, stringsAsFactors = FALSE)
#read_counts <- read_counts[rownames(read_counts) %in% piRNAdb$V1 | !(rownames(read_counts) %in% ncRNA$V1),]
read_counts <- read_counts[!(rownames(read_counts) %in% ncRNA$V1),]


#unique mapping reads
unique_mapping <- read.delim("../data/small_sa_og_norRNA_25_33_list_unique_mapping_genome_v0_m1.txt", header=FALSE, stringsAsFactors=FALSE)
#read_counts <- read_counts[rownames(read_counts) %in% unique_mapping$V1 | (rownames(read_counts) %in% piRNAdb$V1),]
#read_counts <- read_counts[rownames(read_counts) %in% unique_mapping$V1 ,]







#loading in the bowtie alignments from piRNA against cDNA
piRNA_gene_alignments <- read.delim("../alignments/biogenesis_alignments/sa_og_norRNA_25_33_noncRNA_cDNA_v0.txt", header=FALSE, stringsAsFactors=FALSE)[,1:3]
colnames(piRNA_gene_alignments) <- c("piRNA", "gene", "bp_start")

#piRNA_gene_alignments <- piRNA_gene_alignments[piRNA_gene_alignments$piRNA %in% rownames(read_counts),]
piRNA_gene_alignments <- unique(piRNA_gene_alignments)

piRNA_gene_alignments$hgnc_symbol <- sapply(piRNA_gene_alignments$gene, FUN = function(x){unlist(strsplit(x, split = " "))[2]})
piRNA_gene_alignments$length <- sapply(piRNA_gene_alignments$gene, FUN = function(x){unlist(strsplit(x, split = " "))[3]})
piRNA_gene_alignments$length <-  as.numeric(sapply(piRNA_gene_alignments$length, sub, pattern =  "transcript_length:", replacement =  ""))


#find the distance alignment occured from 3 prime end of transcript
piRNA_gene_alignments$distance_end <- piRNA_gene_alignments$length - piRNA_gene_alignments$bp_start



test <- read_counts


test$piRNA <- rownames(test)

test <- merge(test, piRNA_gene_alignments[,c(1,4)], by = 'piRNA', all = TRUE)

test$hgnc_symbol <- sapply(FUN = function(row){
  if(is.na(test$hgnc_symbol[row])){
    return(test$piRNA[row])
  }
  else{
    return(test$hgnc_symbol[row])
  }
}, X = 1:nrow(test))

test <- subset(test, select = -piRNA)

test <- aggregate(data = test, .~hgnc_symbol, sum)

rownames(test) <- test$hgnc_symbol

test <- subset(test, select = -hgnc_symbol)

test <- test[rowSums(test > 1) >= 6,]


read_counts <- test
#keep alignments from DE sRNA
#piRNA_gene_alignments <- piRNA_gene_alignments[piRNA_gene_alignments$piRNA %in% voom_top_sRNA$sRNA,]






#load in the differentially expressed gene

#voom analysis
#cell lines of cancer and differentiated
design.array <- c(rep("A",7), rep("B",7))
cellline <- factor(c(1:7, 1:7), levels = 1:7)
celltype <- factor(design.array)
design <-model.matrix(~cellline + celltype)

#converting to DGEList and keeping only expressions significantly above 0
dge_counts <- DGEList(counts=read_counts)
#dge_counts <- dge_counts[rowSums(cpm(dge_counts) >=1) >=7,]
dge_counts <- calcNormFactors(dge_counts, method = c("TMM"))

#voom normalization
voom_counts <- voomWithQualityWeights(dge_counts,design, plot = TRUE, normalize.method = "quantile")

voom_counts <- cpm(dge_counts,log=TRUE, prior.count = 3)
#MDSplot
plotMDS(dge_counts, labels = colnames(read_counts), cex = 1, main = "MDS Plot for piRNA Expressions", dim.plot = c(1,2))

#create the fit for the counts
fit <- eBayes(lmFit(voom_counts,design))


#voom fold changesmiRNA_alignments
voom_sRNA <- topTable(fit, coef = "celltypeB", sort.by="logFC", number =nrow(voom_counts) )
voom_top <- topTable(fit,coef = "celltypeB", sort.by="logFC", number =nrow(voom_counts),lfc = 2,p.value = 0.05)
#voomNormExp <- voomCounts$E
rm(celltype);rm(design.array); rm(design); rm(dge_counts)
#my name function

voom_sRNA$sRNA <- rownames(voom_sRNA)


voom_top_sRNA <- voom_top
voom_top_sRNA$sRNA <- rownames(voom_top_sRNA)

#derp <- voom_counts

#voom_counts <- list()
voom_counts$E <- as.data.frame(voom_counts)
#voom_counts$E <- read_counts




