



#loading required packages
require(knitr)
require(limma)
require(edgeR)
require(biomaRt)
require(ggplot2)
require(NMF)
require(DESeq2)
require(statmod)


library(fgsea)
library(data.table)
library(ggplot2)
require("reactome.db")



#OLD NO LONGER IN USE function to annotate ensembl gene IDs to HGNC gene symbols
annotate_ENSG_TopTable<-function(tops){
  # library(biomaRt)
  #tops<-voomCounts$E
  
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
                     filters="ensembl_gene_id",  values=row.names(tops),
                     mart = ensembl)
  #head(annotation)
  #   if(length(unique(annotation$ensembl_gene_id)) < length(annotation$ensembl_gene_id)){
  #     annotation<-annotation[- grep("ENSG00000207704",annotation$ensembl_gene_id)[2],]  
  #   }
  row.names(annotation)<-make.names(annotation$ensembl_gene_id,unique = TRUE)
  
  #save.image(file="top.genes.Rdata")
  #load("top.genes.Rdata")
  m<-match(row.names(tops),annotation$ensembl_gene_id)
  #topGenes<-merge(annotation,tops,by.x="row.names",by.y="row.names")
  topGenes<-cbind(annotation,tops[m,])
  #topGenes$hgnc_symbol
  return(topGenes);
}


#function to annotate ensembl gene IDs to HGNC gene symbols
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
  countsData <- merge(countsData, annotation, by = "ensembl_gene_id", all.x = T)
  #countsData <- countsData[countsData$hgnc_symbol != "",]
  #rownames(countsData) <- countsData$hgnc_symbol
  #countsData <- countsData[-c(1,22)]
  return(countsData)
  
}


#Function to import star counts from separate output folders
#requires the folder containing all the star counts and the column for use. Almost always column 2
importStarCounts <- function(folder, column){
  
  temp <- list.files(folder, full.names = TRUE)
  counts <- lapply(temp, read.delim, header=FALSE, stringsAsFactors = FALSE)
  file_names <- list.files(folder)
  
  counts <- lapply(counts, `[`, c(1,column))
  #turn each express output into a dataframe and merge them
  final_frame <- as.data.frame(counts[1])
  colnames(final_frame) <- c("ensembl_id", file_names[1])
  for(n in 2:length(counts)){
    temp_frame <- as.data.frame(counts[n])
    colnames(temp_frame) <- c("ensembl_id", file_names[n])
    final_frame <- merge(final_frame, temp_frame, by = "ensembl_id", all =TRUE)
    
  }
  rownames(final_frame) <- final_frame$ensembl_id
  final_frame <- subset(final_frame, select = -ensembl_id)
  #final_frame[is.na(finalFrame)] <- 0
  
  #final_frame <- final_frame[!final_frame$gene=="N_unmapped",]
  #final_frame <- final_frame[!final_frame$gene=="N_multimapping",]
  #final_frame <- final_frame[!final_frame$gene=="N_noFeature",]
  #final_frame <- final_frame[!final_frame$gene=="N_ambiguous",]
  
  
  return(final_frame)
}




#load in the star counts
read_counts <- importStarCounts("../SRP040525_gene_counts/", column = 4)
colnames(read_counts) <- gsub(pattern = "ReadsPerGene.out.tab", "", colnames(read_counts))
#write.table(test, "../data/read_counts.txt", quote = FALSE, row.names = F)


#remove reads that did not get counted due to multi mapping and unalignment
read_counts <- read_counts[!rownames(read_counts)=="N_unmapped",]
read_counts <- read_counts[!rownames(read_counts)=="N_multimapping",]
read_counts <- read_counts[!rownames(read_counts)=="N_noFeature",]
read_counts <- read_counts[!rownames(read_counts)=="N_ambiguous",]


#select only the relevant sample sets
read_counts <- read_counts[,c(2,5,8,11)]

#factor of the cell types, knockdown or wild type
kd <- factor(c(rep("WT",2),   #D stands for differentiated cell or non-GSCs
               rep("KD",2)    #S stands for stem cell of GSC
               #rep("C",2),
               # rep("D",2)
), levels = c("WT","KD")
)
celltype <- factor(c("S","D", "S","D"), levels = c("D", "S"))  #Factor for the samples used. 8 paired samples
#create the design matrix based on the factors above
design <- model.matrix(~  kd + celltype)


#converting to DGEList
dge_counts <- DGEList(counts=read_counts[rowSums(cpm(read_counts) >= 1) >= 4,])

#normalization
dge_counts <- calcNormFactors(dge_counts, method = c("TMM"))

#voom for creating converting RNA seq reads to microaray data for the pipeline
voom_counts <-voom(dge_counts, design, plot=T)

#voom_counts <- cpm(dge_counts,log=T, prior.count = 3)

#MDS plot of the data
plotMDS(dge_counts, cex = 1, main = "MDS of Stem vs Diff", labels = colnames(read_counts), top = 500, dim.plot = c(1,2), gene.selection = 'common')


#linear fit and eBayes funtions for differential expression
fit <-  eBayes(lmFit(voom_counts,design))

#voom fold changes and p values for all genes
voom_all <- topTable(fit, coef = "celltypeS", sort.by="p", number =nrow(voom_counts))
#voom fold changes and p values based on significant hits
voom_top <- topTable(fit,coef = "celltypeS", sort.by="p", number =nrow(voom_counts),lfc = 2,p.value = 0.05)


#name conversion to HGNC symbols
voom_all_HGNC <- EnsembltoHGNC(voom_all)
voom_top_HGNC <- EnsembltoHGNC(voom_top)



#write the data
write.table(x = voom_all_HGNC, "../data/StemVsDiff_all_DE_genes.txt", quote = FALSE, sep = ",", row.names = F)
write.table(x = voom_top_HGNC, "../data/StemVsDiff_top_DE_genes.txt", quote = FALSE, sep = ",", row.names = F)



