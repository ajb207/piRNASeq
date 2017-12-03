
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


removeRepeats <- function(read_counts){
  
  YRNA <- read.table("../alignments/biogenesis_small_sa_og_norRNA/small_sa_og_norRNA_alignments_YRNA.txt", header = FALSE, stringsAsFactors = FALSE)
  snoRNA <- read.table("../alignments/biogenesis_alignments_norRNA_25_33/small_sa_og_norRNA_alignments_snoRNA.txt", header = FALSE, stringsAsFactors = FALSE)
  
  rmsk_plus40 <- read.table("../alignments/small_sa_og_norRNA_alignments_repeatmasker_plus40.txt", header = FALSE, stringsAsFactors = FALSE)
  rmsk_plus40_tRNA <- rmsk_plus40[grepl(pattern = "tRNA", rmsk_plus40$V2),]
  rmsk <- read.table("../alignments/biogenesis_small_sa_og_norRNA/small_sa_og_norRNA_alignments_repeatmasker.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  
  rmsk$Family <- ifelse(grepl("tRNA", rmsk$V2), "tRNA", 
                        ifelse(grepl(")n", rmsk$V2), "Simple", 
                               ifelse(grepl("Alu", rmsk$V2), "Alu", 
                                      ifelse(rmsk$V2 %in% c("U1", "U2", "U3", "U4", "U5", "U6", "U7", "U8"), "U", 
                                             ifelse(grepl("HY", rmsk$V2), "Y-RNA", 
                                                    ifelse(grepl(pattern =  paste(c("L1", "L2", "L3", "L4", "L5", "L6"), collapse = "|"), rmsk$V2), "LINE", 
                                                           ifelse(grepl("MIR", rmsk$V2), "MIR", "other" )))))))
  
  rmsk$Family[rmsk$V2 == "5S"] <- "5S"
  rmsk$Family[rmsk$V2 == "7SLRNA"] <- "7SLRNA"
  rmsk$Family[rmsk$V2 == "7SK"] <- "7SK"
  
  
  
  
  
  
  #read_counts <- read_counts[!(rownames(read_counts) %in% tRNA$V1),]
  #read_counts <- read_counts[!(rownames(read_counts) %in% lncRNA$V1),]
  read_counts <- read_counts[!(rownames(read_counts) %in% YRNA$V1),]
  #read_counts <- read_counts[!(rownames(read_counts) %in% snoRNA$V1),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "Simple"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "Alu"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "5S"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "7SLRNA"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "7SK"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "other"]),]
  
  return(read_counts)
}

