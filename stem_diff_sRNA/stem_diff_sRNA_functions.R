
#compute various information about the alignments.
alignmentInfo <- 
  function(alignments){
    
    
    #find the distance alignment occured from 3 prime end of transcript
    #alignments$bp_end <- alignments$length - alignments$bp_start
    #count the number of mismatches
    #alignments$number_mismatches <- str_count(alignments$mismatches, ">")
    
    #finding total targets for each piRNA
    target.number <- 
      alignments %>% 
      subset(select = c(piRNA,hgnc_symbol)) %>%
      table(unique()[,1]) %>%
      as.data.frame()
    colnames(target.number) <- c("piRNA", "targets")
    alignments <- merge(alignments,
                        target.number,
                        by = "piRNA")
    
    #find total time gene is targeted by unique piRNA
    target.number <- 
      alignments %>%
      subset(select = c(piRNA, hgnc_symbol)) %>%
      unique()[,2] %>%
      as.data.frame() %>%
      table() %>%
      as.data.frame()
    
    colnames(target.number) <- c("hgnc_symbol", "gene_frequency")
    alignments <- merge(alignments,
                        target.number, 
                        by = "hgnc_symbol")
    
    #denote weather the piRNA and/or the gene is differentially expressed
    alignments$DE_piRNA <- FALSE
    alignments$DE_piRNA[alignments$piRNA %in% rownames(voom_top)] <- TRUE
    
    alignments$DE_gene <- FALSE
    alignments$DE_gene[alignments$hgnc_symbol %in% top_genes$hgnc_symbol] <- TRUE
    
    #denote whether the 1T 10A bias is seen
    #alignments$`1T` <-  ifelse(sapply(alignments$piRNA, substr,1,1) == "T",TRUE,FALSE) 
    #alignments$`10A` <-  ifelse(sapply(alignments$piRNA, substr,10,10) == "A",TRUE,FALSE)
    #alignments$`1T10A` <- ifelse((alignments$`1T` & alignments$`10A`), TRUE, FALSE)
    
    
    target.number <- 
      alignments %>%
      subset(select = c(piRNA,transcript_id)) %>% 
      unique()[,2] %>%
      table() %>%
      as.data.frame()
    colnames(target.number) <- c("transcript", "targets")
    #rpkm <- merge(subset(alignments, select = c(transcript, length)) , target.number, by = "transcript")
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


removeRepeats <- function(read_counts){
  
  YRNA <- read.table("../alignments/biogenesis_small_sa_og_norRNA/small_sa_og_norRNA_alignments_YRNA.txt", header = FALSE, stringsAsFactors = FALSE)
  snoRNA <- read.table("../alignments/biogenesis_alignments_norRNA_25_33/small_sa_og_norRNA_alignments_snoRNA.txt", header = FALSE, stringsAsFactors = FALSE)
  
  rmsk_plus40 <- read.table("../alignments/small_sa_og_norRNA_alignments_repeatmasker_plus40.txt", header = FALSE, stringsAsFactors = FALSE)
  rmsk_plus40_tRNA <- rmsk_plus40[grepl(pattern = "tRNA", rmsk_plus40$V2),]
  #rmsk <- read.table("../alignments/biogenesis_small_sa_og_norRNA/small_sa_og_norRNA_alignments_repeatmasker.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  rmsk <- read.table("../alignments/small_sa_og_norRNA_alignments_repeatmasker_plus40.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  
  
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
  
  return(read_counts)
}

removeRepeats <- function(read_counts){
  
  tRNA <- read.table("../alignments/biogenesis_small_sa_og_norRNA/small_sa_og_norRNA_alignments_tRNA.txt", header = FALSE, stringsAsFactors = FALSE)
  lncRNA <- read.table("../alignments/biogenesis_small_sa_og_norRNA/small_sa_og_norRNA_alignments_lncRNA.txt", header = FALSE, stringsAsFactors = FALSE)
  piRNA <- read.table("../alignments/biogenesis_small_sa_og_norRNA/small_sa_og_norRNA_alignments_piRNA.txt", header = FALSE, stringsAsFactors = FALSE)
  YRNA <- read.table("../alignments/biogenesis_small_sa_og_norRNA/small_sa_og_norRNA_alignments_YRNA.txt", header = FALSE, stringsAsFactors = FALSE)
  snoRNA <- read.table("../alignments/biogenesis_small_sa_og_norRNA/small_sa_og_norRNA_alignments_snoRNA.txt", header = FALSE, stringsAsFactors = FALSE)
  miRNA <- read.table("../alignments/biogenesis_small_sa_og_norRNA/small_sa_og_norRNA_alignments_miRNA.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
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
  read_counts <- read_counts[!(rownames(read_counts) %in% snoRNA$V1),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "Simple"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "Alu"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "5S"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "7SLRNA"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "7SK"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "other"]),]
  
  return(read_counts)
}




#create a fasta file from a list of sequences
createpiRNAFasta <- 
  function(fileName, extension){
    
    #read in the old file
    fileText <- readLines(con  =fileName)
    
    newText <-vector(length = length(fileText)*2)
    
    #change name
    n <- 1
    for(line in seq(1,length(newText),2)){
      #the fasta sequence to be inverted
      sequence <- fileText[n]
      #go to each letter and switch it
      newText[line] <- paste(">",sequence, sep = "") 
      n <- n + 1
    }
    
    n <- 1
    for(line in seq(2,length(newText),2)){
      sequence <- fileText[n]
      #go to each letter and switch it
      newText[line] <- sequence
      n <- n + 1
    }
    
    return(newText)
    #write.table(derp, "test.fa", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }



#trim reads to only the desired length
trimFastaTo <- function(fileName, from, to){
  
  require(Biostrings)
  #read in the old file
  fasta <- readDNAStringSet(fileName)
  
  fasta <- stackStrings(fasta,from=from,to=to,Lpadding.letter = "A", Rpadding.letter = "A")
  
  return(fasta)
  
}


removeSRP040525Repeats <- function(read_counts){
  
  
  snoRNA <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_snoRNA_list.txt", header = FALSE, stringsAsFactors = FALSE)
  
  rmsk <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_repeatmasker.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  
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
  read_counts <- read_counts[!(rownames(read_counts) %in% snoRNA$V1),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "HY"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "Simple"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "Alu"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "5S"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "7SLRNA"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "7SK"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "other"]),]
  
  return(read_counts)
}




plotCoexpression <- 
  function(coexpression){
    
    # create and save the coexpression plot
    print(
      
      #create the correlation plot 
      ggplot(data = coexpression, aes(y = gene_logFC, 
                                      x = piRNA_logFC)) +
        # geom_rect(aes(xmin = 2, xmax =ceiling(max(coexpression$piRNA_logFC)), ymax = -2, ymin = floor(min(coexpression$gene_logFC))), fill = "white", color = "blue", alpha = 1/20000, size = .7) + 
        # geom_rect(aes(xmax = -2, xmin = floor(min(coexpression$piRNA_logFC)), ymin = 2, ymax = ceiling(max(coexpression$gene_logFC))),fill = "white", color = "blue", alpha = 1/20000, size = .7) + 
        geom_point(aes(color = gene_pval, 
                       y = gene_logFC,
                       x = piRNA_logFC, 
                       alpha = 1-piRNA_pval), 
                   size = .9) +
        scale_color_gradientn(name ="Gene Pval" , 
                              colors = c("darkblue","lightblue","red"),
                              breaks = c( .01,.1,1), 
                              limits = c(0,1), 
                              values = c(.0000001, .1, 1)) +
        scale_alpha_continuous(labels = c(1,.1,.05,.001), 
                               limits = c(.001,  1),
                               breaks = c( .001,.05,.1,1),
                               name = "piRNA Pval") + 
        geom_hline(aes(yintercept = 2), 
                   linetype ="dotted") +
        geom_hline(aes(yintercept = -2), 
                   linetype = "dotted") +
        geom_vline(aes(xintercept = 2), 
                   linetype ="dotted") +
        geom_vline(aes(xintercept = -2), 
                   linetype ="dotted") +
        scale_y_continuous(breaks = c(seq((floor(min(coexpression$gene_logFC))), ceiling(max(coexpression$gene_logFC)), 1))) +
        scale_x_continuous(breaks = c(seq((floor(min(coexpression$piRNA_logFC))), ceiling(max(coexpression$piRNA_logFC)), 1))) +
        #scale_linetype_manual(guide = FALSE) +
        
        guides(linetype = FALSE) +   
        #ggtitle("3UTR Targeting: Gene LFC ~ piRNA LFC ") +
        #scale_color_manual(values = c("black", "blue")) + 
        xlab("piRNA logFC") +
        ylab("Target logFC") +
        #scale_size(range =c(.5,1.5)) + 
        # labs(  color = "Gene Pval") +
        theme_bw() 
      #theme(panel.grid.major.x = element_blank()) 
      #theme(legend.position ="right")
      
    )
  }


# construct the coepression data aframe for use later
# personal use code
#
buildCoexpression <- 
  function(alignments, #converted before hand,
           piRNA.expression,
           gene.expression){
    
    coexpression <- 
      unique(
        merge(x = piRNA.expression,
              y = alignments,
              by = "piRNA"))
    
    coexpression <-merge(x = coexpression,
                         y = gene.expression,
                         by = "gene")
    return(coexpression)
    
  }
