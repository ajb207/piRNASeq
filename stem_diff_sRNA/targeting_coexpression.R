

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



conversion <- read.table("../../reference_sets/ensembl_conversion.txt", header = T, stringsAsFactors = F, sep = "\t")[,c(1,2)]
piRNA_expression <- read.table("../data/targeting_sRNA.csv", stringsAsFactors = F, header  = T, sep = ",")
target_expression <- read.csv("../../RNAseq/data/StemVsDiff_all_DE_genes.txt", sep = ",", stringsAsFactors = F, header = T)

alignments <- read.delim("../data/cDNA_temp_alignments_bt2v1score2.txt", header=FALSE, stringsAsFactors=FALSE)
alignments[,2] <- gsub(pattern = "\\..*", replacement = "", alignments[,2])
alignments[,2] <- conversion$gene_id[match(alignments[,2] , table = conversion$transcript_id)]



createCoexpression <- function(alignments,piRNA_expression, target_expression,fileName = NA){
  
  
  colnames(alignments) <- c("piRNA", "target", "bp_start", "alignment_info", "mismatches")
  #alignments$transcript_id <- gsub(pattern = "\\..*", replacement = "", alignments$transcript_id)
  #alignments$transcript_id <- sapply(alignments$target, FUN = function(x){unlist(strsplit(x, split = "\\|"))[2]})
  colnames(piRNA_expression) <- paste("piRNA.", colnames(piRNA_expression), sep = "")
  colnames(target_expression) <- paste("target.", colnames(target_expression), sep = "")
  
  coexpression <- unique(merge(alignments, piRNA_expression, by.y = 'row.names', by.x = 'piRNA'))
  
  coexpression <-merge(coexpression,target_expression, by.x = "target", by.y = 'row.names')

  
  if(is.na(fileName)==F){
    write.table(coexpression, fileName, col.names = T, row.names = T, sep = ",")
  }
  return(coexpression)
}

test <- createCoexpression(alignments, piRNA_expression, target_expression)



#pdf(paste("../../figures/3UTR_coexpression_",name ,".pdf"), 6,4)
print(
  
  #create the correlation plot 
  ggplot(data = coexpression,aes( y = target.logFC, x = piRNA.logFC)) +
    # geom_rect(aes(xmin = 2, xmax =ceiling(max(coexpression$piRNA_logFC)), ymax = -2, ymin = floor(min(coexpression$gene_logFC))), fill = "white", color = "blue", alpha = 1/20000, size = .7) + 
    # geom_rect(aes(xmax = -2, xmin = floor(min(coexpression$piRNA_logFC)), ymin = 2, ymax = ceiling(max(coexpression$gene_logFC))),fill = "white", color = "blue", alpha = 1/20000, size = .7) + 
    geom_point(aes(color = target.adj.P.Val, y = target.logFC, x = piRNA.logFC, alpha = 1-piRNA.adj.P.Val), size = .9) +
    scale_color_gradientn(name ="Target Pval" , colors = c("darkblue","lightblue","red"),breaks = c( .01,.1,1), limits = c(0,1), values = c(.0000001, .1, 1)) +
    scale_alpha_continuous( labels = c(1,.1,.05,.001), limits = c(.001,  1), breaks = c( .001,.05,.1,1), name = "piRNA Pval") + 
    geom_hline(aes(yintercept = 2), linetype ="dotted") +
    geom_hline(aes(yintercept = -2), linetype = "dotted") +
    geom_vline(aes(xintercept = 2), linetype ="dotted") +
    geom_vline(aes(xintercept = -2), linetype ="dotted") +
    scale_y_continuous(breaks = c(seq((floor(min(coexpression$target.logFC))), ceiling(max(coexpression$target.logFC)), 1))) +
    scale_x_continuous(breaks = c(seq((floor(min(coexpression$piRNA.logFC))), ceiling(max(coexpression$piRNA.logFC)), 1))) +
    #scale_linetype_manual(guide = FALSE) +
    
    guides(linetype = FALSE) +   
    #ggtitle("3UTR Targeting: Gene LFC ~ piRNA LFC ") +
    #scale_color_manual(values = c("black", "blue")) + 
    xlab("piRNA LFC") +
    ylab("Target LFC") +
    #scale_size(range =c(.5,1.5)) + 
    # labs(  color = "Gene Pval") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank()) 
  
  
)
#dev.off()

