



#load required libraries

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
require(stringr)
require(scales)
require(knitr)



#annotation of all small RNA

#loading in the read counts
read_counts <- read.delim("../sRNA_counts/SRP040525_counts.txt", sep=",", stringsAsFactors=FALSE, row.names = 1)[,]
colnames(read_counts) <- gsub(colnames(read_counts), pattern ="_sa_collapsed.fq.", replacement =  "") #changing the names of the reads
og <- read.table("../alignments/SRP040525_sa_alignments_genome_list.txt", header = FALSE, stringsAsFactors = FALSE) #reads found on the genome
rRNA <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_rRNA_list.txt", header = FALSE, stringsAsFactors = FALSE) #rRNA reads
#YRNA <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_YRNA_list.txt", header = FALSE, stringsAsFactors = FALSE) #YRNA reads

#remove unwanted reads as defined above
read_counts <- read_counts[rownames(read_counts) %in% og$V1,]
read_counts <- read_counts[!rownames(read_counts) %in% rRNA$V1,]
#read_counts <- read_counts[!rownames(read_counts) %in% YRNA$V1,]
#keep reads in length of piRNAs
read_counts <- read_counts[nchar(rownames(read_counts)) >= 25 & nchar(rownames(read_counts)) <= 33,]
read_counts <- cpm(read_counts)


#the filters are what sequences the piRNA align to
filters <- list.files("../alignments/biogenesis_SRP040525/", full.names = T)
titles <- c("mRNA", "lincRNA", "lncRNA", "miRNA","Extended\npiRNA\nBase",  "piRNA\nBase", "Extended\ntRNA", "Extended\nRepeat\nElements", "tRNA", "Repeat\nElements", "rRNA", "snoRNA",  "Y-RNA")


#create a dataframe of the biogenesis
df <- as.data.frame(matrix(nrow = ncol(read_counts), ncol = 0, byrow = TRUE))

df$cell_line <- colnames(read_counts)

df$cell_type <- factor(rep(c(rep(c("Stem"),2),rep(c("Diff"), 4)),2), c("Stem", "Diff"))
df$cell_line <- rep(c("0", "10", "17FB", "17MN", "35FB", "35MN"),2)

df$total_sRNA <- colSums(read_counts)




all <- c(NA)
#all_filters <- c(NA)
#join all the biogenesis counts together
for(n in c(3,4,10,12)){
  
  #import the desired read counts
  filter <- read.table(filters[n], header = FALSE, stringsAsFactors = FALSE)
  df[,titles[n]] <- colSums(read_counts[rownames(read_counts) %in% filter$V1,])
  
  all <- c(all, filter$V1)
}

df$'No Previous\nAnnotation'  <- colSums(read_counts[!rownames(read_counts) %in% all,])




pdf("../figures/SRP040525_25_33_biogenesis.pdf", 6,4)

  
#melt the data frame into a usable format by ggplot
temp <- melt(data = df[,-c(3)], id.vars = c(1,2))

#temp <- aggregate(data = temp, value ~ cell_line + cell_type+variable, sum)

  ggplot(data = temp ) +
    #ggplot(data = temp) + 
    geom_boxplot(aes(y = value, x = variable, fill = cell_type)) +
    #labs(title = "Small RNA Anotation of 25-33 nt Reads") +
    ylab("Counts Per Million") +
    xlab("Small RNA Annotation") +
    scale_y_continuous(labels = comma) +
    scale_fill_brewer(name = "Phenotype", labels = c("Stem", "Differentiated"), breaks = c("Stem", "Diff")) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank() 
    )
  

dev.off()




#piRNA piRBase only


read_counts <- read.delim("../sRNA_counts/SRP040525_counts.txt", sep=",", stringsAsFactors=FALSE, row.names = 1)[,]
#read_counts <- read_counts[,4:19]
colnames(read_counts) <- gsub(colnames(read_counts), pattern ="_sa_collapsed.fq.", replacement =  "")
og <- read.table("../alignments/SRP040525_sa_alignments_genome_list.txt", header = FALSE, stringsAsFactors = FALSE)
rRNA <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_rRNA_list.txt", header = FALSE, stringsAsFactors = FALSE)
YRNA <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_YRNA_list.txt", header = FALSE, stringsAsFactors = FALSE)
read_counts <- read_counts[rownames(read_counts) %in% og$V1,]
read_counts <- read_counts[!rownames(read_counts) %in% rRNA$V1,]
read_counts <- read_counts[!rownames(read_counts) %in% YRNA$V1,]
read_counts <- read_counts[nchar(rownames(read_counts)) >= 25 & nchar(rownames(read_counts)) <= 33,]

#read_counts

read_counts <- cpm(read_counts)


#the filters are what sequences the piRNA align to
filters <- list.files("../alignments/biogenesis_SRP040525/", full.names = T)
titles <- c("mRNA", "lincRNA", "lncRNA", "miRNA","Extended\npiRBase\nTranscriptome",  "Original\npiRBase", "Extended\ntRNA", "Extended\nRepeat\nElements", "tRNA", "Repeat\nElements", "rRNA", "snoRNA",  "Y-RNA")


#create a dataframe of the biogenesis
df <- as.data.frame(matrix(nrow = ncol(read_counts), ncol = 0, byrow = TRUE))

df$cell_line <- colnames(read_counts)

df$cell_type <- factor(rep(c(rep(c("Stem"),2),rep(c("Diff"), 4)),2), c("Stem", "Diff"))
df$cell_line <- rep(c("0", "10", "17FB", "17MN", "35FB", "35MN"),2)

df$Total <- colSums(read_counts)


#unique <- read.table("../alignments/Sr", header = FALSE, stringsAsFactors = FALSE)
temp <- df
df$piRNA <- "Unique\nMapping piRNA"
#temp$piRNA <- "Repeat\nMapping piRNA"
#df$Total <- colSums(read_counts[rownames(read_counts) %in% unique$V1,])
#temp$Total <- colSums(read_counts[!rownames(read_counts) %in% unique$V1,])

all_unique <- c(NA)
all_repeat <- c(NA)
#join all the biogenesis counts together
for(n in c(5,6)){
  
  #temp_counts <- read_counts[(rownames(read_counts) %in% unique$V1),]
  temp_counts <- read_counts
  
  filter <- read.table(filters[n], header = FALSE, stringsAsFactors = FALSE)
  #print(file)
  column_sums <- as.data.frame(colSums(temp_counts[rownames(temp_counts) %in% filter$V1,]))
  rownames(column_sums) <- colnames(read_counts)
  colnames(column_sums) <- titles[n]
  df <- cbind(df,column_sums)
  all_unique <- c(all_unique, filter$V1)
}

# 
# for(n in c(1,5,10,11)){
#   
#   temp_counts <- read_counts[!(rownames(read_counts) %in% unique$V1),]
#   #temp_counts <- read_counts
#   
#   filter <- read.table(filters[n], header = FALSE, stringsAsFactors = FALSE)
#   #print(file)
#   column_sums <- as.data.frame(colSums(temp_counts[rownames(temp_counts) %in% filter$V1,]))
#   rownames(column_sums) <- colnames(read_counts)
#   colnames(column_sums) <- titles[n]
#   temp <- cbind(df,column_sums)
#   all_repeat <- c(all_unique, filter$V1)
# }
# 
# 
# 
# 
df$'No Previous\nAnnotation' <- colSums(temp_counts[!rownames(temp_counts) %in% all_unique,])
# temp$'No Previous\nAnnotation' <- colSums(temp_counts[!rownames(temp_counts) %in% all_unique,])
# 
# df <- rbind(df, temp)


#  Counts of piRNA from each predicted biogenesis location based on PHENOTYPE
# 
pdf("../figures/SRP040525_Boxplot_piRNA_biogenesis_phenotype_piRBaseonly.pdf", 6, 4)
#



temp <- melt(data = df[,-c(3,4)], id.vars = c(1,2))
#temp <- aggregate(data = temp, value ~ cell_line + cell_type+variable, mean)


ggplot(data = temp) +
  #ggplot(data = temp) + 
  geom_boxplot(aes(y = value, x = variable, fill = cell_type)) +
  #geom_point(aes(shape = cell_line, x = variable, y = value, color = cell_type), position = position_jitterdodge()) +
  # labs(title = "Putative piRNA Secondary Annotation") +
  ylab("Counts Per Million of piRNA") +
  xlab("piRNA Annotation") +
  scale_y_continuous(labels = comma) +
  scale_x_discrete( limits = c("Original\npiRBase", "Extended\npiRBase\nTranscriptome", "No Previous\nAnnotation")) +
  scale_fill_brewer(name = "Phenotype", labels = c("Stem", "Differentiated"), breaks = c("Stem", "Diff")) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()




#piRNA biogenesis without piRBase


read_counts <- read.delim("../sRNA_counts/SRP040525_counts.txt", sep=",", stringsAsFactors=FALSE, row.names = 1)[,]
#read_counts <- read_counts[,4:19]
colnames(read_counts) <- gsub(colnames(read_counts), pattern ="_sa_collapsed.fq.", replacement =  "")
og <- read.table("../alignments/SRP040525_sa_alignments_genome_list.txt", header = FALSE, stringsAsFactors = FALSE)
rRNA <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_rRNA_list.txt", header = FALSE, stringsAsFactors = FALSE)
YRNA <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_YRNA_list.txt", header = FALSE, stringsAsFactors = FALSE)
read_counts <- read_counts[rownames(read_counts) %in% og$V1,]
read_counts <- read_counts[!rownames(read_counts) %in% rRNA$V1,]
read_counts <- read_counts[!rownames(read_counts) %in% YRNA$V1,]
read_counts <- read_counts[nchar(rownames(read_counts)) >= 25 & nchar(rownames(read_counts)) <= 33,]

#read_counts

read_counts <- cpm(read_counts)


#the filters are what sequences the piRNA align to
filters <- list.files("../alignments/biogenesis_SRP040525/", full.names = T)
titles <- c("mRNA", "lincRNA", "lncRNA", "miRNA","Extended\npiRBase",  "piRNA\nBase", "Extended\ntRNA", "Extended\nRepeat\nElements", "tRNA", "Repeat\nElements", "rRNA", "snoRNA",  "Y-RNA")


#create a dataframe of the biogenesis
df <- as.data.frame(matrix(nrow = ncol(read_counts), ncol = 0, byrow = TRUE))

df$cell_line <- colnames(read_counts)

df$cell_type <- factor(rep(c(rep(c("Stem"),2),rep(c("Diff"), 4)),2), c("Stem", "Diff"))
df$cell_line <- rep(c("0", "10", "17FB", "17MN", "35FB", "35MN"),2)

df$Total <- colSums(read_counts)


#unique <- read.table("../alignments/Sr", header = FALSE, stringsAsFactors = FALSE)
temp <- df
df$piRNA <- "Unique\nMapping piRNA"
#temp$piRNA <- "Repeat\nMapping piRNA"
#df$Total <- colSums(read_counts[rownames(read_counts) %in% unique$V1,])
#temp$Total <- colSums(read_counts[!rownames(read_counts) %in% unique$V1,])

all_unique <- c(NA)
all_repeat <- c(NA)
#join all the biogenesis counts together
for(n in c(1,9,12)){
  
  #temp_counts <- read_counts[(rownames(read_counts) %in% unique$V1),]
  temp_counts <- read_counts
  
  filter <- read.table(filters[n], header = FALSE, stringsAsFactors = FALSE)
  #print(file)
  column_sums <- as.data.frame(colSums(temp_counts[rownames(temp_counts) %in% filter$V1,]))
  rownames(column_sums) <- colnames(read_counts)
  colnames(column_sums) <- titles[n]
  df <- cbind(df,column_sums)
  all_unique <- c(all_unique, filter$V1)
}

# 
# for(n in c(1,5,10,11)){
#   
#   temp_counts <- read_counts[!(rownames(read_counts) %in% unique$V1),]
#   #temp_counts <- read_counts
#   
#   filter <- read.table(filters[n], header = FALSE, stringsAsFactors = FALSE)
#   #print(file)
#   column_sums <- as.data.frame(colSums(temp_counts[rownames(temp_counts) %in% filter$V1,]))
#   rownames(column_sums) <- colnames(read_counts)
#   colnames(column_sums) <- titles[n]
#   temp <- cbind(df,column_sums)
#   all_repeat <- c(all_unique, filter$V1)
# }
# 
# 
# 
# 
df$'No Previous\nAnnotation' <- colSums(temp_counts[!rownames(temp_counts) %in% all_unique,])
# temp$'No Previous\nAnnotation' <- colSums(temp_counts[!rownames(temp_counts) %in% all_unique,])
# 
# df <- rbind(df, temp)


#  Counts of piRNA from each predicted biogenesis location based on PHENOTYPE
# 
pdf("../figures/SRP040525_Boxplot_piRNA_biogenesis_phenotype_nopiRBase.pdf", 6, 4)
#



temp <- melt(data = df[,-c(3,4)], id.vars = c(1,2))
#temp <- aggregate(data = temp, value ~ cell_line + cell_type+variable, mean)


ggplot(data = temp) +
  #ggplot(data = temp) + 
  geom_boxplot(aes(y = value, x = variable, fill = cell_type)) +
  #geom_point(aes(shape = cell_line, x = variable, y = value, color = cell_type), position = position_jitterdodge()) +
  # labs(title = "Putative piRNA Secondary Annotation") +
  ylab("Counts Per Million of piRNA") +
  xlab("piRNA Annotation") +
  scale_y_continuous(labels = comma) +
  scale_x_discrete( limits = c("mRNA", "tRNA", "snoRNA", "No Previous\nAnnotation")) +
  scale_fill_brewer(name = "Phenotype", labels = c("Stem", "Differentiated"), breaks = c("Stem", "Diff")) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()



#based on piRNA type and phenotype
# 
# 
# temp <- melt(data = df[,-c(3)], id.vars = c(1,2,3))
# temp <- aggregate(data = temp, value ~ piRNA +cell_line + cell_type+variable, sum)
# temp$piRNA <- factor(temp$piRNA, levels = c("Unique\nMapping piRNA", "Repeat\nMapping piRNA"))
# 
# pdf("../../figures/Boxplot_piRNA_biogenesis_phenotype_and_piRNAtype.pdf", 12, 4)
# 
# 
# ggplot(data = temp) +
#   #ggplot(data = temp) + 
#   geom_boxplot(aes(y = value, x = variable, fill = cell_type)) +
#   facet_grid(.~piRNA) +
#   # labs(title = "Putative piRNA Secondary Annotation") +
#   ylab("Counts Per Million of piRNA") +
#   xlab("piRNA Annotation") +
#   scale_y_continuous(labels = comma) +
#   scale_x_discrete(limits = c("Clusters","mRNA", "tRNA", "snoRNA", "piRBase", "No Previous\nAnnotation")) +
#   scale_fill_brewer(name = "Phenotype", labels = c("GSC", "Non-GSC"), breaks = c("GSC", "Non-GSC")) +
#   #scale_color_manual(name = "piRNA", labels = c("Unique\nMapping", "Repeat\nMapping"), breaks = c("UM", "RM"), values = c("grey","black")) +
#   #scale_linetype_manual(name = "piRNA", labels = c("Unique\nMapping", "Repeat\nMapping"), breaks = c("UM", "RM"), values = c(5,1)) +
#   #scale_size_manual(name = "piRNA", labels = c("Unique\nMapping", "Repeat\nMapping"), breaks = c("UM", "RM"), values = c(.5,1)) +
#   
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 0),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank()
#   )
# 
# dev.off()
# 
# 




#repeat elements

#import the desired read counts

read_counts <- read.delim("../sRNA_counts/SRP040525_counts.txt", sep=",", stringsAsFactors=FALSE, row.names = 1)[,]
#read_counts <- read_counts[,4:19]
colnames(read_counts) <- gsub(colnames(read_counts), pattern ="_sa_collapsed.fq.", replacement =  "")
og <- read.table("../alignments/SRP040525_sa_alignments_genome_list.txt", header = FALSE, stringsAsFactors = FALSE)
rRNA <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_rRNA_list.txt", header = FALSE, stringsAsFactors = FALSE)
#YRNA <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_YRNA_list.txt", header = FALSE, stringsAsFactors = FALSE)
#read_counts <- read_counts[(rownames(read_counts) %in% filter$V1),]
read_counts <- read_counts[rownames(read_counts) %in% og$V1,]
read_counts <- read_counts[!rownames(read_counts) %in% rRNA$V1,]
#read_counts <- read_counts[!rownames(read_counts) %in% YRNA$V1,]
read_counts <- read_counts[nchar(rownames(read_counts)) >= 25 & nchar(rownames(read_counts)) <= 33,]

read_counts <- as.data.frame(cpm(read_counts))

rmsk <- read.table("../alignments/biogenesis_SRP040525/SRP040525_sa_alignments_repeatmasker.txt", header = FALSE, stringsAsFactors = FALSE, sep = "\t")



read_counts$sRNA <- rownames(read_counts)
read_counts <- merge(read_counts, rmsk[,c(1,2)], by.x = "sRNA", by.y = "V1")
read_counts <- aggregate(data =read_counts[,-1], .~V2, FUN = sum)
#denote the family of the repeat element
read_counts$Family <- ifelse(grepl("tRNA", read_counts$V2), "tRNA", 
                             ifelse(grepl(")n", read_counts$V2), "Simple", 
                                    ifelse(grepl("Alu", read_counts$V2), "SINE", 
                                           ifelse(read_counts$V2 %in% c("U1", "U2", "U3", "U4", "U5", "U6", "U7", "U8"), "AU-rich", 
                                                  ifelse(grepl("HY", read_counts$V2), "Y-RNA", 
                                                         ifelse(grepl(pattern =  paste(c("L1", "L2", "L3", "L4", "L5", "L6"), collapse = "|"), read_counts$V2), "LINE", 
                                                                ifelse(grepl("MIR", read_counts$V2), "SINE", "other" )))))))

read_counts$Family[read_counts$V2 == "5S"] <- "other"
read_counts$Family[read_counts$V2 == "7SLRNA"] <- "7SLRNA"
read_counts$Family[read_counts$V2 == "7SK"] <- "other"



#read_counts[,2:17] <- cpm(read_counts[,2:17])

#write.table(rownames(read_counts), "../fastx_files/sa_og_norRNA_25_33_alignments_repeat.txt", quote=F, row.names = F, col.names = T, sep = ",")

repeats <- as.data.frame(matrix(nrow = nrow(read_counts), ncol =0, byrow = T ))
repeats$Family <- read_counts$Family
repeats$Stem <- rowMeans(read_counts[,c(2,3,9,10)])
repeats$Diff <- rowMeans(read_counts[,c(4,5,6,7,9,10,11,12)])

repeats <- aggregate(data =repeats, .~Family, FUN = sum)


pdf("../figures/SRP040525_25_33_repeatelements_bar.pdf", 6,4)
  

temp <- repeats
temp <- melt(temp, 1)
  
  ggplot(data = temp) +
    geom_bar(aes(y = value, x = variable, fill = Family), stat = "identity") +
    #geom_text(aes(x = variable, y = value, label = Family)) +
    #coord_polar(theta = "y") + 
    labs(title = "") +
    ylab("Counts Per Million") +
    xlab("Phenotype") +
    scale_fill_brewer(palette="Spectral") +
    scale_x_discrete(labels = c("Stem", "Differentiated")) +
    scale_y_continuous(labels = comma) +
    #scale_fill_brewer(name = "Phenotype", labels = c("Stem", "Diff"), breaks = c("Stem", "Diff")) +
    guides(fill = guide_legend(reverse=F)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 0),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  

dev.off()


