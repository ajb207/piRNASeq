




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




#clustering reads
read_counts <- read.delim("../sRNA_counts/small_sa_og_norRNA_25_33_counts_table.txt", sep=",", stringsAsFactors=FALSE, row.names = 1)[,]
read_counts <- read_counts[,4:19]
read_counts <- read_counts[,-c(3,11)]



#filter to keep only piRNA that meet requirments  
filter <- read.table("../data/small_sa_og_norRNA_25_33_nodust_noncRNA_list.txt", header = FALSE, stringsAsFactors = FALSE)
read_counts <- read_counts[rownames(read_counts) %in% filter$V1,]
sum(read_counts)




unique_mapping <- read.delim("../data/small_sa_og_norRNA_25_33_list_unique_mapping_genome_v1_m5.txt", header=FALSE, stringsAsFactors=FALSE)
read_counts <- read_counts[rownames(read_counts) %in% unique_mapping$V1,]




#shortening reads

read_counts$short <- strtrim(rownames(read_counts),21)

read_counts <- aggregate(. ~ short, data = read_counts, FUN = sum, na.rm = TRUE)
rownames(read_counts) <- read_counts$short

read_counts <- subset(read_counts, select = -short)

write.table(read_counts, "../sRNA_counts/small_sa_og_norRNA_25_33_21bp_unique_mapping_genome_v1_m5_counts_table.txt", quote = FALSE, sep = ",")


dist_frame <- matrix(nrow = length(read_counts$short), ncol = length(read_counts$short))
dist_frame <- data.frame(dist_frame); rownames(dist_frame) <- read_counts$short; colnames(dist_frame) <- read_counts$short

for(read in read_counts$short){
  dist_frame[,read] <- stringdist(read, read_counts$short, method = "lv",nthread = 1)
  
}


dist_frame[dist_frame > 2] <- 100


clusters <- apcluster(as.matrix(dist_frame),q = 1)
#plot(aggExCluster(s = as.matrix(dist_frame)), ticks = 10)

for(cluster in clusters@clusters){
  if(length(cluster) > 1){
    print(cluster)
  }
}




load("clusters_99q.object")
clusters <- clusters_99q
rm(clusters_95q)

for(read in read_counts$short){
  for(cluster_number in 1:length(clusters@clusters)){
    if(read %in% names(clusters@clusters[[cluster_number]])){
      read_counts$short[read_counts$short == read] <- names(clusters@exemplars[cluster_number])
    }
  }
}


clusterReplace <- function(read){
  for(cluster_number in 1:length(clusters@clusters)){
    if(read %in% names(clusters@clusters[[cluster_number]])){
      #read_counts$short[read_counts$short == read] <- names(clusters@exemplars[[cluster_number]])
      return(names(clusters@clusters[[cluster_number]][1]))
    }
    else{
      return(read)
    }
  }
}

test <- sapply(read_counts$short, FUN = clusterReplace)

read_counts$short <- test

read_counts <- aggregate(. ~ short, data = read_counts, FUN = sum, na.rm = TRUE)
rownames(read_counts) <- read_counts$short

read_counts <- subset(read_counts, select = -short)


write.table(read_counts, "../sRNA_counts/small_sa_og_norRNA_25_33_21bp_unique_mapping_table.txt", sep = ",", quote = FALSE)
