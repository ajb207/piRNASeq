
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
                                                           ifelse(grepl("MIR", rmsk$V2), "MIR", "other" )
                                                           )
                                                    )
                                             )
                                      )
                               )
                        )
  
  rmsk$Family[rmsk$V2 == "5S"] <- "5S"
  rmsk$Family[rmsk$V2 == "7SLRNA"] <- "7SLRNA"
  rmsk$Family[rmsk$V2 == "7SK"] <- "7SK"
  
  
  
  
  
  
  #read_counts <- read_counts[!(rownames(read_counts) %in% tRNA$V1),]
  #read_counts <- read_counts[!(rownames(read_counts) %in% lncRNA$V1),]
  #read_counts <- read_counts[!(rownames(read_counts) %in% YRNA$V1),]
  #read_counts <- read_counts[!(rownames(read_counts) %in% snoRNA$V1),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "Simple"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "Alu"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "5S"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "7SLRNA"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "7SK"]),]
  read_counts <- read_counts[!(rownames(read_counts) %in% rmsk$V1[rmsk$Family == "other"]),]
  
  return(read_counts)
}


#import the desired read counts
read_counts <- read.delim("../../data_files/sRNA_counts/small_sa_og_norRNA_25_33_counts.txt", sep=",", stringsAsFactors=FALSE, row.names = 1)

rRNA <- read.table(file = "../../data_files/stem_diff_alignments/biogenesis_small_sa_og_norRNA/small_sa_og_norRNA_alignments_rRNA.txt", 
                   header = FALSE, 
                   stringsAsFactors = FALSE)
read_counts <- read_counts[!(rownames(read_counts) %in% rRNA$V1),]
YRNA <- read.table(file = "../../data_files/stem_diff_alignments/biogenesis_small_sa_og_norRNA/small_sa_og_norRNA_alignments_YRNA.txt", 
                   header = FALSE, 
                   stringsAsFactors = FALSE)
read_counts <- read_counts[!(rownames(read_counts) %in% YRNA$V1),]
read_counts <- removeRepeats(read_counts)


read_counts <- read_counts[,4:19]
#sort the samples so
read_counts <- read_counts[,c(9:16, 1:8)]

read_counts <- read_counts[,-c(3,11)]


#read_counts <- read_counts[nchar(rownames(read_counts)) >= 28,]

read_counts <- read_counts[rowSums(cpm(read_counts) > 1) >= 7,]

#shorten the reads to just the predicted targeting sequence
read_counts$short <- strtrim(rownames(read_counts),25)
read_counts <- aggregate(. ~ short, 
                         data = read_counts,
                         FUN = sum, 
                         na.rm = TRUE)
rownames(read_counts) <- read_counts$short
read_counts <- subset(read_counts, select = -short)


write.table(x = rownames(read_counts), 
            file = "../data/temp_read_names.txt", 
            quote=F,
            row.names = F,
            col.names = F)

targeting_sRNA <- read.table(file = "../data/targeting_sRNA.csv",
                             stringsAsFactors = F,
                             header  = T, 
                             sep = ",")
#targeting_sRNA <- removeRepeats(targeting_sRNA)
write.table(x = rownames(targeting_sRNA), 
            file = "../data/temp_read_names.txt", 
            quote=F,
            row.names = F,
            col.names = F)


system("grep -f ../data/temp_read_names.txt ../alignments/list_targeting/small_sa_og_norRNA_25_33_25bp_targeting_promoter_bt2v1score2_trimmed.txt > ../data/promoter_temp_alignments_bt2v1score2.txt")
system("grep -f ../data/temp_read_names.txt ../alignments/list_targeting/small_sa_og_norRNA_25_33_25bp_targeting_promoter_bt2v1score4_trimmed.txt > ../data/promoter_temp_alignments_bt2v1score4.txt" )
system("grep -f ../data/temp_read_names.txt ../alignments/list_targeting/small_sa_og_norRNA_25_33_25bp_targeting_promoter_bt2v1score6_trimmed.txt > ../data/promoter_temp_alignments_bt2v1score6.txt" )



#system("grep -f ../data/temp_read_names.txt ../data/repeat_temp_alignments_bt2v1score2.txt > ../data/temp.txt; cat ../data/temp.txt > ../data/repeat_temp_alignments_bt2v1score2.txt")
system("grep -f ../data/temp_read_names.txt ../alignments/list_targeting/small_sa_og_norRNA_25_33_25bp_targeting_repeat_bt2v1score2_trimmed.txt > ../data/repeat_temp_alignments_bt2v1score2.txt")
system("grep -f ../data/temp_read_names.txt ../alignments/list_targeting/small_sa_og_norRNA_25_33_25bp_targeting_repeat_bt2v1score4_trimmed.txt > ../data/repeat_temp_alignments_bt2v1score4.txt" )
system("grep -f ../data/temp_read_names.txt ../alignments/list_targeting/small_sa_og_norRNA_25_33_25bp_targeting_repeat_bt2v1score6_trimmed.txt > ../data/repeat_temp_alignments_bt2v1score6.txt" )



system("grep -f ../data/temp_read_names.txt ../alignments/list_targeting/small_sa_og_norRNA_25_33_25bp_targeting_cDNA_bt2v1score2_trimmed.txt > ../data/cDNA_temp_alignments_bt2v1score2.txt" )
system("grep -f ../data/temp_read_names.txt ../alignments/list_targeting/small_sa_og_norRNA_25_33_25bp_targeting_cDNA_bt2v1score4_trimmed.txt > ../data/cDNA_temp_alignments_bt2v1score4.txt" )
system("grep -f ../data/temp_read_names.txt ../alignments/list_targeting/small_sa_og_norRNA_25_33_25bp_targeting_cDNA_bt2v1score6_trimmed.txt > ../data/cDNA_temp_alignments_bt2v1score6.txt" )
system("grep -f ../data/temp_read_names.txt ../alignments/list_targeting/DE_targeting_cDNA_bt2v1score8.sam | cut -f 1,3,4,6,19 >  ../data/cDNA_temp_alignments_bt2v1score8.txt" )
system("grep -f ../data/temp_read_names.txt ../alignments/list_targeting/DE_targeting_cDNA_bt2v1score10.sam | cut -f 1,3,4,6,19 >  ../data/cDNA_temp_alignments_bt2v1score10.txt" )


system("grep -f ../data/temp_read_names.txt ../alignments/list_targeting/small_sa_og_norRNA_25_33_25bp_targeting_3UTR_bt2v1score2_trimmed.txt > ../data/3UTR_temp_alignments_bt2v1score2.txt" )
system("grep -f ../data/temp_read_names.txt ../alignments/list_targeting/small_sa_og_norRNA_25_33_25bp_targeting_3UTR_bt2v1score4_trimmed.txt > ../data/3UTR_temp_alignments_bt2v1score4.txt" )
system("grep -f ../data/temp_read_names.txt ../alignments/list_targeting/small_sa_og_norRNA_25_33_25bp_targeting_3UTR_bt2v1score6_trimmed.txt > ../data/3UTR_temp_alignments_bt2v1score6.txt" )

