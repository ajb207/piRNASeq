




loadmiRDeep <- function(file){
  
  read_file <- file
  
  mirdeep <- readLines(file)
  start <- grep(mirdeep, pattern = "provisional id")
  mirdeep <- mirdeep[start:length(mirdeep)]
  end <-  which(mirdeep == "")[1]
  mirdeep <- mirdeep[1:end]
  mirdeep <- unlist(strsplit(mirdeep, split = "\t"))
  
  novel <- matrix(mirdeep[-c(1:17)], byrow = TRUE, ncol = 17 , dimnames = list(NULL,mirdeep[1:17]))
  
  novel <- subset(novel, select = c(`mature read count`, `consensus mature sequence`))
  
   
  
  
  mirdeep <- readLines(file)
  start <- grep(mirdeep, pattern = "tag id")
  mirdeep <- mirdeep[start:length(mirdeep)]
  end <-  which(mirdeep == "")[1]
  mirdeep <- mirdeep[1:end]
  mirdeep <- unlist(strsplit(mirdeep, split = "\t"))
  
  known <- matrix(mirdeep[-c(1:17)], byrow = TRUE, ncol = 17 , dimnames = list(NULL,mirdeep[1:17]))
  
  known <- subset(known, select = c(`mature read count`, `mature miRBase miRNA`))
  
  
  colnames(novel) <- c("count", 'miRNA')
  colnames(known) <- c("count", 'miRNA')
  
  counts <- rbind(novel, known)
  counts <- as.data.frame(counts)
  counts$count <- as.numeric(as.character(counts$count))
  counts$miRNA <- as.character(counts$miRNA)
  counts <- aggregate(count ~ miRNA, data = counts, FUN = mean, na.rm = TRUE)
  return(counts)
}

files <-  list.files("../mirdeep_files", pattern = ".csv", full.names = TRUE)

mirdeep <- lapply(X = files, FUN =loadmiRDeep)


names(mirdeep) <- c("387++","456++","3359++","3691++","3832++","4121++","4302++","H2S++","387-","456-", "3359-","3691-","3832-","4121-","4302-","H2S-", "1123++","528++","1123-","528-","hNP1_old","hNP1_new","17231","16157")

test <- loadmiRDeep()

Reduce(function(...) merge(..., all=T, by = "miRNA"), mirdeep)


df <- mirdeep[[1]]
colnames(df) <- c("miRNA", "387++")

for(n in 2:24){
  temp <- mirdeep[[n]]
  colnames(temp) <- c("miRNA", names(mirdeep)[n])
  df <- merge(df, temp, by = "miRNA", all = T)
}


rownames(df) <- df$miRNA
df[is.na(df)] <- 0

df <- subset(df, select = -miRNA)

write.table(df, "../sRNA_counts/small_sa_og_norRNA_18_24_miRDeep_counts_table.txt", quote = FALSE, sep = ",")
