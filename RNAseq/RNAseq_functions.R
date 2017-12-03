

source('global.R')

#line count to get number of reads


lineCount <- function(linecounts){
  #rename columns
  colnames(linecounts) <- c("lines", "file")
  #every forth line is a sequence, so divide by 4
  linecounts[,1] <- linecounts[,1]/4
  linecounts[,2] <- gsub("_.*", "", linecounts[,2])
  return(linecounts)
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
importStarCounts <- 
  function(folder, column){
    
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
    
    
    
    return(final_frame)
  }
