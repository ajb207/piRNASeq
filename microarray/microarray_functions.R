
#creates an expression set from microarray data for analysis with limma
createExpressionSet <- function(projectFolder ){
  # Load data
  ## location of your probe and qc files
  dataFile <- paste(projectFolder, "Data/SampleProbeProfile.txt", sep="")
  qcFile <-   paste(projectFolder,"Data/ControlProbeProfile.txt", sep="")
  ## create an expressionset
  eset <- readBeadSummaryData(dataFile=dataFile, qcFile=qcFile,
                              ProbeID="PROBE_ID", controlID="ProbeID",
                              skip=0, qc.skip=0,
                              annoCols=c("SYMBOL", "DEFINITION", "SYNONYMS", "CHROMOSOME", "ILMN_GENE", "SEARCH_KEY"))
  
  #plot(eset, what = "lumi::boxplot", logMode = FALSE)
  
  # Normalize, transform
  eset <- normaliseIllumina(eset, method="quantile", transform= "log2")
  
  # Some of limma's stuff downstream doesn't work with whatever kind
  # of object that you get out of normaliseIllumina().
  # I coerce it to an "ExpressionSet" and everything seems to work fine.
  class(eset) <- "ExpressionSet"
  
  return(eset)
  
  # Analyze with limma
  ## Make a design matrix
  ## Make a contrast matrix
  ## analyze the normal way: lmFit(), contrasts.fit(), eBayes(), topTable()
}



#probe_conversion <- read.delim("../gene_symbol_probe_id_gene_name.txt", stringsAsFactors=FALSE)

#probe_to_gene <- function(top_table){

#  top_table$probes <- rownames(top_table)
#  return(merge(top_table, probe_conversion, by.x = "probes", by.y = "probe_id", all.x=TRUE))
#}

getGeneResults <- function(results, fit, line){
  
  significant_rows <- apply(abs(results), MARGIN =  1, FUN = sum) >= 1 
  probe_conversion <- subset(select = c(ProbeID, SYMBOL),
                             x = topTable(fit[significant_rows,],
                                          number = "all",
                                          lfc = lfc, 
                                          p.value = .05)
  )
  results <- as.data.frame(results)
  results$ProbeID <- rownames(results)
  results_table <- merge(results, probe_conversion, by = "ProbeID")
  rownames(results_table) <- results_table$ProbeID
  results_table <- subset(results_table, select = -ProbeID)
  results <- subset(results, select = -ProbeID)
  colnames(results_table) <- c(paste(colnames(results), line, sep ="_"), "SYMBOL")
  return(results_table)
}



