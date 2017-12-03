#######################################################################
#
# import-beadstudio.R
# Stephen Turner, December 2014
#
# Ask the core to export text file data with the info below.
# Assumes control probe file has for each sample:
# ProbeID, AVG_Signal, BEAD_STDERR, Avg_NBEADS, Detection Pval.
# Assumes sample probe file has for each sample:
# ProbeID, Symbol, AVG_Signal, BEAD_STDERR, Avg_NBEADS, Detection Pval
# And other annotation columns:
# SEARCH_KEY, ILMN_GENE, CHROMOSOME, DEFINITION, SYNONYMS
#######################################################################

# Load libraries
library(knitr)
library(edgeR)
library(biomaRt)
library(ggplot2)
library(gplots)
library(beadarray)
library(limma)
library(illuminaHumanv4.db)
library("simpleaffy")
library("affyQCReport")
library(annotate)



createExpressionSet <- function(projectFolder){
# Load data
## location of your probe and qc files
dataFile <- paste(projectFolder, "Data/SampleGeneProfile.txt", sep="")
qcFile <-   paste(projectFolder,"Data/ControlGeneProfile.txt", sep="")
## create an expressionset
eset <- readBeadSummaryData(dataFile=dataFile, qcFile=qcFile,
                            ProbeID="SYMBOL", controlID="TargetID",
                            skip=0, qc.skip=0,
                            annoCols=c("SYMBOL", "DEFINITION", "SYNONYMS", "CHROMOSOME", "ILMN_GENE", "SEARCH_KEY"))

# Optional: Annotate the samples (example)
## Manually (bad)
#pData(eset)$condition <- factor(rep(c("ctl", "trt"), each=3))
## Better / more reproducible to do this by importing a csv/table than doing it manually. You've been warned.
#pData <- read.csv("data/metadata.csv", header=TRUE, row.names=1)

# Optional: I use Illumina's annotation. You can annotate the probes yourself if you want. 
# See http://www.bioconductor.org/help/workflows/annotation/annotation/

# Optional: Remove probes that aren't annotated with a gene
#annotated <- !is.na((fData(eset)$SYMBOL))
#table(annotated)
#eset <- eset[annotated,]
#rm(annotated)

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




makeVoomFitWQW <- function(counts, designArray){
  
  celltype <- factor(designArray)
  design <-model.matrix(~celltype)
  
  #converting to DGEList and keeping only expressions significantly above 0
  dgeCounts <- DGEList(counts=counts)
  
  #annotate_ENSG_TopTable(dgeCounts$counts)
  
  #dge.isexpr <- rs > quantile(rs,probs=.25)
  dge.isexpr <- rowSums(cpm(dgeCounts)>1) >=1
  dgeCounts <- dgeCounts[dge.isexpr,]
  
  #normalization
  dgeCountsNormalization <- calcNormFactors(dgeCounts, method = c("TMM"))
  #norm.dge.1187 <- dge.1187
  

  #voom with wuailty weights. default below
  voomCounts <-vooma(dgeCountsNormalization, design, plot=TRUE)
  fit <- eBayes(lmFit(voomCounts,design))
  
  return(fit)
  
}


annotateProbeToSymbol <- function(probeDF){

  probeDF$probe_id <- rownames(probeDF)
  return(merge(probeDF, toTable(illuminaHumanv4SYMBOL[]),  by = "probe_id"))
}




designArray=c("A","A","A","A", "B","B","B","B", "B","B","B","B")
#designArray=c("A","A","A","A", "B","B","B","B", "C","C","C","C")

fit_c3691L <- makeVoomFitWQW(createExpressionSet("Cell_Line_3691/"), designArray)
fit_c3832L <- makeVoomFitWQW(createExpressionSet("Cell_Line_3832/"), designArray)
fit_c4121L <- makeVoomFitWQW(createExpressionSet("Cell_Line_4121/"), designArray)




#voom fold changes
voom_all_genes_c3691L <- topTable(fit_c3691L, sort.by="P", number =nrow(fit_c3691L$p.value) )
voom_topFC_c3691L <- topTable(fitWQWc3691L, sort.by="P", number =nrow(fitWQWc3691L$p.value),lfc = 0.5, p.value = 0.05)

voom_all_genes_c3832L <-annotateProbeToSymbol(topTable(fitWQWc3832L, sort.by="P", number =nrow(fitWQWc3832L$p.value) ))
VoomWQWTopFoldChangesc3832L <- annotateProbeToSymbol(topTable(fitWQWc3832L, sort.by="P", number =nrow(fitWQWc3832L$p.value),lfc = 0.05, p.value = 0.05))

voom_all_genes_c4121L <- annotateProbeToSymbol(topTable(fitWQWc4121L, sort.by="P", number =nrow(fitWQWc4121L$p.value) ))
VoomWQWTopFoldChangesc4121L <- annotateProbeToSymbol(topTable(fitWQWc4121L, sort.by="P", number =nrow(fitWQWc4121L$p.value),lfc = 0.05, p.value = 0.05))


#voomNormExp <- voomCountsWQW$E

 

