## METHODS
# *chipversion* data was imported into GenomeStudio (v , Illumina, ... , CA).
# Missing bead types were imputed using the Gene Expression module in Genome Studio.
# The control probe and sample probe data was exported before normalization and background correction steps.
# 
#
#
#
## The HumanHT v4.0 chip design includes:
#   
#   A: WT serum starved = baseline
# B: WT serum released
# C: Mutant serum starved = baseline
# D: Mutant serum released
# 
# Comparisons include:
#   
#   1- (A vs. C)
# 2- (B vs. D)
# 3- (A vs. B) vs. (C vs. D)
# 4- (A vs. C) vs. (B vs. D)
# 5- (A+B) vs. (C+D)


#
#### FOLLOW THE NUMBERS #1,2,3..


# source("http://bioconductor.org/biocLite.R")
# biocLite("beadarray")
# biocLite("lumi")
# biocLite("MetabricExpression")
rm(list=ls())
library(limma);
library(lumi) #1
library(beadarray)

# ?beadarray
# limmaUsersGuide()

setwd('/Users/yehial/Dropbox/Lamis_Lumi_1')

# http://www.gettinggeneticsdone.com/2014/12/importing-illumina-beadarray-data-into-r.html
#######################################################################
# import-beadstudio.R
# Ask the core to export text file data with the info below.
# Assumes control probe file has for each sample:
# ProbeID, AVG_Signal, BEAD_STDERR, Avg_NBEADS, Detection Pval.
# Assumes sample probe file has for each sample:
# ProbeID,  Symbol, AVG_Signal, BEAD_STDERR, Avg_NBEADS, Detection Pval
# And other annotation columns:
# SEARCH_KEY, ILMN_GENE, CHROMOSOME, DEFINITION, SYNONYMS
#######################################################################

# Load data
# library(beadarray)
## location of your probe and qc files
dataFile <- "data/SampleProbeProfile.txt"
qcFile <- "data/ControlProbeProfile.txt"  # select all that goes under control probe profile on illumina
## create an expressionset
eset <- readBeadSummaryData(dataFile=dataFile, qcFile=qcFile,ProbeID="ProbeID", 
                            controlID="ProbeID",skip=8, qc.skip=8,
                            annoCols=c("SYMBOL", "DEFINITION", "SYNONYMS", "CHROMOSOME", "ILMN_GENE", "SEARCH_KEY"))

#class(eset) <- "ExpressionSet"
#plot(eset, what='boxplot' )


# Optional: Annotate the samples (example)
## Manually (bad)
#pData(eset)$condition <- factor(rep(c("ctl", "trt"), each=3))
## Better / more reproducible to do this by importing a csv/table than doing it manually. You've been warned.
pData <- read.csv("data/humanht-12_sample sheet template.csv", header=TRUE, row.names=1,skip=8)

# remove these samples / optional
 pData<-pData[!row.names(pData) %in% c("WT R2","WT SS3","MUT SS3"),]
# Optional: I use Illumina's annotation. You can annotate the probes yourself if you want.
# See http://www.bioconductor.org/help/workflows/annotation/annotation/

# Optional: Remove probes that aren't annotated with a gene
annotated <- !is.na((fData(eset)$SYMBOL))
table(annotated)
eset <- eset[annotated,]
rm(annotated)

#remove these / optional
eset<- eset[, !row.names(pData(eset)) %in% c("WT R2","WT SS3","MUT SS3")]

# class(eset) <- "ExpressionSet"
# pdf("Figure1.boxplot.of.microarray.intensity.no_norm.pdf")
# plot(eset, what='boxplot',col="yellow" )
# dev.off()

# Normalize, transform
eset <- normaliseIllumina(eset, method="quantile", transform= "log2")

# Some of limma's stuff downstream doesn't work with whatever kind
# of object that you get out of normaliseIllumina().
# I coerce it to an "ExpressionSet" and everything seems to work fine.
class(eset) <- "ExpressionSet"

# Analyze with limma
## Make a design matrix
## Make a contrast matrix
## analyze the normal way: lmFit(), contrasts.fit(), eBayes(), topTable()
pdf("Figure2.boxplot.of.microarray.intensity.pdf",6,6)
#par(mar=c(6,2,2,2)) 
plot(eset, what='boxplot',col=blues9[2])
dev.off()
# LIMMA

pdf("Figure3.Multidimensional.scaling.of.normalized.data.pdf")
plotMDS(eset,col=c(rep("green",3),rep("red",3),rep("navy",3),rep("black",3)))
dev.off()



# OUTPUT Expression values 
summary(eset)
exp_profiles<-exprs(eset)
dim(exp_profiles)
write.csv(exp_profiles,"data/normalized_exp_profiles.csv")


# THBS1	5810685	396.1	1	NaN	10.974	25	0.00000	356.1	1	NaN	14.073	13	0.00000	295.1	1	NaN	5.916	31	0.00000	1233.4	1	NaN	31.425	25	0.00000	1635.9	1	NaN	44.766	17	0.00000	1434.6	1	NaN	30.613	22	0.00000	403.0	1	NaN	9.598	13	0.00000	413.5	1	NaN	12.733	22	0.00000	386.4	1	NaN	11.205	25	0.00000	1249.6	1	NaN	44.741	25	0.00000	1410.9	1	NaN	57.710	16	0.00000	1033.7	1	NaN	20.824	27	0.00000		7057	NM_003246.2	THBS1	ILMN_1686116	
# DKK1	4880609	89.1	1	NaN	2.355	13	0.98442	91.3	1	NaN	2.987	16	0.72468	91.5	1	NaN	3.425	21	0.18701	92.9	1	NaN	2.505	19	0.12727	117.9	1	NaN	4.917	14	0.00260	101.2	1	NaN	4.116	16	0.03636	98.1	1	NaN	4.132	16	0.82597	92.7	1	NaN	3.537	18	0.96494	94.5	1	NaN	2.360	19	0.87273	107.4	1	NaN	4.050	16	0.05974	116.2	1	NaN	4.999	19	0.00130	109.9	1	NaN	3.147	22	0.00130		22943	NM_012242.2	DKK1	ILMN_2174189	
# DKK1	5570102	510.8	1	NaN	12.130	36	0.00000	484.9	1	NaN	15.939	24	0.00000	380.0	1	NaN	7.895	25	0.00000	1409.7	1	NaN	45.512	29	0.00000	1833.5	1	NaN	53.483	24	0.00000	1311.8	1	NaN	43.292	21	0.00000	536.7	1	NaN	12.763	35	0.00000	440.0	1	NaN	10.551	26	0.00000	446.5	1	NaN	13.979	27	0.00000	1753.2	1	NaN	57.399	22	0.00000	1633.0	1	NaN	53.783	24	0.00000	1438.0	1	NaN	38.094	27	0.00000		22943	NM_012242.2	DKK1	ILMN_1773337	
# EGR1	870338	793.0	1	NaN	19.689	27	0.00000	873.9	1	NaN	29.494	30	0.00000	1141.6	1	NaN	25.176	31	0.00000	853.4	1	NaN	28.611	31	0.00000	895.6	1	NaN	28.158	24	0.00000	1080.1	1	NaN	25.975	35	0.00000	2166.5	1	NaN	37.095	27	0.00000	1752.6	1	NaN	40.656	31	0.00000	1349.5	1	NaN	37.844	24	0.00000	1054.6	1	NaN	27.111	31	0.00000	985.3	1	NaN	26.501	24	0.00000	930.7	1	NaN	21.176	34	0.00000		1958	NM_001964.2	EGR1	ILMN_1762899



barplot(exp_profiles["870338", ],beside=TRUE)



genes_of_interest<-exp_profiles[c("5810685", "4880609", "5570102","870338"), ]
row.names(genes_of_interest)<-c("THBS1", "DKK1 (p1)", "DKK1", "EGR1")
pdf("Figure6_Serum_response_diff.pdf",13,8)
barplot(genes_of_interest[c(4),],beside=TRUE)
legend(1,10.5,row.names(genes_of_interest)[c(4)],fill=c("black","grey"))
dev.off()


design <- model.matrix(~0+factor(pData$Sample_Group))
colnames(design)<-sub(" ","_",substr(colnames(design),27,35))
colnames(design)

#An appropriate design matrix can be created and a linear model fitted using
fit <- lmFit(eset, design)


# The HumanHT v4.0 chip design includes:
#   
# A: WT serum starved = baseline
# B: WT serum released
# C: Mutant serum starved = baseline
# D: Mutant serum released
# 
# Comparisons include:
#   
# 1- (A vs. C)
# 2- (B vs. D)
# 3- (A vs. B) vs. (C vs. D)
# 4- (A vs. C) vs. (B vs. D)
# 5- (A+B) vs. (C+D)

contrast.SS       <- makeContrasts(WT_SS-MUT_SS,levels=design) # 1

contrast.R        <- makeContrasts(WT_R-MUT_R,levels=design) # 2

contrast.R_SS     <- makeContrasts(WT_R-MUT_R, WT_SS-MUT_SS,levels=design) # 4 
contrast.WT2_MUT2 <- makeContrasts( (WT_R+WT_SS) - (MUT_R+MUT_SS),levels=design)  #5
contrast.WT_MUT      <- makeContrasts( (WT_R-WT_SS) , (MUT_R - MUT_SS),levels=design) #3


# 3. which genes respond differently in mutant compared to wild-type cells.
contrast.2by2      <- makeContrasts( (WT_R-WT_SS) - (MUT_R - MUT_SS),levels=design) #3


fit2.SS       <- contrasts.fit(fit,contrast.SS       )
fit2.R        <- contrasts.fit(fit,contrast.R        )
fit2.WT_MUT   <- contrasts.fit(fit,contrast.WT_MUT   )
fit2.R_SS     <- contrasts.fit(fit,contrast.R_SS     )
fit2.2by2     <- contrasts.fit(fit,contrast.2by2  )



fit2.SS       <- eBayes( fit2.SS      )
fit2.R        <- eBayes( fit2.R       )
fit2.WT_MUT   <- eBayes( fit2.WT_MUT  )
fit2.R_SS     <- eBayes( fit2.R_SS    )
fit2.WT2_MUT2 <- eBayes( fit2.WT2_MUT2)
fit2.2by2     <- eBayes( fit2.2by2 )


#results <- decideTests(fit2.R_SS,adjust.method = "BH",lfc = 0.5)

results <- decideTests(fit2.WT_MUT,adjust.method = "BH",lfc = 1)
pdf("Figure4.venndiagram.of.wt_R.vs.wt_ss.vs.mut_R.vs.mut_ss_lfc1.pdf",10,10)
vennDiagram(results)
dev.off()

#(WT_R - WT_SS) , (MUT_R - MUT_SS)
set1 <-topTable(fit2.WT_MUT ,coef = 1, adjust="BH",number=nrow(fit2.WT_MUT))
set2 <-topTable(fit2.WT_MUT ,coef = 2, adjust="BH",number=nrow(fit2.WT_MUT))

write.csv(  set1        ,"topTable_WT_R.vs.WT_SS.csv")
write.csv(  set2        ,"topTable_MUT_R.vs.MUT_SS.csv")

set1[setdiff(x=set1$ProbeID,y=set2$ProbeID),] # set1 only
set2[setdiff(y=set1$ProbeID,x=set2$ProbeID),] # set2 only

set1[intersect(x=set1$ProbeID,y=set2$ProbeID),] # set1 only


 

# A list of top genes differential expressed in group2 versus group1 can be obtained from
topTable_SS       <-topTable(fit2.SS       , adjust="BH",number=nrow(fit2.SS      ))
topTable_R        <-topTable(fit2.R        , adjust="BH",number=nrow(fit2.R       ))
topTable_WT_MUT   <-topTable(fit2.WT_MUT   , adjust="BH",number=nrow(fit2.WT_MUT  ))
topTable_R_SS     <-topTable(fit2.R_SS     , adjust="BH",number=nrow(fit2.R_SS    ))
topTable_WT2_MUT2 <-topTable(fit2.WT2_MUT2 , adjust="BH",number=nrow(fit2.WT2_MUT2))

#topTable_WT2_MUT2 <-
#
topTable(fit2.2by2, adjust="BH")

topTable_fit2.2by2 <-topTable(fit2.2by2      , adjust="BH",number=nrow(fit2.2by2))




#write.csv(  topTable_SS       ,"topTable_SS.csv")
write.csv(  topTable_SS       ,"topTable_SS_without3rdReplicate.csv")

#write.csv(  topTable_R        ,"topTable_R.csv")
write.csv(  topTable_R        ,"topTable_R_withoutWTR2.csv")
write.csv(  topTable_WT_MUT   ,"topTable_WT_MUT.csv")
write.csv(  topTable_R_SS     ,"topTable_R_SS.csv")
write.csv(  topTable_WT2_MUT2 ,"topTable_WT2_MUT2.csv")
write.csv(  topTable_fit2.2by2 ,"topTable_fit2.2by2_9samples.csv")


topTable_fit2.2by2



# The outcome of each hypothesis test can be assigned using
results <- decideTests(fit2)
summary(results)
# A Venn diagram showing numbers of genes significant in each comparison can be obtained from
vennDiagram(results)


 
 
dim(NOR)
setwd("~/Research/Herceptin/")

write.csv(NOR,     "rma_limmaNOR.csv")
write.csv(POST,    "rma_limmaPOST.csv")
write.csv(PCR,     "rma_limmaPCR.csv")
write.csv(BASELINE,"rma_limmaBASELINE.csv")

## get significant gene list with FDR adjusted p.values less than 0.01
write.csv(NOR[abs(NOR[,4])>2,],     "rma_limmaNOR_FC2.csv") # 11
write.csv(POST[abs(POST[,4])>2,],    "rma_limmaPOST_FC2.csv") # 231
write.csv(PCR[abs(PCR[,4])>2,],     "rma_limmaPCR_FC2.csv") # 108
write.csv(BASELINE[abs(BASELINE[,4])>2,],"rma_limmaBASELINE_FC2.csv") # 39

#39 /Users/gurkan/Research/Herceptin/rma_limmaBASELINE_FC2.csv
#11 /Users/gurkan/Research/Herceptin/rma_limmaNOR_FC2.csv
#108 /Users/gurkan/Research/Herceptin/rma_limmaPCR_FC2.csv
#231 /Users/gurkan/Research/Herceptin/rma_limmaPOST_FC2.csv
library(VennDiagram)
#pdf("Venn_NOR_BASE_POST_PCR.pdf",9,9)

venn.plot <- venn.diagram(list(NOR=NOR[abs(NOR[,4])>2,1], POST=POST[abs(POST[,4])>2,1],
                               PCR=PCR[abs(PCR[,4])>2,1], BASELINE=BASELINE[abs(BASELINE[,4])>2,1]),
                          filename="Venn_NOR_BASE_POST_PCR.png",
                          height = 750, width = 750,cex=1,resolution = 72,);
#dev.off()
#####
gNOR=NOR[abs(NOR[,4])>2,1]
gPOST=POST[abs(POST[,4])>2,1]
gPCR=PCR[abs(PCR[,4])>2,1]
gBASELINE=BASELINE[abs(BASELINE[,4])>2,1]

library(BiocGenerics)
setdiff(gNOR,union(gPOST,union(gBASELINE,gPCR)))

tT <- NOR
tT <- BASELINE
tT <- POST
tT <- PCR
head(tT)
#Then plot
x <- -log10(tT$adj.P.Val)
pdf("sig.line.fit2NOR.pdf",10,10)
pdf("sig.line.fit2POST.pdf",10,10)
pdf("sig.line.fit2PCR.pdf",10,10)
pdf("sig.line.fit2BASELINE.pdf",10,10)
plot(x, type="l")
sigline <- c(.05, .01, .005, .001,.0005, .0001)
sigline <- c(-log10(sigline),50)
sigcolors <- c("red", "blue", "green", "yellow","pink","purple","black")
sapply(1:length(sigline), function(x){abline(h=sigline[x], col=sigcolors[x])})
dev.off()

# GSEA:
# write(dataMatrix2, file = "dataMatrix2.txt", append = FALSE, sep = "\t")
#   
# write.table(as.matrix(t(s)),"dataMatrix2_annotation.cls",sep="\t")
# write(dataMatrix2, "dataMatrix2.txt" ,sep="\t")
# 
# write.table(cbind(row.names(dataMatrix2),row.names(dataMatrix2),row.names(dataMatrix2)),
#             "dataMatrix2.chp" ,sep="\t")



#### Part tres ######
#GSEA Results showed significant 



# BASELINE GENES:
baseline_genes<-BASELINE[abs(BASELINE[,4])>2,]$ProbeID
#find them in the matrix
table(dataMatrix2probes[,2] %in% baseline_genes) #57 probes

library(gplots)
baseline_gene_expression<- dataMatrix2[dataMatrix2probes[,2] %in% baseline_genes,] 
row.names(baseline_gene_expression)<-a[row.names(baseline_gene_expression),2]
pdf("baseline_genes_expression.pdf",10,10)
heatmap.2(baseline_gene_expression[,grep("PCR", colnames(baseline_gene_expression))],
          cexCol=0.9,cexRow=0.9,margins=c(10,10),trace="none",col=bluered)
dev.off()

pcr_nor_baseline_genes<-c("TFF1","KRT14","KRT17","KRT6B","SERPINA5","SERPINA3","AGT","ESR1",
                          "C4ORF7","KRT5","ACTG2","ENPP5","KLK5","GFRA1","LOC650517","ATP6V1B1",
                          "TMPRSS3","KLK5","HS.388347","C20ORF114","LAMC2","KCNF1","C11ORF70","FAM3B",
                          "ANXA8","GSTA1","KRT16","SERHL2","TRPV6","PHLDA2","TRPV6","FAM5C","ABCC11",
                          "CNTNAP2","TMEM45B","SOX11","ZP2","BMP7","DCD")
# PATHWAY GENES:
pathway_genes <- pcr_nor_baseline_genes
length(pathway_genes)
pathway_genes<-c("DHDDS","IDI2","PDSS2","GGPS1","HMGCS2","ACAT1","HMGCS1",
                 "MVK","IDI1","PMVK","HMGCR","PDSS1","ACAT2","MVD","FDPS")
table(dataMatrix2probes[,2] %in% pathway_genes) #57 probes
pathway_genes_expression<- dataMatrix2[dataMatrix2probes[,2] %in% pathway_genes,] 
#row.names(pathway_genes_expression)<-annot[row.names(pathway_genes_expression),2]
toHeatMap<-aggregate(pathway_genes_expression[,grep("PCR_baseline|NOR_baseline", 
                                                    colnames(pathway_genes_expression))],
                     by=list(row.names(pathway_genes_expression)),FUN=mean)
row.names(toHeatMap)<-toHeatMap$Group.1
toHeatMap$Group.1<-NULL
dim(toHeatMap)
pdf("PSR_NOR_genes_expression_All.pdf",12,11)
heatmap.2(as.matrix(pathway_genes_expression),cexCol=0.5,cexRow=0.9,margins=c(12,10),trace="none",col=bluered)
dev.off()

# TODO:
1. Focus on Baseline and post only datapoints. 
2. look at diff. expr. genes.
3. find groups that can distinguish groups pcr/objr/nor




### FOCUS
source("http://bioconductor.org/biocLite.R")
biocLite("ReactomePA")
library(ReactomePA)

post_enriched<-enrichPathway(gene=gPOST,pvalueCutoff=0.05,qvalueCutoff=0.05,readable=T);
summary(post_enriched)




