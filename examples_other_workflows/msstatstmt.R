## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(MSstatsTMT)

## -----------------------------------------------------------------------------
# read in PD PSM sheet
# raw.pd <- read.delim("161117_SILAC_HeLa_UPS1_TMT10_5Mixtures_3TechRep_UPSdB_Multiconsensus_PD22_Intensity_PSMs.txt")
head(raw.pd)

# Read in annotation including condition and biological replicates per run and channel.
# Users should make this annotation file. It is not the output from Proteome Discoverer.
annotation.pd <- read.csv(file="PD_Annotation.csv", header=TRUE)
head(annotation.pd) 

# do not remove PSM with missing values within one run
input.pd <- PDtoMSstatsTMTFormat(raw.pd, annotation.pd)
head(input.pd)

# remove PSM with missing values within one run
input.pd.no.miss <- PDtoMSstatsTMTFormat(raw.pd, annotation.pd,
                                         rmPSM_withMissing_withinRun = TRUE)
head(input.pd.no.miss)

## -----------------------------------------------------------------------------
# Read in MaxQuant files
# proteinGroups <- read.table("proteinGroups.txt", sep="\t", header=TRUE)

# evidence <- read.table("evidence.txt", sep="\t", header=TRUE)

# Users should make this annotation file. It is not the output from MaxQuant.
# annotation.mq <- read.csv(file="MQ_Annotation.csv", header=TRUE)

input.mq <- MaxQtoMSstatsTMTFormat(evidence, proteinGroups, annotation.mq)
head(input.mq)

## -----------------------------------------------------------------------------
# Read in SpectroMine PSM report
# raw.mine <- read.csv('20180831_095547_CID-OT-MS3-Short_PSM Report_20180831_103118.xls', sep="\t")

# Users should make this annotation file. It is not the output from SpectroMine
# annotation.mine <- read.csv(file="Mine_Annotation.csv", header=TRUE)

input.mine <- SpectroMinetoMSstatsTMTFormat(raw.mine, annotation.mine)
head(input.mine)

## -----------------------------------------------------------------------------
# read in MSstatsTMT report from OpenMS
# raw.om <- read.csv("OpenMS_20200222/20200225_MSstatsTMT_OpenMS_Export.csv")
head(raw.om)

# the function only requries one input file
input.om <- OpenMStoMSstatsTMTFormat(raw.om)
head(input.om)

## ----message = FALSE , warning = FALSE----------------------------------------
# use MSstats for protein summarization
quant.msstats <- proteinSummarization(input.pd,
                                      method="msstats",
                                      global_norm=TRUE,
                                      reference_norm=TRUE,
                                      remove_norm_channel = TRUE,
                                      remove_empty_channel = TRUE)
head(quant.msstats)

# use Median for protein summarization
# since median method doesn't impute missing values, 
# we need to use the input data without missing values
quant.median <- proteinSummarization(input.pd.no.miss,
                                     method="Median",
                                     global_norm=TRUE,
                                     reference_norm=TRUE,
                                     remove_norm_channel = TRUE,
                                     remove_empty_channel = TRUE)

head(quant.median)

## -----------------------------------------------------------------------------
## Profile plot without norm channnels and empty channels
dataProcessPlotsTMT(data.peptide = input.pd,
                    data.summarization = quant.msstats,
                    type = 'ProfilePlot',
                    width = 21, # adjust the figure width since there are 15 TMT runs.
                    height = 7)

# ## Profile plot with all the channels
# quant.msstats.all <- proteinSummarization(input.pd,
#                                       method="msstats",
#                                       normalization=TRUE,
#                                       remove_norm_channel=FALSE,
#                                       remove_empty_channel=FALSE)
# 
# dataProcessPlotsTMT(data.peptide = input.pd,
#                      data.summarization = quant.msstats.all,
#                      type = 'ProfilePlot',
#                      width = 21, # adjust the figure width since there are 15 TMT runs.
#                      height = 7)

## Quality control plot 
# dataProcessPlotsTMT(data.peptide=input.pd,
# data.summarization=quant.msstats, 
# type='QCPlot',
# width = 21, # adjust the figure width since there are 15 TMT runs. 
# height = 7)

## ----message = FALSE, warning = FALSE-----------------------------------------
# test for all the possible pairs of conditions
test.pairwise <- groupComparisonTMT(quant.msstats)
head(test.pairwise)

# Check the conditions in the protein data
levels(quant.msstats$Condition)
# Only compare condition 0.125 and 1
comparison<-matrix(c(-1,0,0,1),nrow=1)
# Set the names of each row
row.names(comparison)<-"1-0.125"
# Set the column names
colnames(comparison)<- c("0.125", "0.5", "0.667", "1")
comparison

test.contrast <- groupComparisonTMT(data = quant.msstats, contrast.matrix = comparison)
head(test.contrast)
