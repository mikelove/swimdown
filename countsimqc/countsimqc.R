load("countsimqc_data.rda")
load("countsimqc_geuvadis.rda")
#pkg <- c("rmarkdown", "edgeR", "DESeq2", "dplyr", "tidyr", "ggplot2",
#         "SummarizedExperiment", "genefilter", "DT", "GenomeInfoDbData",
#         "caTools", "randtests")
#pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#devtools::install_github("csoneson/countsimQC")

library(countsimQC)
suppressPackageStartupMessages(library(DESeq2))
cts <- round(cts)
dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=samps,
                              design=~condition)
countsimQCReport(ddsList = list("sim"=dds, "GEUVADIS"=geuvadis),
                 outputFile = "countsimReport.html",
                 outputDir = "./",
                 description = "Comparison of simulated transcript-level counts",
                 maxNForCorr = 50,
                 maxNForDisp = 50,
                 seed = 1,
                 calculateStatistics=FALSE)
