suppressPackageStartupMessages(library(DESeq2))
# http://duffel.rail.bio/recount/v2/ERP001942/rse_gene.Rdata
load("rse_gene.Rdata")
source("my_scale_counts.R")
x <- my_scale_counts(rse_gene)
x <- x[,!duplicated(x$sample)]

# https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/
mdata <- read.table("E-GEUV-1.sdrf.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)
mdata <- mdata[,c("Comment.ENA_RUN.","Characteristics.population.","Performer")]
colnames(mdata) <- c("run","population","center")
table(x$run %in% mdata$run)

x <- x[,x$run %in% mdata$run]
mdata <- mdata[mdata$run %in% x$run,]
mdata <- mdata[!duplicated(mdata$run),]
mdata$center <- factor(mdata$center)
mdata$population <- factor(mdata$population)
cd2 <- merge(colData(x), mdata, by="run")
colData(x) <- cd2

dds <- DESeqDataSet(x, ~center + population)
dds <- dds[rowMeans(counts(dds)) > 5,]
dds <- estimateSizeFactors(dds)

library(parallel)
nworkers <- 6
options(mc.cores=nworkers)
idx <- factor(sort(rep(1:nworkers, length.out=nrow(dds))))

# takes 20 min
dds.list <- mclapply(1:nworkers, function(i) {
  estimateDispersionsGeneEst(dds[idx == i,])
})

dds <- do.call(rbind, dds.list)
#save(dds, file="recount_geuvadis_dds.rda")

#

if (FALSE) {
  suppressPackageStartupMessages(library(DESeq2))
  load("recount_geuvadis_dds.rda")
  with(subset(mcols(dds), dispGeneEst > 1e-4),
       plot(baseMean, dispGeneEst, log="xy"))
  abline(h=10^c(-1.5,-1,-.5),col="red")
  meanDispPairs <- subset(mcols(dds), dispGeneEst > 1e-4, c(baseMean, dispGeneEst))
  save(meanDispPairs, file="meanDispPairs.rda")
}
