txi.list <- list()
n <- 12
sample <- paste0(rep(1:n,2),"_",rep(1:2,each=n))
dir <- "quant"
library(tximport)
meths <- c("salmon","kallisto")
for (m in meths) {
  filename <- ifelse(m == "salmon","quant.sf","abundance.h5")
  files <- file.path(dir,m,sample,filename)
  names(files) <- sample
  # for DGE and DTE, some methods don't accept offsets
  txi.list[[m]] <- tximport(files,
                            type=m,
                            txOut=TRUE,
                            countsFromAbundance="lengthScaledTPM",
                            varReduce=TRUE)
}
save(txi.list, file="txi.rda")

txi.list <- list()
for (m in meths) {
  filename <- ifelse(m == "salmon","quant.sf","abundance.h5")
  files <- file.path(dir,m,sample,filename)
  names(files) <- sample
  # use scaledTPM counts for DTU
  txi.list[[m]] <- tximport(files,
                            type=m,
                            txOut=TRUE,
                            countsFromAbundance="scaledTPM",
                            varReduce=TRUE)
}
save(txi.list, file="txi_scaled.rda")

txi.list <- list()
for (m in meths) {
  filename <- ifelse(m == "salmon","quant.sf","abundance.h5")
  files <- file.path(dir,m,sample,filename)
  names(files) <- sample
  # DESeq2 can take raw counts and offset
  txi.list[[m]] <- tximport(files,
                            type=m,
                            txOut=TRUE,
                            varReduce=TRUE)
}
save(txi.list, file="txi_raw.rda")

# small EDA of counts

load("txi.rda")
load("data/simulate.rda")
outdir <- "out"
quantdir <- "quant"

n <- 12
samps <- paste0(rep(1:n,2),"_",rep(1:2,each=n))

ngene <- nrow(fold_changes)
cmat <- matrix(NA, nrow=ngene, ncol=2*n)
for (i in 1:n) {
  print(i)
  # these are counts _before_ fragment GC bias
  load(paste0(outdir,"/out_",i,"/sim_counts_matrix.rda"))
  cmat[,c(i,i+n)] <- counts_matrix
}
rownames(cmat) <- rownames(counts_matrix)
colnames(cmat) <- samps
#save(cmat, file="data/cmat.rda")

counts.lfc <- log2(rowMeans(cmat[,(n+1):(2*n)])+1) -
  log2(rowMeans(cmat[,1:n])+1)
simulated <- log2(fold_changes[,2]+.01) - log2(fold_changes[,1]+.01)

plot(simulated,
     counts.lfc,
     col=rgb(0,0,0,.3), cex=.5, xlim=c(-6,6), ylim=c(-6,6))
abline(0,1, col="red", lwd=2)
abline(v=c(-1,1)*log2(6), lty=2)
