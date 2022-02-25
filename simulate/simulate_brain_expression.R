#quants <- read.delim("ERR188297/quant.sf", string=FALSE)
cortex_tpm <- read.delim("frontal_cortex.tsv.gz", skip=2)
# average over the 10 samples
quants <- data.frame(Name = cortex_tpm$transcript_id,
                     TPM = rowMeans(quants[,-c(1,2)]))
# rounding difference
quants$TPM <- quants$TPM * 1e6 / sum(quants$TPM)

suppressPackageStartupMessages(library(GenomicFeatures))
#gtf <- "gencode.v26.annotation.gtf.gz"
txdb.filename <- "gencode.v26.annotation.sqlite"
#txdb <- makeTxDbFromGFF(gtf)
#saveDb(txdb, txdb.filename)
txdb <- loadDb(txdb.filename)

txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]
all(quants$Name %in% txdf$TXNAME)
quants$GENEID <- txdf$GENEID[match(quants$Name, txdf$TXNAME)]
quants <- quants[order(quants$GENEID),]

ebt <- exonsBy(txdb, by="tx", use.names=TRUE)
txp.len <- sum(width(ebt))
quants$Length <- txp.len[quants$Name]
quants$NumReads <- quants$TPM * quants$Length
quants$NumReads <- round(quants$NumReads * 50e6 / sum(quants$NumReads))
table(quants$NumReads > 10)

quant.tpm <- quants$TPM
names(quant.tpm) <- quants$Name
# threshold < 10 reads to 0 TPM
quant.tpm[quants$NumReads < 10] <- 0
quant.len <- quants$Length

# http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.transcripts.fa.gz
fasta0 <- "gencode.v26.transcripts.fa"
fasta <- "gencode.v26.transcripts.short.names.fa"
suppressPackageStartupMessages(library(Biostrings))
txseq <- readDNAStringSet(fasta0)
names(txseq) <- sub("\\|.*","",names(txseq))
writeXStringSet(txseq, file=fasta)
#gc <- as.numeric(letterFrequency(txseq, "GC", as.prob=TRUE))


### simulate DTU and DTE

set.seed(1)
# ratio of expressed genes with DGE, DTE or DTU
dge.ratio <- 0.1
dte.ratio <- 0
dtu.ratio <- 0.1
exprs.txs <- names(quant.tpm)[quant.tpm > 0]
exprs.genes <- unique(txdf$GENEID[match(exprs.txs, txdf$TXNAME)])
dge.genes <- sample(exprs.genes, length(exprs.genes) * dge.ratio)
dte.genes <- sample(setdiff(exprs.genes, dge.genes), length(exprs.genes) * dte.ratio)
# exclude single isoform genes for DTU
exprs.genes.multi <- setdiff(exprs.genes, txdf$GENEID[txdf$ntx == 1])
dtu.genes <- sample( setdiff(exprs.genes.multi, union(dge.genes, dte.genes)),
                    length(exprs.genes) * dtu.ratio)
intersect(dge.genes, dte.genes) # should be null set
intersect(dge.genes, dtu.genes) # should be null set
intersect(dte.genes, dtu.genes) # should be null set

# re-order txdf for speed
txdf <- txdf[match(names(quant.tpm),txdf$TXNAME),]
all.equal(names(quant.tpm), txdf$TXNAME)

# this will keep the isoform TPMs for the two groups
tpms <- matrix(nrow=length(quant.tpm), ncol=2)
rownames(tpms) <- names(quant.tpm)
tpms[,1] <- quant.tpm
tpms[,2] <- quant.tpm
iso.dge <- logical(length(quant.tpm))
iso.dte <- logical(length(quant.tpm))
iso.dte.only <- logical(length(quant.tpm))
iso.dtu <- logical(length(quant.tpm))
names(iso.dge) <- names(quant.tpm)
names(iso.dte) <- names(quant.tpm)
names(iso.dte.only) <- names(quant.tpm)
names(iso.dtu) <- names(quant.tpm)

sampleOne <- function(x) if (length(x) == 1) x else sample(x,1)
sampleUpToTwo <- function(x) if (length(x) <= 2) x else sample(x,2)

if (FALSE) {
  
  # DGE
  for (i in seq_along(dge.genes)) {
    if (i %% 100 == 0) cat(i)
    this.txs <- which(txdf$GENEID == dge.genes[i])
    tpm.this.txs <- quant.tpm[this.txs]
    iso.dge[this.txs] <- TRUE
    coinflip <- sample(c(FALSE,TRUE),1)
    fc <- runif(1,2,6)
    if (coinflip) {
      tpms[this.txs,2] <- fc * tpm.this.txs
    } else {
      tpms[this.txs,1] <- fc * tpm.this.txs
    }
  }

  # DTE
  for (i in seq_along(dte.genes)) {
    if (i %% 100 == 0) cat(i)
    this.txs <- which(txdf$GENEID == dte.genes[i])
    tpm.this.txs <- quant.tpm[this.txs]
    exprs.tx <- sampleOne(which(tpm.this.txs > 0))
    iso.dte[this.txs][exprs.tx] <- TRUE
    # if only one expressed txp, then it's DTE and also not DTU
    if (sum(tpm.this.txs > 0) == 1) {
      iso.dte.only[this.txs][exprs.tx] <- TRUE
    }
    coinflip <- sample(c(FALSE,TRUE),1)
    #fc <- runif(1,2,6)
    fc <- runif(1,2,3)
    if (coinflip) {
      tpms[this.txs,2][exprs.tx] <- fc * tpm.this.txs[exprs.tx]
    } else {
      tpms[this.txs,1][exprs.tx] <- fc * tpm.this.txs[exprs.tx]
    }
  }

}
  
# DTU
for (i in seq_along(dtu.genes)) {
  if (i %% 100 == 0) cat(i)
  this.txs <- which(txdf$GENEID == dtu.genes[i])
  tpm.this.txs <- quant.tpm[this.txs]
  high.tx <- unname(which.max(tpm.this.txs))
  # try to pick a non-zero TPM txp
  low.tx <- if (all(tpm.this.txs[-high.tx] == 0)) {
    sampleUpToTwo(which(tpm.this.txs == 0))
  } else {
    #sampleOne(which(tpm.this.txs > 0 & tpm.this.txs < max(tpm.this.txs)))
    sampleUpToTwo(which(tpm.this.txs > 0 & tpm.this.txs < max(tpm.this.txs)))
  }
  iso.dtu[this.txs][c(low.tx,high.tx)] <- TRUE
  for (lt in low.tx) {
    tpms[this.txs,2][lt] <- tpm.this.txs[high.tx] / length(low.tx)
  }
  tpms[this.txs,2][high.tx] <- sum(tpm.this.txs[low.tx])
}

table(iso.dte,iso.dtu,iso.dge)

length(dte.genes)
2 * length(dtu.genes)
table(diff=tpms[,1] != tpms[,2], exprs=quant.tpm > 0)
table(diff=tpms[,1] != tpms[,2], iso.dte)
table(diff=tpms[,1] != tpms[,2], iso.dtu)

# MA plot of true, underlying diffs
cols <- ifelse(iso.dge, "purple", ifelse(iso.dte, "blue", ifelse(iso.dtu, "red", "black")))
cex <- .2
plot(log2((tpms[,1] + 1) * (tpms[,2] + 1)),
     log2((tpms[,2] + 1) / (tpms[,1] + 1)) + runif(nrow(tpms),-.1,.1),
     col=cols, cex=cex)
abline(h=log2(c(1/6,1/2,2,6)),lty=2)

# make counts from TPMs
fraglen <- 200
lib.size <- 40e6
per.nuc <- tpms * pmax(quant.len - fraglen, 1)
sim.counts.mat <- t(t(per.nuc) / colSums(per.nuc) * lib.size)
# now scale the counts so the median ratio is again 1
med.ratio <- median(sim.counts.mat[,2]/sim.counts.mat[,1], na.rm=TRUE)
sim.counts.mat[,2] <- sim.counts.mat[,2] / med.ratio
median(sim.counts.mat[,2]/sim.counts.mat[,1], na.rm=TRUE)

# only simulate when expected count 5 or more in one group
keep <- sim.counts.mat[,1] >= 5 | sim.counts.mat[,2] >= 5
sim.counts.mat <- sim.counts.mat[keep,]

# subset the transcripts
all(rownames(sim.counts.mat) %in% names(txseq))
txseq <- txseq[ rownames(sim.counts.mat) ]
stopifnot(all(names(txseq) == rownames(sim.counts.mat)))
writeXStringSet(txseq, "data/transcripts.fa")

# the format for polyester
fold_changes <- matrix(1,nrow=nrow(sim.counts.mat),ncol=2)
fold_changes[,2] <- sim.counts.mat[,2]/sim.counts.mat[,1]
# numerical tolerance -- if fold_changes[,2] is very close to 1, set it to 1
fold_changes[ abs(fold_changes[,2] - 1) < 1e-6 ,2] <- 1
rownames(fold_changes) <- rownames(sim.counts.mat)

# dealing with zero denominator
zero.denom <- sim.counts.mat[,1] == 0
fold_changes[zero.denom,1] <- 0
fold_changes[zero.denom,2] <- 1
sim_counts <- sim.counts.mat[,1]
sim_counts[zero.denom] <- sim.counts.mat[zero.denom,2]

alfc <- abs(log2(fold_changes[,2]))
hist(log2(fold_changes[alfc < 3,2]), breaks=100, col="grey")
hist(log2(fold_changes[alfc < 3 & alfc > 0,2]), breaks=100, col="grey")

# draw the dispersion from the joint distribution of GEUVADIS data
sim.means <- rowMeans(sim.counts.mat)
load("meanDispPairs/meanDispPairs.rda")
match.idx <- sapply(sim.means, function(mu) {
 which.min(abs(mu - meanDispPairs$baseMean))
})
disps <- meanDispPairs$dispGeneEst[match.idx]
plot(sim.means, disps, log="xy", cex=.1) # looks right

# log10 disp ~ N(-1,0.5)
#set.seed(1)
#disps <- 10^rnorm(nrow(sim.counts.mat),-1,0.5)

summary(disps)
summary(disps[sim.means > 1000])

library(alpine)
load("data/fitpar_all.rda")
fitpar[[1]][["model.params"]]$gc.knots <- seq(from=.4, to=.6, length=3)
fitpar[[1]][["model.params"]]$gc.bk <- c(0,1)
# plotGC(fitpar, "all", col=rep(1:2,each=15))
frag_GC_bias <- plotGC(fitpar, "all", return.type=2)
frag_GC_bias <- frag_GC_bias[,16:30]
#plot(0:100/100, frag_GC_bias[,1], type="l", ylim=c(0,1))
#for (i in 2:15) {
#  lines(0:100/100, frag_GC_bias[,i])
#}

save(tpms, dge.genes, dte.genes, dtu.genes,
     iso.dge, iso.dte, iso.dte.only, iso.dtu,
     txdf, sim.counts.mat, sim_counts,
     fold_changes, disps, frag_GC_bias,
     file="data/simulate.rda")

library(devtools)
si <- session_info()$packages
write.table(si, file="session_info.txt", quote=FALSE, sep="\t")

