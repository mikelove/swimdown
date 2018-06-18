load("txi.rda")
load("data/simulate.rda")

suppressPackageStartupMessages({
  library(tximport)
  library(EBSeq)
})

txi <- txi.list[["salmon"]]
txi.g <- summarizeToGene(txi, tx2gene=txdf[,2:1], countsFromAbundance="lengthScaledTPM")

n <- 12
n.sub <- 12
idx <- c(1:n.sub, n + 1:n.sub)

st <- system.time({
  ebseq.padj <- ebseq(txi$counts[,idx], n.sub, gene=FALSE)
})

st.gene <- system.time({
  ebseq.padj.gene <- ebseq(txi.g$counts[,idx], n.sub, gene=TRUE)
})

save(ebseq.padj, ebseq.padj.gene, st, st.gene,
     file=paste0("ebseq/ebseq_",n.sub,".rda"))

ebseq <- function(cts, n.sub, gene=TRUE) {
  keep <- rowSums(cts >= 10) >= n.sub
  x <- round(cts[keep,])
  condition <- factor(rep(1:2,each=n.sub))
  sizes <- MedianNorm(x)
  if (gene) {
    suppressMessages({
      res <- EBTest(Data=x, Conditions=condition,
                    sizeFactors=sizes, maxround=5, Print=FALSE)
    })
  } else {
    suppressMessages({
      Ng <- GetNg(txdf[,2], txdf[,1])$IsoformNgTrun
      res <- EBTest(Data=x, NgVector=Ng, Conditions=condition,
                    sizeFactors=sizes, maxround=5, Print=FALSE)
    })
  }
  padj <- rep(1, nrow(x))
  names(padj) <- rownames(x)
  padj[match(rownames(res$PPMat), rownames(x))] <- res$PPMat[,"PPEE"]
  padj
}
