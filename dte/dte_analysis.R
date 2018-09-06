load("txi.rda")
load("data/simulate.rda")

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(limma)
  library(edgeR)
  library(samr)
  library(iCOBRA)
  library(ggplot2)
})

# all methods other than sleuth will use Salmon quants
txi <- txi.list[["salmon"]]
txi.g <- summarizeToGene(txi, tx2gene=txdf[,2:1],
                         countsFromAbundance="lengthScaledTPM")

# DESeq2 will use counts + offset
# (edgeR performs similarly with lengthScaledTPM or counts)
load("txi_raw.rda")
txi.raw <- txi.list[["salmon"]]
txi.raw.g <- summarizeToGene(txi.raw, tx2gene=txdf[,2:1])

# total sample size
n <- 12

source("dte_plots.R")

# subset to, e.g. 1-3 vs 13-15 for n.sub=3, etc.
for (n.sub in c(3,6,9,12)) {

  message(n.sub)
  
  idx <- c(1:n.sub, n + 1:n.sub)

  padj <- data.frame(row.names=rownames(txi$counts))
  padj.gene <- data.frame(row.names=rownames(txi.g$counts))

  ### sleuth

  load(paste0("sleuth/sleuth_",n.sub,".rda"))
  
  padj$sleuth <- sleuth.res$qval[match(rownames(padj), sleuth.res$target_id)]
  padj.gene$sleuth <- sleuth.res.gene$qval[match(rownames(padj.gene), sleuth.res.gene$target_id)]

  ### EBSeq

  load(paste0("ebseq/ebseq_",n.sub,".rda"))
  padj$EBSeq <- ebseq.padj[match(rownames(padj), names(ebseq.padj))]
  padj.gene$EBSeq <- ebseq.padj.gene[match(rownames(padj.gene), names(ebseq.padj.gene))]

  ### the rest are faster, can be run in loop:

  source("dte_funcs.R")

  ### limma-voom, edgeR, and edgeR QL

  message("limma, edgeR, edgeR QL")

  st <- system.time({ 
    tt <- limmavoom(txi$counts[,idx], n.sub)
  })
  st.gene <- system.time({ 
    tt.gene <- limmavoom(txi.g$counts[,idx], n.sub)
  })

  padj$limma <- tt$adj.P.Val[match(rownames(padj), rownames(tt))]
  padj.gene$limma <- tt.gene$adj.P.Val[match(rownames(padj.gene), rownames(tt.gene))]

  st <- system.time({ 
    tt <- edger(txi$counts[,idx], n.sub)
  })
  st.gene <- system.time({ 
    tt.gene <- edger(txi.g$counts[,idx], n.sub)
  })

  padj$edgeR <- tt$FDR[match(rownames(padj), rownames(tt))]
  padj.gene$edgeR <- tt.gene$FDR[match(rownames(padj.gene), rownames(tt.gene))]

  st <- system.time({ 
    tt <- edgerQL(txi$counts[,idx], n.sub)
  })
  st.gene <- system.time({ 
    tt.gene <- edgerQL(txi.g$counts[,idx], n.sub)
  })

  padj$edgeR.QL <- tt$FDR[match(rownames(padj), rownames(tt))]
  padj.gene$edgeR.QL <- tt.gene$FDR[match(rownames(padj.gene), rownames(tt.gene))]

  ### DESeq2

  message("DESeq2")

  st <- system.time({
    res <- deseq2(txi.raw, n.sub)
  })
  st.gene <- system.time({
    res.gene <- deseq2(txi.raw.g, n.sub)
  })

  padj$DESeq2 <- res$padj[match(rownames(padj), rownames(res))]
  padj.gene$DESeq2 <- res.gene$padj[match(rownames(padj.gene), rownames(res.gene))]
  
  ### SAMseq

  message("SAMseq")

  st <- system.time({
    sam.padj <- samseq(txi$counts[,idx], n.sub)
  })
  st.gene <- system.time({
    sam.padj.gene <- samseq(txi.g$counts[,idx], n.sub)
  })

  padj$SAMseq <- sam.padj[match(rownames(padj), names(sam.padj))]
  padj.gene$SAMseq <- sam.padj.gene[match(rownames(padj.gene), names(sam.padj.gene))]

  ### iCOBRA

  # gene-level
  status <- as.numeric(rownames(padj.gene) %in% union(dge.genes, dte.genes))
  cd.gene <- COBRAData(padj=padj.gene,
                       truth=data.frame(status, row.names=rownames(padj.gene)))
  pdf(file=paste0("res_dte/dge_",n.sub,".pdf"))
  print(myicobra(cd.gene, "Gene"))
  dev.off()

  # txp-level
  iso.any <- iso.dtu | iso.dte | iso.dge
  status <- as.numeric(rownames(padj) %in% names(iso.dtu)[iso.any])
  cd <- COBRAData(padj=padj,
                  truth=data.frame(status, row.names=rownames(padj)))
  pdf(file=paste0("res_dte/dte_",n.sub,".pdf"))
  print(myicobra(cd, "Transcript"))
  dev.off()

  pdf(file=paste0("res_dte/dge_fdrbkd_",n.sub,".pdf"), width=9)
  fdrBreakdownGene(padj.gene, alpha=.05)
  dev.off()

  pdf(file=paste0("res_dte/dte_fdrbkd_",n.sub,".pdf"), width=9)
  fdrBreakdownTxp(padj, alpha=.05)
  dev.off()
  
}
