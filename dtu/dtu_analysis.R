load("txi_scaled.rda") # scaledTPM
cts <- txi.list[["salmon"]]$counts
load("data/simulate.rda")
all(rownames(cts) %in% txdf$TXNAME)
txdf <- txdf[match(rownames(cts),txdf$TXNAME),]
all(rownames(cts) == txdf$TXNAME)

# define the truth

txdf <- txdf[match(rownames(tpms), txdf$TXNAME),]
txdf$dtu.genes <- iso.dtu | iso.dte & !iso.dte.only
full.dtu.genes <- unique(txdf$GENEID[txdf$dtu.genes])

txp.exprs <- rowSums(tpms) > 0
dtu.dte.genes <- unique(txdf$GENEID[iso.dte & !iso.dte.only])
txdf$full.dtu <- iso.dtu | (txdf$GENEID %in% dtu.dte.genes & txp.exprs)
dtu.txps <- txdf$TXNAME[txdf$full.dtu]

suppressPackageStartupMessages({
  library(DRIMSeq)
  library(DEXSeq)
  library(stageR)
  library(iCOBRA)
  library(ggplot2)
  library(RColorBrewer)
})

source("dtu_plots.R")

dir <- "res_scaled"

for (n.sub in c(3,6,9,12)) {

  print(n.sub)

  load(paste0("drim_scaled/drim_",n.sub,".rda"))
  res <- DRIMSeq::results(d)
  res$pvalue <- ifelse(is.na(res$pvalue), 1, res$pvalue)
  res.txp <- DRIMSeq::results(d, level="feature")
  res.txp$pvalue <- ifelse(is.na(res.txp$pvalue), 1, res.txp$pvalue)

  # filtering on the SD of proportions
  smallProportionSD <- function(d, filter=.1) {
    cts <- as.matrix(subset(counts(d), select=-c(gene_id, feature_id)))
    gene.cts <- rowsum(cts, counts(d)$gene_id)
    total.cts <- gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
    props <- cts/total.cts
    propSD <- sqrt(rowVars(props))
    propSD < filter
  }
  filt <- smallProportionSD(d)
  res.txp.filt <- res.txp
  res.txp.filt$pvalue[filt] <- 1
  res.txp.filt$adj_pvalue[filt] <- 1
  
  # stageR analysis for DRIMSeq

  strp <- function(x) substr(x,1,15)
  pScreen <- res$pvalue
  names(pScreen) <- strp(res$gene_id)
  pConfirmation <- matrix(res.txp$pvalue, ncol=1)
  rownames(pConfirmation) <- strp(res.txp$feature_id)
  tx2gene <- res.txp[, c("feature_id", "gene_id")]
  for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])
  stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                        pScreenAdjusted=FALSE, tx2gene=tx2gene)
  stageRObj <- stageWiseAdjustment(object=stageRObj, method="dtu", alpha=0.05)
  suppressWarnings({
    drim.padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=TRUE)
  })
  drim.padj$dtu <- drim.padj$txID %in% strp(dtu.txps)

  # stageR analysis for DRIMSeq-filt

  strp <- function(x) substr(x,1,15)
  pScreen <- res$pvalue
  names(pScreen) <- strp(res$gene_id)
  pConfirmation <- matrix(res.txp.filt$pvalue, ncol=1)
  rownames(pConfirmation) <- strp(res.txp.filt$feature_id)
  tx2gene <- res.txp[, c("feature_id", "gene_id")]
  for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])
  stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                        pScreenAdjusted=FALSE, tx2gene=tx2gene)
  stageRObj <- stageWiseAdjustment(object=stageRObj, method="dtu", alpha=0.05)
  suppressWarnings({
    drim.filt.padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=TRUE)
  })
  drim.filt.padj$dtu <- drim.filt.padj$txID %in% strp(dtu.txps)
  
  # DEXSeq results

  load(paste0("dex_scaled/dex_",n.sub,".rda"))

  dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
  qval <- perGeneQValue(dxr)
  dxr.g <- data.frame(gene=names(qval),qval)

  # stageR analysis for DEXSeq

  library(stageR)
  strp <- function(x) substr(x,1,15)
  pConfirmation <- matrix(dxr$pvalue,ncol=1)
  dimnames(pConfirmation) <- list(strp(dxr$featureID),"transcript")
  pScreen <- qval
  names(pScreen) <- strp(names(pScreen))
  tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])
  for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])
  stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                        pScreenAdjusted=TRUE, tx2gene=tx2gene)
  stageRObj <- stageWiseAdjustment(object=stageRObj, method="dtu", alpha=0.05)
  suppressWarnings({
    dex.padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=TRUE)
  })
  dex.padj$dtu <- dex.padj$txID %in% strp(dtu.txps)

  # SUPPA

  suppa <- read.delim(paste0("suppa/diff_empirical_threshold_",n.sub,".dpsi"))
  names(suppa) <- c("txp.gene","dpsi","pval")
  suppa$gene <- sub(";.*", "", suppa$txp.gene)
  suppa$txp <- sub(".*;", "", suppa$txp.gene)
  # subset to those filtered by DRIMSeq to help with multiple testing for SUPPA?
  suppa <- suppa[match(res.txp$feature_id, suppa$txp),]
  suppa$padj <- p.adjust(suppa$pval, method="BH")
  suppa.dxr <- as(DataFrame(groupID=suppa$gene,
                            pvalue=suppa$pval,
                            padj=rep(1, nrow(suppa))), "DEXSeqResults")
  qval <- perGeneQValue(suppa.dxr)
  suppa.g <- data.frame(gene=names(qval), qval=qval)

  # stageR analysis for SUPPA2
  
  suppa$pval <- ifelse(is.na(suppa$pval), 1, suppa$pval)
  pConfirmation <- matrix(suppa$pval,ncol=1)
  dimnames(pConfirmation) <- list(strp(suppa$txp),"transcript")
  pScreen <- suppa.qval
  names(pScreen) <- strp(names(pScreen))
  tx2gene <- as.data.frame(suppa[,c("txp", "gene")])
  for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])
  stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                        pScreenAdjusted=TRUE, tx2gene=tx2gene)
  stageRObj <- stageWiseAdjustment(object=stageRObj, method="dtu", alpha=0.05)
  suppressWarnings({
    suppa.padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=TRUE)
  })
  suppa.padj$dtu <- suppa.padj$txID %in% strp(dtu.txps)

  ## RATs

  load(paste0("rats/rats_boot_",n.sub,".rda"))
  rats.b.txp <- rats_boot$Transcripts
  rats.b <- rats_boot$Genes
  # remove the NA rows by filtering to DRIMSeq genes and txps (same as given to RATs)
  # (RATs puts rows of NA for genes/txps in the `annot` but not in counts)
  rats.b.txp <- rats.b.txp[match(res.txp$feature_id, rats.b.txp$target_id),]
  rats.b <- rats.b[match(res$gene_id, rats.b$parent_id),]
  # some txps and genes are post-hoc filtered by RATS based on 3 criteria.
  # we need to generalize the 'DTU' column for any padj cutoff to use iCOBRA.
  # so set the 'pval' and 'pval_corr' column to 1 for the RATs post-hoc filtered txps and genes
  rats.b.txp$pval[!(rats.b.txp$elig_fx & rats.b.txp$quant_reprod & rats.b.txp$rep_reprod)] <- 1
  rats.b.txp$pval_corr[!(rats.b.txp$elig_fx & rats.b.txp$quant_reprod & rats.b.txp$rep_reprod)] <- 1
  rats.b$pval[!(rats.b$elig_fx & rats.b$quant_reprod & rats.b$rep_reprod)] <- 1
  rats.b$pval_corr[!(rats.b$elig_fx & rats.b$quant_reprod & rats.b$rep_reprod)] <- 1
  
  # stageR analysis for RATs

  pConfirmation <- matrix(rats.b.txp$pval,ncol=1)
  dimnames(pConfirmation) <- list(strp(rats.b.txp$target_id),"transcript")
  pScreen <- rats.b$pval_corr
  names(pScreen) <- strp(rats.b$parent_id)
  tx2gene <- as.data.frame(rats.b.txp[,c("target_id", "parent_id")])
  for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])
  stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                        pScreenAdjusted=TRUE, tx2gene=tx2gene)
  stageRObj <- stageWiseAdjustment(object=stageRObj, method="dtu", alpha=0.05)
  suppressWarnings({
    rats.b.padj <- getAdjustedPValues(stageRObj, order=FALSE, onlySignificantGenes=TRUE)
  })
  rats.b.padj$dtu <- rats.b.padj$txID %in% strp(dtu.txps)

  print("making plots")
  
  # gene-level screening
  padj <- data.frame(row.names=txdf$TXNAME)
  padj.gene <- data.frame(row.names=unique(txdf$GENEID))
  
  padj.gene$DRIMSeq <- res$adj_pvalue[match(rownames(padj.gene), res$gene_id)]
  padj.gene$DEXSeq <- dxr.g$qval[match(rownames(padj.gene), dxr.g$gene)]
  padj.gene$SUPPA2 <- suppa.g$qval[match(rownames(padj.gene), suppa.g$gene)]
  padj.gene$RATs <- rats.b$pval_corr[match(rownames(padj.gene), rats.b$parent_id)]
  truth <- data.frame(status=as.numeric(rownames(padj.gene) %in% full.dtu.genes),
                      row.names=rownames(padj.gene))

  nas <- apply(padj.gene, 1, function(x) all(is.na(x)))
  padj.gene <- padj.gene[!nas,]
  truth <- truth[!nas,,drop=FALSE]

  cd.gene <- COBRAData(padj=padj.gene, truth=truth)
  pdf(file=paste0(dir,"/dtu_gene_",n.sub,".pdf"))
  print(myicobra(cd.gene, "Gene", c(1,2,4)))
  dev.off()

  # txp-level
  padj$DRIMSeq <- res.txp$adj_pvalue[match(rownames(padj), res.txp$feature_id)]
  padj$DRIMSeq.filt <- res.txp.filt$adj_pvalue[match(rownames(padj), res.txp.filt$feature_id)]
  padj$DEXSeq <- dxr$padj[match(rownames(padj), dxr$featureID)]
  padj$SUPPA2 <- suppa$padj[match(rownames(padj), suppa$txp)]
  padj$RATs <- rats.b.txp$pval_corr[match(rownames(padj), rats.b.txp$target_id)]
  truth <- data.frame(status=as.numeric(rownames(padj) %in% dtu.txps),
                      row.names=rownames(padj))

  nas <- apply(padj, 1, function(x) all(is.na(x)))
  padj <- padj[!nas,]
  truth <- truth[!nas,,drop=FALSE]

  cd <- COBRAData(padj=padj, truth=truth)
  pdf(file=paste0(dir,"/dtu_txp_",n.sub,".pdf"))
  print(myicobra(cd, "Transcript"))
  dev.off()
  
  # stageR txp-level confirmation
  alpha <- .05
  tot.true <- sum(rownames(padj) %in% dtu.txps)
  tprs <- 1/tot.true * c(sum(drim.padj$transcript[drim.padj$dtu] < alpha, na.rm=TRUE),
                         sum(drim.filt.padj$transcript[drim.padj$dtu] < alpha, na.rm=TRUE),
                         sum(dex.padj$transcript[dex.padj$dtu] < alpha, na.rm=TRUE),
                         sum(rats.b.padj$transcript[rats.b.padj$dtu] < alpha, na.rm=TRUE),
                         sum(suppa.padj$transcript[suppa.padj$dtu] < alpha, na.rm=TRUE))
  ofdrs <- c(ofdr(drim.padj,alpha=alpha),
             ofdr(drim.filt.padj,alpha=alpha),
             ofdr(dex.padj,alpha=alpha),
             ofdr(rats.b.padj,alpha=alpha),
             ofdr(suppa.padj,alpha=alpha))
  methods <- c("DRIMSeq", "DRIMSeq-filt", "DEXSeq", "RATs", "SUPPA2")
  (df <- data.frame(method=methods,
                    TPR=tprs, OFDR=round(ofdrs,3)))
  write.table(df,row.names=FALSE, quote=FALSE, sep="\t",
              file=paste0(dir,"/OFDR_",n.sub,".tsv"))

  # FDR breakdown plots

  pdf(file=paste0(dir,"/dtu_gene_fdrbrk_",n.sub,".pdf"), width=9)
  fdrBreakdownGene(padj.gene, alpha=.05)
  dev.off()

  pdf(file=paste0(dir,"/dtu_txp_fdrbrk_",n.sub,".pdf"), width=9)
  fdrBreakdownTxp(padj, alpha=.05)
  dev.off()

  # TP breakdown plot

  pdf(file=paste0(dir,"/dtu_gene_tpbrk_",n.sub,".pdf"), width=9)
  tpBreakdownGene(padj.gene, alpha=.05)
  dev.off()
  
}

### OFDR plots

pdf(file=file.path(dir,"ofdr.pdf"), width=7, height=5)
plotOFDR(with2=FALSE)
dev.off()
