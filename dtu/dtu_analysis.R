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
})

# DRIMSeq and DEXSeq analyses now saved in respective dirs

dir <- "res_scaled"

for (n.sub in c(3,6,9,12)) {

  print(n.sub)

  load(paste0("drim_scaled/drim_",n.sub,".rda"))
  res <- DRIMSeq::results(d)
  res$pvalue <- ifelse(is.na(res$pvalue), 1, res$pvalue)
  res.txp <- DRIMSeq::results(d, level="feature")
  res.txp$pvalue <- ifelse(is.na(res.txp$pvalue), 1, res.txp$pvalue)

  # try filtering on the SD of proportions
  getSampleProportions <- function(d) {
    cts <- as.matrix(subset(counts(d), select=-c(gene_id, feature_id)))
    gene.cts <- rowsum(cts, counts(d)$gene_id)
    total.cts <- gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
    cts/total.cts
  }
  res.txp.filt <- res.txp
  prop.d <- getSampleProportions(d)
  res.txp.filt$prop.sd <- sqrt(rowVars(prop.d))
  res.txp.filt$pvalue[res.txp.filt$prop.sd < .1] <- 1
  res.txp.filt$adj_pvalue <- p.adjust(res.txp.filt$pvalue, method="BH")

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

  print("making plots")
  
  # iCOBRA
  
  ofdr <- function(padj, alpha=.05) {
    all.not.dtu <- sapply(with(padj, split(dtu, geneID)), function(x) ! any(x))
    sig.not.dtu <- sapply(with(padj, split(transcript < alpha & ! dtu, geneID)), any)
    mean(all.not.dtu | sig.not.dtu)
  }
  
  myicobra <- function(cd, lvl="Gene") {
    cp <- calculate_performance(cd,
                                binary_truth="status",
                                aspect=c("fdrtpr","fdrtprcurve"),
                                thrs=c(.01,.05,.1))
    cobraplot <- prepare_data_for_plot(cp)
    plot_fdrtprcurve(cobraplot, plottype="points",
                     xaxisrange=c(0,max(fdrtpr(cp)$FDR)),
                     yaxisrange=c(0,max(fdrtpr(cp)$TPR)),
                     title=paste0(lvl,"-level screening, n=",n.sub," vs ",n.sub))
  }  
  
  # gene-level screening
  padj <- data.frame(row.names=txdf$TXNAME)
  padj.gene <- data.frame(row.names=unique(txdf$GENEID))
  
  padj.gene$DRIMSeq <- res$adj_pvalue[match(rownames(padj.gene), res$gene_id)]
  padj.gene$DEXSeq <- dxr.g$qval[match(rownames(padj.gene), dxr.g$gene)]
  padj.gene$SUPPA2 <- suppa.g$qval[match(rownames(padj.gene), suppa.g$gene)]
  truth <- data.frame(status=as.numeric(rownames(padj.gene) %in% full.dtu.genes),
                      row.names=rownames(padj.gene))

  nas <- apply(padj.gene, 1, function(x) all(is.na(x)))
  padj.gene <- padj.gene[!nas,]
  truth <- truth[!nas,,drop=FALSE]

  cd.gene <- COBRAData(padj=padj.gene, truth=truth)
  pdf(file=paste0(dir,"/dtu_gene_",n.sub,".pdf"))
  print(myicobra(cd.gene, "Gene", c(1,2,4)))
  dev.off()

  # catching DGE?
  res$dge <- res$gene_id %in% dge.genes
  dxr.g$dge <- dxr.g$gene %in% dge.genes
  suppa.g$dge <- suppa.g$gene %in% dge.genes
  sum(res$adj_pvalue[res$dge] < .05, na.rm=TRUE)
  sum(dxr.g$qval[dxr.g$dge] < .05)
  sum(suppa.g$qval[suppa.g$dge] < .05)

  # stageR txp-level confirmation 
  alpha <- .05
  tps <- c(sum(drim.padj$transcript[drim.padj$dtu] < alpha, na.rm=TRUE),
           sum(drim.filt.padj$transcript[drim.padj$dtu] < alpha, na.rm=TRUE),
           sum(dex.padj$transcript[dex.padj$dtu] < alpha, na.rm=TRUE))
  ofdrs <- c(ofdr(drim.padj,alpha=alpha),
             ofdr(drim.filt.padj,alpha=alpha),
             ofdr(dex.padj,alpha=alpha))
  (df <- data.frame(method=c("DRIMSeq","DRIMSeq-filt","DEXSeq"), TP=tps, OFDR=round(ofdrs,3)))
  write.table(df,row.names=FALSE, quote=FALSE, sep="\t",
              file=paste0(dir,"/OFDR_",n.sub,".tsv"))

  # txp-level
  padj$DRIMSeq <- res.txp$adj_pvalue[match(rownames(padj), res.txp$feature_id)]
  padj$DRIMSeq.filt <- res.txp.filt$adj_pvalue[match(rownames(padj), res.txp.filt$feature_id)]
  padj$DEXSeq <- dxr$padj[match(rownames(padj), dxr$featureID)]
  padj$SUPPA2 <- suppa$padj[match(rownames(padj), suppa$txp)]
  truth <- data.frame(status=as.numeric(rownames(padj) %in% dtu.txps),
                      row.names=rownames(padj))

  nas <- apply(padj, 1, function(x) all(is.na(x)))
  padj <- padj[!nas,]
  truth <- truth[!nas,,drop=FALSE]

  cd <- COBRAData(padj=padj, truth=truth)
  pdf(file=paste0(dir,"/dtu_txp_",n.sub,".pdf"))
  print(myicobra(cd, "Transcript"))
  dev.off()
  
}

if (FALSE) {
  x <- read.table("res_scaled/OFDR_3.tsv", header=TRUE)
  for (i in c(6,9,12)) {
    x <- rbind(x, read.table(paste0("res_scaled/OFDR_",i,".tsv"), header=TRUE))
  }
  x$n <- factor(rep(1:4*3, each=3), levels=1:4*3)
  library(ggplot2)
  library(RColorBrewer)
  pdf(file="res_scaled/ofdr.pdf", width=7, height=5)
  ggplot(x, aes(OFDR, TP, color=method, label=n)) +
    geom_vline(xintercept=0.05, color="gray", linetype=2) + 
    geom_point(size=4) + geom_path(show.legend=FALSE) +
    geom_text(nudge_x = 0.02, size=8, show.legend = FALSE) +
    plot_theme() + xlim(0,0.4) +
    scale_colour_manual(values=gg_color_hue(4)[1:3])
  dev.off()

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  plot_theme <- function() {
    theme_grey() +
      theme(legend.position = "right",
            panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            strip.background = element_rect(fill = NA, colour = "black"),
            axis.text.x = element_text(angle = 90, vjust = 0.5,
                                       hjust = 1, size = 15),
            axis.text.y = element_text(size = 15),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20))
  }
  
}
