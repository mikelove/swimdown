suppa <- read.delim("suppa/diff_empirical_12.dpsi")
names(suppa) <- c("txp.gene","dpsi","pval")
suppa$gene <- sub(";.*", "", suppa$txp.gene)
suppa$txp <- sub(".*;", "", suppa$txp.gene)
suppa <- suppa[!is.nan(suppa$dpsi),]
suppa$padj <- p.adjust(suppa$pval, method="BH")

load("data/simulate.rda")

# SUPPA

txdf <- txdf[match(rownames(tpms), txdf$TXNAME),]
txdf$dtu.genes <- iso.dtu | iso.dte & !iso.dte.only
full.dtu.genes <- unique(txdf$GENEID[txdf$dtu.genes])

txp.exprs <- rowSums(tpms) > 0
txdf$full.dtu <- iso.dtu | (txdf$GENEID %in% dte.genes & txp.exprs)
dtu.txps <- txdf$TXNAME[txdf$full.dtu]

suppa$dtu <- suppa$txp %in% dtu.txps
with(suppa, table(sig=padj < .1, dtu))
with(suppa, table(sig=padj < .05, dtu))
with(suppa, table(sig=padj < .01, dtu))

library(DEXSeq)
suppa.dxr <- as(DataFrame(groupID=suppa$gene,
                          pvalue=suppa$pval,
                          padj=rep(1, nrow(suppa))), "DEXSeqResults")
qval <- perGeneQValue(suppa.dxr)
suppa.g <- data.frame(gene=names(qval), qval=qval)
suppa.g$dtu <- suppa.g$gene %in% full.dtu.genes
with(suppa.g, table(sig=qval < .1, dtu))
with(suppa.g, table(sig=qval < .05, dtu))
with(suppa.g, table(sig=qval < .01, dtu))
