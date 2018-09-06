library(DRIMSeq)
library(DEXSeq)
library(iCOBRA)
n <- 3
sample <- paste0("s",rep(1:n,2),"_",rep(1:2,each=n))
samps <- data.frame(sample_id=sample, condition=factor(rep(1:2,each=n)))

# the Soneson et al (2016) counts for the Homo sapiens simulation have been deposited here:
# https://doi.org/10.5281/zenodo.1409201

library(readr)
cts <- matrix(nrow=145342,ncol=6)
for (i in 1:6) {
  x <- read_delim(paste0("isoform_prefilt_data/sim1hsapiens_original_kallisto/kallisto",i,".txt"), delim="\t", col_names=FALSE)
  names(x) <- c("transcript_id", "count")
  if (i == 1) {
    rownames(cts) <- x$transcript_id
  }
  cts[,i] <- x$count
}
colnames(cts) <- sample
truth <- read_delim("isoform_prefilt_data/E-MTAB-3766.additional/Hs_truth.txt", delim=" ")
counts <- data.frame(gene_id=truth[[1]], feature_id=truth[[2]], cts)
true.genes <- as.character(unique(truth$gene_id[truth$gene_ds_status == 1]))

#
d <- dmDSdata(counts=counts, samples=samps)
# filter05_1samp
d <- dmFilter(d,
              min_samps_feature_expr=0, min_feature_expr=0,
              min_samps_feature_prop=1, min_feature_prop=0.05,
              min_samps_gene_expr=0,  min_gene_expr=0)
# filter05_3samp
d <- dmFilter(d,
              min_samps_feature_expr=0, min_feature_expr=0,
              min_samps_feature_prop=n, min_feature_prop=0.05,
              min_samps_gene_expr=0,  min_gene_expr=0)
# filter10_3samp_cts
d <- dmFilter(d,
              min_samps_feature_expr=n, min_feature_expr=10,
              min_samps_feature_prop=n, min_feature_prop=0.10,
              min_samps_gene_expr=2*n,  min_gene_expr=10)
d
design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
system.time({
  dex.res <- dex(d)
})
system.time({
  drim.res <- drim(d)
})

#save(dex.res, drim.res, file="isoform_prefilt_results/filter05_1samp.rda")
#save(dex.res, drim.res, file="isoform_prefilt_results/filter05_3samp.rda")
#save(dex.res, drim.res, file="isoform_prefilt_results/filter10_3samp_cts.rda")

padj.gene <- data.frame(row.names=unique(counts(d)$gene_id))
load("isoform_prefilt_results/filter05_1samp.rda")
padj.gene$Prop05.any.DRIMSeq <- drim.res$adj_pvalue[match(rownames(padj.gene),drim.res$gene_id)]
padj.gene$Prop05.any.DEXSeq <- dex.res[rownames(padj.gene)]
load("isoform_prefilt_results/filter05_3samp.rda")
padj.gene$Prop05.half.DRIMSeq <- drim.res$adj_pvalue[match(rownames(padj.gene),drim.res$gene_id)]
padj.gene$Prop05.half.DEXSeq <- dex.res[rownames(padj.gene)]
load("isoform_prefilt_results/filter10_3samp_cts.rda")
padj.gene$Prop10.half.cts.DRIMSeq <- drim.res$adj_pvalue[match(rownames(padj.gene),drim.res$gene_id)]
padj.gene$Prop10.half.cts.DEXSeq <- dex.res[rownames(padj.gene)]

padj.gene <- padj.gene[apply(as.matrix(padj.gene), 1, function(x) !all(is.na(x))),]
truth <- data.frame(status=as.numeric(rownames(padj.gene) %in% true.genes),
                    row.names=rownames(padj.gene))
cd.gene <- COBRAData(padj=padj.gene, truth=truth)

pdf("isoform_prefilt_results/isoform_prefilt.pdf", height=6)
print(myicobra(cd.gene, addit="Soneson (2016)"))
dev.off()

# functions

drim <- function(d) {
  set.seed(1)
  suppressMessages({
    d <- dmPrecision(d, design=design_full)
  })
  d <- dmFit(d, design=design_full)
  d <- dmTest(d, coef="condition2")
  res <- DRIMSeq::results(d)
  res$pvalue <- ifelse(is.na(res$pvalue), 1, res$pvalue)
  res
}

dex <- function(d) {
  sample.data <- DRIMSeq::samples(d)
  count.data <- round(counts(d)[,-c(1:2)])
  suppressMessages({dxd <- DEXSeqDataSet(countData=count.data,
                       sampleData=sample.data,
                       design=~sample + exon + condition:exon,
                       featureID=counts(d)$feature_id,
                       groupID=counts(d)$gene_id)})
  set.seed(1)
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, quiet=TRUE)
  dxd <- testForDEU(dxd, reducedModel=~sample + exon)
  dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
  qval <- perGeneQValue(dxr)
  qval
}

myicobra <- function(cd, lvl="Gene", addit) {
  cp <- calculate_performance(cd,
                              binary_truth="status",
                              aspect=c("fdrtpr","fdrtprcurve"),
                              thrs=c(.01,.05,.1))
  cobraplot <- prepare_data_for_plot(cp, colorscheme="Paired")
  plot_fdrtprcurve(cobraplot, plottype="points",
                   xaxisrange=c(0,.3),
                   yaxisrange=c(.5, .8),
                   stripsize=0,
                   title=paste0(lvl,"-level screening, n=",n," vs ",n,", ",addit))
}
