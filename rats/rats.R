library(DRIMSeq)
n <- 12
sample <- paste0("s",rep(1:n,2),"_",rep(1:2,each=n))
samps <- data.frame(sample_id=sample, condition=factor(rep(1:2,each=n)))
load("txi_scaled.rda")
# these are scaledTPM counts from abundance,
# this is recommended in the RATs vignette:
# https://github.com/bartongroup/RATS/blob/master/vignettes/input.Rmd#read-counts-tpms-etc
cts <- txi.list[["salmon"]]$counts
colnames(cts) <- sample
load("data/simulate.rda")

print("all tximport counts in the TXDF")
all(rownames(cts) %in% txdf$TXNAME)
txdf <- txdf[match(rownames(cts),txdf$TXNAME),]

print("counts and TXDF are aligned")
all(rownames(cts) == txdf$TXNAME)
counts <- data.frame(gene_id=txdf$GENEID,
                     feature_id=txdf$TXNAME,
                     cts)

tx2gene <- counts[,2:1]
names(tx2gene) <- c("target_id","parent_id")

library(rats)

for (n.sub in c(2,3,6,9,12)) {
  print(n.sub)
  idx <- c(1:n.sub, n + 1:n.sub)

  # provide DRIMSeq filtered counts to RATs:
  # using the same minimal counts and proportion filters
  # as used with the other methods
  d <- dmDSdata(counts=counts[,c(1:2,2+idx)], samples=samps[idx,])
  d <- dmFilter(d,
                min_samps_feature_expr=n.sub, min_feature_expr=10,
                min_samps_feature_prop=n.sub, min_feature_prop=0.1,
                min_samps_gene_expr=2*n.sub, min_gene_expr=10)
  d

  countA <- data.table(counts(d)[,c("feature_id",paste0("s",1:n.sub,"_1"))])
  countB <- data.table(counts(d)[,c("feature_id",paste0("s",1:n.sub,"_2"))])
  
  # 7 minutes for counts for 6 vs 6
  st_count <- system.time({
    rats_count <- call_DTU(annot=tx2gene, count_data_A=countA, count_data_B=countB, qboot=FALSE)
  })
  
  save(rats_count, st_count, file=paste0("rats/rats_",n.sub,".rda"))
}
