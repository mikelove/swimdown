library(sleuth)
n <- 12
sample <- paste0(rep(1:n,2),"_",rep(1:2,each=n))
files <- paste0("quant/kallisto/",sample)
samps <- data.frame(sample,
                    condition=factor(rep(1:2,each=n)),
                    path=files, stringsAsFactors=FALSE)

suppressPackageStartupMessages(library(GenomicFeatures))
txdb.filename <- "gencode.v28.annotation.sqlite"
txdb <- loadDb(txdb.filename)
txdf <- select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tx2gene <- txdf[,2:1]
colnames(tx2gene) <- c("target_id","gene_id")

# here comes 'n.sub'
n.sub <- 12
set.seed(1)
idx <- c(1:n.sub, n + 1:n.sub)

st <- system.time({
  so <- sleuth_prep(samps[idx,], ~condition, num_cores=1)
  so <- sleuth_fit(so)
  so <- sleuth_fit(so, ~1, 'reduced')
  so <- sleuth_lrt(so, 'reduced', 'full')
  sleuth.res <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
})

st.gene <- system.time({
  so <- sleuth_prep(samps[idx,], ~condition,
                    target_mapping=tx2gene,
                    aggregation_column="gene_id",
                    num_cores=1)
  so <- sleuth_fit(so)
  so <- sleuth_fit(so, ~1, 'reduced')
  so <- sleuth_lrt(so, 'reduced', 'full')
  sleuth.res.gene <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
})

save(sleuth.res, sleuth.res.gene, st, st.gene,
     file=paste0("sleuth/sleuth_",n.sub,".rda"))
