library(DRIMSeq)
n <- 12
sample <- paste0("s",rep(1:n,2),"_",rep(1:2,each=n))
samps <- data.frame(sample_id=sample, condition=factor(rep(1:2,each=n)))
load("txi_scaled.rda")
cts <- txi.list[["salmon"]]$counts
colnames(cts) <- sample
load("data/simulate.rda")

## n <- 12
## sample <- paste0(rep(1:n,2),"_",rep(1:2,each=n))
## dir <- "quant"
## library(tximport)
## files <- file.path(dir,"salmon",sample,"quant.sf")
## names(files) <- sample
## txi.boot <- tximport(files, type="salmon", txOut=TRUE)
## save(txi.boot, file="txi_boot.rda")
load("txi_boot.rda")

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

  # the DRIMSeq object 'd' is used to pick transcripts from
  # the bootstrap count matrices
  m <- match(counts(d)$feature_id, rownames(txi.boot$counts))

  boot <- list() # list of lists of bootstrap matrices (each one is txps x bootstraps)
  for (j in 1:2) {
    boot[[j]] <- lapply(1:n.sub, function(i) {
      # pull out bootstrap estimated counts
      # j is the group (1 or 2)
      # i is the sample within the group
      mat <- txi.boot$infReps[[paste0(i,"_",j)]]
      # the following lines generates scaledTPM counts from abundance
      # and is recommended from the RATs vignette:
      # https://github.com/bartongroup/RATS/blob/master/vignettes/input.Rmd#read-counts-tpms-etc
      # divide out effective length to get something propto TPM
      abund <- mat / txi.boot$length[,paste0(i,"_",j)]
      # scale abundances up to total mapped counts ('cts' is over all transcripts)
      scaledTPM <- t(t(abund) / colSums(abund) * colSums(cts)[paste0("s",i,"_",j)])
      # subset to the filtered transcripts/genes
      scaledTPM <- scaledTPM[m,]
      data.table(feature_id=counts(d)$feature_id, scaledTPM)
    })
  }

  st_boot <- system.time({
    rats_boot <- call_DTU(annot=tx2gene, boot_data_A=boot[[1]], boot_data_B=boot[[2]], qboot=TRUE)
  })
  
  save(rats_boot, st_boot, file=paste0("rats/rats_boot_",n.sub,".rda"))
}
