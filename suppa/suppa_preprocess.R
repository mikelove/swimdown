load("txi.rda")
x <- txi.list[["salmon"]]$abundance
# eliminate 1e-9's
x[x < 0.01] <- 0

n <- 12
write.table(x[,1:n], file=paste0("suppa/group1_",n,".tpm"), quote=FALSE, sep="\t")
write.table(x[,12 + 1:n], file=paste0("suppa/group2_",n,".tpm"), quote=FALSE, sep="\t")
