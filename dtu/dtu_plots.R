rgb2 <- function(x,y,z) rgb(x,y,z,maxColorValue=255)
cols <- c("DEXSeq"="black",
          "DRIMSeq"=rgb2(230,159,0),
          "DRIMSeq.filt"=rgb2(86,180,233),
          "DRIMSeq-filt"=rgb2(86,180,233),
          "RATs"=rgb2(0,158,115),
          "SUPPA2"=rgb2(0,114,178))

ofdr <- function(padj, alpha=.05) {
  all.not.dtu <- sapply(with(padj, split(dtu, geneID)), function(x) ! any(x))
  sig.not.dtu <- sapply(with(padj, split(transcript < alpha & ! dtu, geneID)), any)
  mean(all.not.dtu | sig.not.dtu)
}

myicobra <- function(cd, lvl="Gene", label=FALSE) {
  cp <- calculate_performance(cd,
                              binary_truth="status",
                              aspect=c("fdrtpr","fdrtprcurve"),
                              thrs=c(.01,.05,.1))
  cobraplot <- prepare_data_for_plot(cp, colorscheme=cols)
  plt <- plot_fdrtprcurve(cobraplot, plottype="points",
                          xaxisrange=c(0,max(max(fdrtpr(cp)$FDR),.21)),
                          yaxisrange=c(0,1),
                          #yaxisrange=c(0,max(fdrtpr(cp)$TPR)),
                          stripsize=0,
                          title=paste0(lvl,"-level screening, n=",n.sub," vs ",n.sub))
  if (label)
    plt <- plt + geom_text(aes(label=method,color=method,vjust=-1))
  plt
}

tpBreakdownGene <- function(padj.gene, alpha) {
  padj.gene <- padj.gene[apply(padj.gene, 1, function(x) !all(is.na(x))),]
  padj.gene$status <- rownames(padj.gene) %in% full.dtu.genes
  padj.gene$type <- ifelse(rownames(padj.gene) %in% dge.genes,
                           "dge",
                    ifelse(rownames(padj.gene) %in% dte.genes,
                           "dte",
                    ifelse(rownames(padj.gene) %in% dtu.genes,
                           "dtu", "null")))
  ref.tab <- table(padj.gene$type[padj.gene$status])
  meths <- setdiff(names(padj.gene), c("status","type"))
  tps <- sapply(meths, function(m) {
    table(TP=padj.gene[[m]] < alpha & padj.gene$status, padj.gene$type)[2,2:3]
  })
  cols <- c("violet","goldenrod")
  main <- paste0("Gene-level, n=", n.sub, " vs ", n.sub)
  par(mfrow=c(1,1))
  barplot(cbind(reference=ref.tab, tps), col=cols, main=main, ylab="True positives")
  legend("topright", rev(c("DTE","DTU")), fill=rev(cols), bg="white", title="Proportion by gene type")
}

fdrBreakdownGene <- function(padj.gene, alpha) {
  padj.gene <- padj.gene[apply(padj.gene, 1, function(x) !all(is.na(x))),]
  padj.gene$status <- rownames(padj.gene) %in% full.dtu.genes
  padj.gene$type <- ifelse(rownames(padj.gene) %in% dge.genes,
                           "dge",
                    ifelse(rownames(padj.gene) %in% dte.genes,
                           "dte",
                    ifelse(rownames(padj.gene) %in% dtu.genes,
                           "dtu", "null")))
  ref.props <- prop.table(table(padj.gene$type[!padj.gene$status]))
  meths <- setdiff(names(padj.gene), c("status","type"))
  fdrs <- sapply(meths, function(m) {
    tot <- sum(padj.gene[[m]] < alpha, na.rm=TRUE)
    table(FP=padj.gene[[m]] < alpha & !padj.gene$status, padj.gene$type)[2,c(1,2,4)]/tot
  })
  cols <- c("seagreen","violet","skyblue")
  main <- paste0("Gene-level, n=", n.sub, " vs ", n.sub)
  par(mfrow=c(2,1))
  ylim <- c(0, max(.4, max(fdrs)))
  barplot(fdrs, col=cols, main=main, ylab="Total height = FDR", ylim=ylim)
  abline(h=alpha, lty=2)
  legend("topright", rev(c("DGE","DTE","null")), fill=rev(cols), bg="white", title="Proportion by gene type")
  barplot(cbind(reference=ref.props, t(t(fdrs)/colSums(fdrs))),
          col=cols, main=main, ylab="Proportion of FP")
}

fdrBreakdownTxp <- function(padj, alpha) {
  padj <- padj[apply(padj, 1, function(x) !all(is.na(x))),]
  padj$status <- rownames(padj) %in% dtu.txps
  padj$gene <- txdf$GENEID[match(rownames(padj), txdf$TXNAME)]
  padj$type <- ifelse(padj$gene %in% dge.genes,
                      "dge",
               ifelse(padj$gene %in% dte.genes,
                      "dte",
               ifelse(padj$gene %in% dtu.genes,
                      "dtu", "null")))
  ref.props <- prop.table(table(padj$type[!padj$status]))
  meths <- setdiff(names(padj), c("status","gene","type"))
  fdrs <- sapply(meths, function(m) {
    tot <- sum(padj[[m]] < alpha, na.rm=TRUE)
    table(FP=padj[[m]] < alpha & !padj$status, padj$type)[2,]/tot
  })
  cols <- c("seagreen","violet","goldenrod","skyblue")
  main <- paste0("Transcript-level, n=", n.sub, " vs ", n.sub)
  par(mfrow=c(2,1))
  ylim <- c(0, max(.4, max(fdrs)))
  barplot(fdrs, col=cols, main=main, ylab="Total height = FDR", ylim=ylim)
  abline(h=alpha, lty=2)
  legend("topright", rev(c("DGE","DTE","DTU","null")), fill=rev(cols), bg="white", title="Proportion by gene type")
  barplot(cbind(reference=ref.props, t(t(fdrs)/colSums(fdrs))),
          col=cols, main=main, ylab="Proportion of FP")
}

plotOFDR <- function(with2=TRUE) {
  if (with2) {
    x <- read.table(file.path(dir,"OFDR_2.tsv"), header=TRUE)
    for (i in c(3,6,9,12)) {
      x <- rbind(x, read.table(paste0(dir,"/OFDR_",i,".tsv"), header=TRUE))
    }
  } else {
    x <- read.table(file.path(dir,"OFDR_3.tsv"), header=TRUE)
    for (i in c(6,9,12)) {
      x <- rbind(x, read.table(paste0(dir,"/OFDR_",i,".tsv"), header=TRUE))
    }
  }
  if (dir == "res_scaled") {
    if (with2) {
      x$n <- factor(rep(c(2,1:4*3), each=5), levels=c(2,1:4*3))
    } else {
      x$n <- factor(rep(1:4*3, each=5), levels=1:4*3)
    }
  } else {
    x$n <- factor(rep(1:4*3, each=3), levels=1:4*3)
  }
  nudgex <- rep(.02, nrow(x))
  nudgey <- rep(0, nrow(x))
  if (dir == "res_scaled") {
    if (with2) {
      nudgex[x$method=="RATS"] <- c(.02, .02, .02, .01, 0)
      nudgex[x$method=="SUPPA2"] <- c(.02, .02, .02, .01, 0)
    } else {
      nudgex[x$method=="RATS"] <- c(.02, .02, .01, 0)
      nudgex[x$method=="SUPPA2"] <- c(.02, .02, .01, 0)
    }
    nudgey[x$method=="RATS"] <- -.02
    nudgey[x$method=="SUPPA2"] <- -.03
  }
  x$control <- ifelse(x$OFDR < .05, 19, 1)
  
  ggplot(x, aes(OFDR, TPR, color=method, label=n)) +
    ylim(0,1) + ylab("TPR (transcripts)") +
    geom_vline(xintercept=0.05, color="gray", linetype=2) + 
    geom_point(aes(shape=control), size=3) + geom_path(show.legend=FALSE) +
    geom_text(nudge_x = nudgex, nudge_y = nudgey, size=4, show.legend = FALSE) +
    plot_theme() + xlim(0,0.4) + scale_shape_identity() +
    scale_color_manual(values=cols)
  
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
