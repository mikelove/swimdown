args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])
library(polyester)
load("data/simulate.rda")
se <- simulate_experiment(fasta="data/transcripts.fa",
                          outdir=paste0("out/out_",i),
                          num_reps=c(1,1),
                          reads_per_transcript=sim_counts,
                          size=1/disps,
                          fold_changes=fold_changes,
                          fraglen=200,
                          fragsd=25,
                          frag_GC_bias=frag_GC_bias[,c(i,i)],
                          readlen=100,
                          seed=i)
