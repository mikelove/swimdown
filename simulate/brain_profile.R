samps <- read.delim("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
idx <- grepl("Frontal Cortex", samps$SMTSD)
table(idx)
head(samps$SAMPID[idx])
tpm <- read.delim("GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz", skip=2, nrows=1, header=FALSE)
tpm_samps <- tpm[1,]
table(tpm_samps %in% samps$SAMPID[idx])
head(which(tpm_samps %in% samps$SAMPID[idx]), 10)
# 206 273 323 344 364 411 512 548 818 865

# zcat GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz | cut -f 1,2,206,273,323,344,364,411,512,548,818,865 > frontal_cortex.tsv
