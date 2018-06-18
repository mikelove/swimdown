DIR='out'
OUT='quant'
INDEX='gencode.v28_salmon_0.10.0'
salmon='salmon-0.10.0_linux_x86_64/bin/salmon'
for ID in `seq 1 4`; do
    echo $ID;
    $salmon quant -p 6 -i $INDEX -l IU \
      --numBootstraps 30 \
      --gcBias \
      -o $OUT/salmon/$ID\_1 \
      -1 $DIR/out_$ID/sample_01_1_shuffled.fa.gz \
      -2 $DIR/out_$ID/sample_01_2_shuffled.fa.gz;
    $salmon quant -p 6 -i $INDEX -l IU \
      --numBootstraps 30 \
      --gcBias \
      -o $OUT/salmon/$ID\_2 \
      -1 $DIR/out_$ID/sample_02_1_shuffled.fa.gz \
      -2 $DIR/out_$ID/sample_02_2_shuffled.fa.gz;
done
