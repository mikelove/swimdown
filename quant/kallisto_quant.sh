DIR='out'
OUT='quant'
INDEX='gencode.v28_kallisto_0.44.0'
kallisto='kallisto_linux-v0.44.0/kallisto'
for ID in `seq 1 4`; do
    echo $ID;
    $kallisto quant -t 6 -b 30 -i $INDEX --bias \
      -o $OUT/kallisto/$ID\_1 \
      $DIR/out_$ID/sample_01_1_shuffled.fa.gz \
      $DIR/out_$ID/sample_01_2_shuffled.fa.gz;
    $kallisto quant -t 6 -b 30 -i $INDEX --bias \
      -o $OUT/kallisto/$ID\_2 \
      $DIR/out_$ID/sample_02_1_shuffled.fa.gz \
      $DIR/out_$ID/sample_02_2_shuffled.fa.gz;
done
