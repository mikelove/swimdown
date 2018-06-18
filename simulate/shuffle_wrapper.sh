CMD='shuffle.sh'
DIR='out'
for i in `seq 1 4`; 
do
    $CMD -l $DIR/out_$i/sample_01_1.fasta -r $DIR/out_$i/sample_01_2.fasta;
    $CMD -l $DIR/out_$i/sample_02_1.fasta -r $DIR/out_$i/sample_02_2.fasta;
    rm $DIR/out_$i/*.fasta
    pigz -1 -p 4 $DIR/out_$i/sample_01_1_shuffled.fa;
    pigz -1 -p 4 $DIR/out_$i/sample_01_2_shuffled.fa;
    pigz -1 -p 4 $DIR/out_$i/sample_02_1_shuffled.fa;
    pigz -1 -p 4 $DIR/out_$i/sample_02_2_shuffled.fa;
done
