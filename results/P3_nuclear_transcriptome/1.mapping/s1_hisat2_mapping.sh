
ref=~/projects/DSN-seq/data/ref/Osa/Osa_genome_chr

source activate DSN-seq

cat ~/projects/DSN-seq/data/seq_data/RNA_editing/UID/sample.UID.list | grep "Osa" | while read sample  R1 R2
do
bash hisat2_mapping.sh $sample $ref $R1 $R2 &
done


