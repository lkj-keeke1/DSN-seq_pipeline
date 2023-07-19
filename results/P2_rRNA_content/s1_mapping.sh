Osa_ref=/home/hj-z2/projects/DSN-seq/data/ref/Osa/Osa_ncbi_rRNA.fasta
Ath_ref=/home/hj-z2/projects/DSN-seq/data/ref/Ath/Ath_ncbi_rRNA.fasta

cat ../../data/seq_data/RNA_editing/UID/sample.UID.list | grep "Osa" | while read sample R1 R2
do
bash bwa_mapping.sh $sample $Osa_ref $R1 $R2 &
done


cat ../../data/seq_data/RNA_editing/UID/sample.UID.list | grep "Ath" | while read sample R1 R2
do
bash bwa_mapping.sh $sample $Ath_ref $R1 $R2 &
done

