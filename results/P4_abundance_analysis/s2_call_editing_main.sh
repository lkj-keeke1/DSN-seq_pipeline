mkdir 1.mpileup 2.vcf 3.filter


# Osa_DSN_merge
sample=Osa_DSN_merge
ref_genome=~/projects/DSN-seq/data/ref/Osa/Osa_organelle_multi-FASTA.fa
length_table=~/projects/DSN-seq/data/ref/Osa/Osa_organelle_gene_length.table

for seed in {11,22,33}
do
bash s2_call_editing_base.sh $sample $ref_genome $length_table $seed &
done


# Osa_Ribo-off_merge
sample=Osa_Ribo-off_merge
ref_genome=~/projects/DSN-seq/data/ref/Osa/Osa_organelle_multi-FASTA.fa
length_table=~/projects/DSN-seq/data/ref/Osa/Osa_organelle_gene_length.table

for seed in {11,22,33}
do
bash s2_call_editing_base.sh $sample $ref_genome $length_table $seed &
done


# Ath_DSN_merge
sample=Ath_DSN_merge
ref_genome=~/projects/DSN-seq/data/ref/Ath/Ath_organelle_multi-FASTA.fa
length_table=~/projects/DSN-seq/data/ref/Ath/Ath_organelle_gene_length.table

for seed in {11,22,33}
do
bash s2_call_editing_base.sh $sample $ref_genome $length_table $seed &
done


