mkdir sample_bam
source activate DSN-seq

# Osa_DSN_merge
sample=Osa_DSN_merge
bam=/home/hj-z2/projects/DSN-seq/results/P1_call_editing/01.mapping/Osa_DSN_merge.filtered.sorted.bam

for seed in {11,22,33}
do
bash s1_bam_sample_base.sh $sample $bam $seed &
done



# Osa_Ribo-off_merge
sample=Osa_Ribo-off_merge
bam=/home/hj-z2/projects/DSN-seq/results/P1_call_editing/01.mapping/Osa_Ribo-off_merge.filtered.sorted.bam

for seed in {11,22,33}
do
bash s1_bam_sample_base.sh $sample $bam $seed &
done



# Ath_DSN_merge
sample=Ath_DSN_merge
bam=/home/hj-z2/projects/DSN-seq/results/P1_call_editing/01.mapping/Ath_DSN_merge.filtered.sorted.bam

for seed in {11,22,33}
do
bash s1_bam_sample_base.sh $sample $bam $seed &
done



