
sample=$1
bam_file=$2
seed=$3

source activate DSN-seq

for frac in {1,2,3,4,5,6,7,8,9}
do
samtools view -@ 2 -b -s ${seed}.${frac} ${bam_file} -o ./sample_bam/${sample}.${seed}_${frac}.filtered.sorted.bam
samtools index -@ 20 ./sample_bam/${sample}.${seed}_${frac}.filtered.sorted.bam
done
