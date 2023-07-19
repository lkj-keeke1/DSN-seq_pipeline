sample=$1
ref_genome=$2
length_table=$3
seed=$4

source activate DSN-seq

selction_script=~/projects/DSN-seq/src/vcf_modify.R
 
#filtering parameters
f_p_value=0.01  # filtering standard for p_value
f_freq=10       # filtering standard for editing frequency/extent (%)
f_depth=8       # filtering standard for sequencing depth

for frac in {1,2,3,4,5,6,7,8,9}
do
samtools mpileup -f ${ref_genome} ./sample_bam/${sample}.${seed}_${frac}.filtered.sorted.bam > ./1.mpileup/${sample}.${seed}_${frac}.mpileup 
java -jar ~/softwares/varscan-2.4.5/VarScan.v2.4.5.jar mpileup2snp ./1.mpileup/${sample}.${seed}_${frac}.mpileup --output-vcf 1 --p-value 1 --min-var-freq 0 --strand-filter 0 > ./2.vcf/${sample}.${seed}_${frac}.vcf 
Rscript ${selction_script} ./2.vcf/${sample}.${seed}_${frac}.vcf ${length_table} ./3.filter/${sample}.${seed}_${frac} ${f_p_value} ${f_freq} ${f_depth} &
done
