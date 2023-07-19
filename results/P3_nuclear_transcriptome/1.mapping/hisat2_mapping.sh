sample=$1
ref=$2
R1=$3
R2=$4

source activate DSN-seq

hisat2 --new-summary -p 10 -x ${ref} -1 ${R1} -2 ${R2} --rna-strandness FR 2>${sample}.hisat2.log |  samtools view -bS - | samtools sort - > ${sample}.sorted.bam
samtools index ${sample}.sorted.bam 
samtools flagstat ${sample}.sorted.bam > ${sample}.flagstat &

