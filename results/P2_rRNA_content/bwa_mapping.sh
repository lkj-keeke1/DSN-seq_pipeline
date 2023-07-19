sample=$1
ref=$2
R1=$3
R2=$4

source activate DSN-seq
bwa mem -t 4 ${ref} ${R1} ${R2}  2>${sample}.bwamem.log | samtools view -bS - | samtools sort - > ${sample}.sorted.bam
samtools flagstat ${sample}.sorted.bam > ${sample}.flagstat &
