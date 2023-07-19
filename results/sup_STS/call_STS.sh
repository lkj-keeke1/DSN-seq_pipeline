source activate DSN-seq

sample=Ath_STS # sample name, the prefix of output files

wkdir=/home/hj-z2/projects/DSN-seq/results/sup_STS # absolute path

## mapping parameters
clean_R1=/home/hj-z2/projects/DSN-seq/results/sup_STS/00.qc/SRR19779049_R1.fq.gz  # clean R1 fastq file, absolute path
clean_R2=/home/hj-z2/projects/DSN-seq/results/sup_STS/00.qc/SRR19779049_R2.fq.gz  # clean R2 fastq file, absolute path
hisat_ref=/home/hj-z2/projects/DSN-seq/data/ref/Ath/Ath_organelle_multi-FASTA      # hisat2 reference, absolute path

hisat2_thread=20   # threads for hisat2 mapping step
samtools_thread=10 # threads for samtools view and sor


## Call editing parameters
ref_genome=/home/hj-z2/projects/DSN-seq/data/ref/Ath/Ath_organelle_multi-FASTA.fa              # multi-FASTA file, absolute path
length_table=/home/hj-z2/projects/DSN-seq/data/ref/Ath/Ath_organelle_gene_length.table                  # fasta sequences length table, absolute path
selction_script=~/projects/DSN-seq/src/vcf_modify.R # vcf file filtering script, absolute path

# samtools parameters
max_depth=100000 # max depth for mpileup files

# filtering parameters
f_p_value=0.01  # filtering standard for p_value
f_freq=10       # filtering standard for editing frequency/extent (%)
f_depth=8       # filtering standard for sequencing depth

cd $wkdir

########## s0:mkdir ##########

if [ -d "01.mapping" ];then
	echo "01.mapping already exists"
	else
	mkdir 01.mapping
fi

if [ -d "02.call_editing" ];then
        echo "02.call_editing already exists"
        else
        mkdir 02.call_editing
fi



########## s1:mapping to reference ##########

## Hisat2 aligns the clean reads (cut first 6 bp) to reference multi-FASTA genome.

cd $wkdir/01.mapping
## Hisat2 mapping
hisat2 --new-summary -p ${hisat2_thread} --mp 2,2 -x ${hisat_ref} -1 ${clean_R1} -2 ${clean_R2} 2>${sample}.hisat2.log | samtools view -@ ${samtools_thread} -bS - | samtools sort -@ ${samtools_thread} - > ${sample}.sorted.bam

## Keep only the properly paired reads
samtools view -@ ${samtools_thread} -f 2 -F 256 ${sample}.sorted.bam -b >${sample}.filtered.sorted.bam 

## samtools index
samtools index ${sample}.sorted.bam
samtools index ${sample}.filtered.sorted.bam

## mapping statistics
samtools flagstat ${sample}.sorted.bam > ${sample}.flagstat & 
samtools flagstat ${sample}.filtered.sorted.bam > ${sample}.filtered.flagstat &
samtools idxstats ${sample}.filtered.sorted.bam > ${sample}.filtered.idxstats &



########## s2: call editing sites ##########

## use Varscan2 to call variants
 
cd ${wkdir}/02.call_editing
## Varscan2 calls variants without filtering
samtools mpileup -d ${max_depth} -f ${ref_genome} ${wkdir}/01.mapping/${sample}.filtered.sorted.bam > ${sample}.mpileup
awk '{print $1" "$2"\t"$4}' ${sample}.mpileup > ${sample}.depth 
java -jar ~/softwares/varscan-2.4.5/VarScan.v2.4.5.jar mpileup2snp ${sample}.mpileup  --output-vcf 1 --p-value 1 --min-var-freq 0 --strand-filter 0  > ${sample}.vcf
 # We adjust the parameters, so that the filtering is done in next step. 

## Obtain confident RNA editing sites
Rscript ${selction_script} ${sample}.vcf ${length_table} ${sample} ${f_p_value} ${f_freq} ${f_depth} &
