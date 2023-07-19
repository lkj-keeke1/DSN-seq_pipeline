source activate DSN-seq

wkdir=$1 # absolute path

sample=$2 # sample name, the prefix of output files



## Call editing parameters
ref_genome=$3              # multi-FASTA file, absolute path
length_table=$4                  # fasta sequences length table, absolute path
selction_script=~/projects/DSN-seq/src/vcf_modify.R # vcf file filtering script, absolute path

# samtools parameters
max_depth=100000000 # max depth for mpileup files

# filtering parameters
f_p_value=0.01  # filtering standard for p_value
f_freq=10       # filtering standard for editing frequency/extent (%)
f_depth=8       # filtering standard for sequencing depth


## use Varscan2 to call variants

cd $wkdir/02.call_editing 
## Varscan2 calls variants without filtering
samtools mpileup -d ${max_depth} -f ${ref_genome} ${wkdir}/01.mapping/${sample}.filtered.sorted.bam > ${sample}.mpileup
awk '{print $1" "$2"\t"$4}' ${sample}.mpileup > ${sample}.depth 
java -jar ~/softwares/varscan-2.4.5/VarScan.v2.4.5.jar mpileup2snp ${sample}.mpileup  --output-vcf 1 --p-value 1 --min-var-freq 0 --strand-filter 0  > ${sample}.vcf
 # We adjust the parameters, so that the filtering is done in next step. 

## Obtain confident RNA editing sites
Rscript ${selction_script} ${sample}.vcf ${length_table} ${sample} ${f_p_value} ${f_freq} ${f_depth} &
