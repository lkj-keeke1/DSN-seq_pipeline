source activate DSN-seq

# fastqc output dir
mkdir bef_qc aft_qc

# fastq file before qc
ls /home/hj-z2/projects/DSN-seq/data/seq_data/RNA_editing/raw/*.gz |xargs fastqc -t 30 -o bef_qc
# fastq file after qc, clean and UMI-dedup
ls /home/hj-z2/projects/DSN-seq/data/seq_data/RNA_editing/UID/*.gz |xargs fastqc -t 30 -o aft_qc

# multiqc summarise 
cd bef_qc
multiqc . &

cd ../aft_qc
multiqc . &
