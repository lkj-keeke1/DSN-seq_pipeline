source activate DSN-seq
# fastqc before qc
ls /data/0.user_data/hj-z2/DSN-seq/GRP23/SRR19779049_* |xargs fastqc -o . -t 2  2>/dev/null &

# qc
fastp -h SRR19779049.html -j SRR19779049.json -i /data/0.user_data/hj-z2/DSN-seq/GRP23/SRR19779049_1.fastq.gz -I /data/0.user_data/hj-z2/DSN-seq/GRP23/SRR19779049_2.fastq.gz  -o SRR19779049_R1.fq.gz -O SRR19779049_R2.fq.gz -w 16 -f 8 -F 8

# fastqc after qc
ls SRR19779049_R*.gz |xargs fastqc -o . -t 2 2>/dev/null &

