
source activate DSN-seq
cat /home/hj-z2/projects/DSN-seq/data/seq_data/RNA_editing/UID/sample.UID.list | while read sample R1 R2
do 
cutadapt -j 2 -u 6 -o ./${sample}.tail_cut.R1.fq.gz  ${R1} 1>/dev/null 2>error.log &
cutadapt -j 2 -u 6 -o ./${sample}.tail_cut.R2.fq.gz  ${R2} 1>/dev/null 2>error.log &
done


