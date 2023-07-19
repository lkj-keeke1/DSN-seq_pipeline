##### seqkit stat #####

 # counts

# read number
 # raw 
 # clean
# base number
 # raw 
 # clean

source activate DSN-seq

 ## raw fq.gz
seqkit stat -a -b -T ~/projects/DSN-seq/data/seq_data/RNA_editing/raw/*.gz -o RNA_editing_raw_seqkit_stat.table &
 ## UID fq.gz
seqkit stat -a -b -T ~/projects/DSN-seq/data/seq_data/RNA_editing/UID/*.gz -o RNA_editing_UID_seqkit_stat.table &

