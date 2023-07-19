wkdir=/home/hj-z2/projects/DSN-seq/results/P1_call_editing
Ath_ref=/home/hj-z2/projects/DSN-seq/data/ref/Ath/Ath_organelle_multi-FASTA
Ath_len=/home/hj-z2/projects/DSN-seq/data/ref/Ath/Ath_organelle_gene_length.table
Ath_ge=/home/hj-z2/projects/DSN-seq/data/ref/Ath/Ath_organelle_multi-FASTA.fa
Osa_ref=/home/hj-z2/projects/DSN-seq/data/ref/Osa/Osa_organelle_multi-FASTA
Osa_len=/home/hj-z2/projects/DSN-seq/data/ref/Osa/Osa_organelle_gene_length.table
Osa_ge=/home/hj-z2/projects/DSN-seq/data/ref/Osa/Osa_organelle_multi-FASTA.fa

ls Osa*.gz | awk -F'.' '{print $1}' | uniq | while read id ;do echo -e "$wkdir\t$id\t`pwd`/$id.tail_cut.R1.fq.gz\t`pwd`/$id.tail_cut.R2.fq.gz\t$Osa_ref\t$Osa_ge\t$Osa_len" >> sample.tail_cut.list;done
ls Ath*.gz | awk -F'.' '{print $1}' | uniq | while read id ;do echo -e "$wkdir\t$id\t`pwd`/$id.tail_cut.R1.fq.gz\t`pwd`/$id.tail_cut.R2.fq.gz\t$Ath_ref\t$Ath_ge\t$Ath_len" >> sample.tail_cut.list;done
