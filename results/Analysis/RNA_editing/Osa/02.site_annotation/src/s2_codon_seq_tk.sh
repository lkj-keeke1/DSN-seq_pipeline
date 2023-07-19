genome=~/projects/DSN-seq/data/ref/Osa/Osa_organelle_multi-FASTA.fa

seqkit subseq ${genome} --bed tmp_codon_tk.bed > tmp_codonseq.fasta
sed 's/^>//' tmp_codonseq.fasta |sed ':a;N;s/\n/\t/g' | sed 's/:. /\t/' > tmp_pre_codon.table
