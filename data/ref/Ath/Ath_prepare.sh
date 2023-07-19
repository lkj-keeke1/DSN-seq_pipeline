source activate DSN-seq

# organelle
hisat2-build Ath_organelle_multi-FASTA.fa Ath_organelle_multi-FASTA >hisat2-build.log 2>&1
samtools faidx Ath_organelle_multi-FASTA.fa
cut -f1,2 Ath_organelle_multi-FASTA.fa.fai > Ath_organelle_gene_length.table

# rRNA
cp /data/1.public_data/genomes/NCBI/Arabidopsis_thaliana_GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_rna_from_genomic.fna.gz .
gunzip -d GCF_000001735.4_TAIR10.1_rna_from_genomic.fna.gz 
grep 'gbkey=rRNA' GCF_000001735.4_TAIR10.1_rna_from_genomic.fna | sed 's/^>\([^ ]*\).*/\1/' > Ath_ncbi_rRNA.idlist
seqtk subseq GCF_000001735.4_TAIR10.1_rna_from_genomic.fna Ath_ncbi_rRNA.idlist > Ath_ncbi_rRNA.fasta
bwa index Ath_ncbi_rRNA.fasta 
