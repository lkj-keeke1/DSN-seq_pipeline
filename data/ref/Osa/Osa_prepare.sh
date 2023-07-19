source activate DSN-seq

# organalle genome
hisat2-build Osa_organelle_multi-FASTA.fa Osa_organelle_multi-FASTA >hisat2-build.log 2>&1
samtools faidx Osa_organelle_multi-FASTA.fa
cut -f1,2 Osa_organelle_multi-FASTA.fa.fai > Osa_organelle_gene_length.table

# rRNA
cp /data/1.public_data/genomes/NCBI/Oryza_sativa_GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_rna_from_genomic.fna.gz .
gunzip -d GCF_001433935.1_IRGSP-1.0_rna_from_genomic.fna.gz 
grep 'gbkey=rRNA' GCF_001433935.1_IRGSP-1.0_rna_from_genomic.fna | sed 's/^>\([^ ]*\).*/\1/' > Osa_ncbi_rRNA.idlist
seqtk subseq GCF_001433935.1_IRGSP-1.0_rna_from_genomic.fna Osa_ncbi_rRNA.idlist > Osa_ncbi_rRNA.fasta
bwa index Osa_ncbi_rRNA.fa

# nuclear genome
## fasta
cp /data/1.public_data/genomes/Ensembl/Oryza_sativa_IRGSP-1.0/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz .
gunzip -d Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz
seqtk subseq Oryza_sativa.IRGSP-1.0.dna.toplevel.fa genome_chr.idlist > Osa_genome_chr.fa
hisat2-build -p 3 Osa_genome_chr.fa Osa_genome_chr > Osa_genome_chr_hisat2-build.log 2>&1 &
## gtf
cp /data/1.public_data/genomes/Ensembl/Oryza_sativa_IRGSP-1.0/Oryza_sativa.IRGSP-1.0.53.gtf.gz .
gunzip Oryza_sativa.IRGSP-1.0.53.gtf.gz 
grep 'transcript_biotype "protein_coding"' Oryza_sativa.IRGSP-1.0.53.gtf > Osa_known_mRNA.gtf
