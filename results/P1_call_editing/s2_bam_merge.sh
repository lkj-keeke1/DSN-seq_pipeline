cd 01.mapping

source activate DSN-seq
samtools merge -@ 10 Osa_DSN_merge.filtered.sorted.bam Osa_DSN_rep*.filtered.sorted.bam
samtools merge -@ 10 Osa_PolyA_merge.filtered.sorted.bam Osa_PolyA_rep*.filtered.sorted.bam
samtools merge -@ 10 Osa_Ribo-off_merge.filtered.sorted.bam Osa_Ribo-off_rep*.filtered.sorted.bamsamtools merge -@ 10 Ath_DSN_merge.filtered.sorted.bam Ath_DSN_rep*.filtered.sorted.bam


samtools index Osa_DSN_merge.filtered.sorted.bam &
samtools index Osa_Ribo-off_merge.filtered.sorted.bam &
samtools index Osa_PolyA_merge.filtered.sorted.bam &
samtools index Ath_DSN_merge.filtered.sorted.bam &
