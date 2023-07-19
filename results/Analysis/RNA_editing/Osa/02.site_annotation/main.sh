source activate DSN-seq
Rscript src/s1_gene_info_annotation.R
bash src/s2_codon_seq_tk.sh
Rscript src/s3_aa_change_annotation.R
Rscript src/s4_canonical_annotation\(optional\).R
bash src/s5.0_rate_depth_annotation_file_list.sh
Rscript src/s5.1_rate_depth_annotation_filtered.R
Rscript src/s5.2_rate_depth_annotation_whole.R
#bash src/s6_rm_tmp.sh
