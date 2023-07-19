source activate DSN-seq
cat ../P0_data_qc/02.cutadapt/sample.tail_cut.list | while read wkdir sample R1 R2 hisat_ref genome length
do 
bash call_editing_sites.sh $wkdir $sample $R1 $R2 $hisat_ref $genome $length &
done

