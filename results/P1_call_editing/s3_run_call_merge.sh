source activate DSN-seq
cat ./sample.merge.list | while read wkdir sample genome length
do 
bash call_only.sh $wkdir $sample $genome $length &
done

