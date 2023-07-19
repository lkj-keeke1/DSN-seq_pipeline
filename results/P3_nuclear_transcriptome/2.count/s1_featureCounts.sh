source activate DSN-seq 
featureCounts -T 32 -s 1 -p -M -g transcript_id  -a ../../../data/ref/Osa/Oryza_sativa.IRGSP-1.0.53.gtf -o all.table ../1.mapping/Osa*.sorted.bam 1>all.log 2>&1 
awk '$1!~/#/' all.table | sed 's/\t\([^\t]*\)\/\([^\/]*\).sorted.bam/\t\2/g' > all.counts
Rscript ~/projects/DSN-seq/src/run_normalization.R all.counts all

