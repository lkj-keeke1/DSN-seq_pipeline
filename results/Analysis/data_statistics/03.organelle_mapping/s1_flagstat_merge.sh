datadir=~/projects/DSN-seq/results/P1_call_editing/01.mapping

ln -s ${datadir}/*.flagstat .
ls *.flagstat | while read id; do awk 'BEGIN{print ARGV[1]}{print $1}' $id > $id.tmp; done
ls *.tmp | while read id; do sed -i s'/.flagstat//' $id; done

paste *.tmp > all_mapping.summary.tmp0
rm *.tmp 
rm *.flagstat

mv all_mapping.summary.tmp0 RNA_editing_organelle_mapping.summary
