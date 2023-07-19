ln -s ../../../../P1_call_editing/01.mapping/Osa*.idxstats .
echo -e "gene\tlength" > annotation.table
cut -f 1,2 Osa_DSN_rep1.filtered.idxstats >> annotation.table
ls *.idxstats | while read id
do
	sid=${id%.*}
	echo ${sid%.*} > ${id}.tmp
	cut -f 3 ${id} >> ${id}.tmp &
done
