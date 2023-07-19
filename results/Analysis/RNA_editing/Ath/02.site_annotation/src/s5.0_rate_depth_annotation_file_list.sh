ls ~/projects/DSN-seq/results/P1_call_editing/02.call_editing/Ath*.out.table > sample.call.list
ls ~/projects/DSN-seq/results/sup_STS/02.call_editing/Ath*.out.table >> sample.call.list
sed -i 's/.out.table//' sample.call.list 
