ls ~/projects/DSN-seq/results/P1_call_editing/02.call_editing/Osa*.out.table > sample.call.list
sed -i 's/.out.table//' sample.call.list 
