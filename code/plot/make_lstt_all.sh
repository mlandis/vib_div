PIDS=()

echo "Generating stochastic mapping tsv files..."
sleep 0.5
echo "fn=\"out.1.t163.f5\";source(\"make_stoch_tsv.Rev\");" | rb &
PIDS+=($!)
echo "fn=\"out.1.t163.f5.mask_fossil_states\";source(\"make_stoch_tsv.Rev\");" | rb &
PIDS+=($!)
echo "fn=\"out.2.t163.f5\";source(\"make_stoch_tsv.Rev\");" | rb &
PIDS+=($!)
echo "fn=\"out.2.t163.f5.mask_fossil_states\";source(\"make_stoch_tsv.Rev\");" | rb &
PIDS+=($!)

for pid in "${PIDS[@]}"
do
    echo "Waiting for ${pid}"
    wait $pid
done


wait
echo "...tsv files done!"

echo "Plotting LSTT bar chars..."
PIDS2=()
Rscript --vanilla plot_lstt.R out.1.t163.f5 &
PIDS2+=($!)
Rscript --vanilla plot_lstt.R out.2.t163.f5 &
PIDS2+=($!)
Rscript --vanilla plot_lstt.R out.1.t163.f5.mask_fossil_states &
PIDS2+=($!)
Rscript --vanilla plot_lstt.R out.2.t163.f5.mask_fossil_states &
PIDS2+=($!)
for pid in "${PIDS2[@]}"
do
    echo "Waiting for ${pid}"
    wait $pid
done

echo "Plotting LSTT curves..."
PIDS2=()
Rscript --vanilla plot_lstt_grid.R out.1.t163.f5 &
PIDS2+=($!)
Rscript --vanilla plot_lstt_grid.R out.2.t163.f5 &
PIDS2+=($!)
Rscript --vanilla plot_lstt_grid.R out.1.t163.f5.mask_fossil_states &
PIDS2+=($!)
Rscript --vanilla plot_lstt_grid.R out.2.t163.f5.mask_fossil_states &
PIDS2+=($!)
for pid in "${PIDS2[@]}"
do
    echo "Waiting for ${pid}"
    wait $pid
done
echo "...done!"
