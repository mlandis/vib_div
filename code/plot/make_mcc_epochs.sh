PIDS=()

echo "Filtering diff. rogue taxon sets..."
Rscript --vanilla drop_fossils.R out.1.t163.f5

echo "Processing MCC trees..."
sleep 0.5
echo "fn=\"181024/out.1.t163.f5\";make_states=false;source(\"make_mcc.Rev\");" | rb &
PIDS+=($!)
echo "fn=\"181024/out.1.t163.f5.mask_fossil_states\";make_states=false;source(\"make_mcc.Rev\");" | rb &
PIDS+=($!)
echo "fn=\"181024/out.2.t163.f5\";make_states=false;source(\"make_mcc.Rev\");" | rb &
PIDS+=($!)
echo "fn=\"181024/out.2.t163.f5.mask_fossil_states\";make_states=false;source(\"make_mcc.Rev\");" | rb &
PIDS+=($!)
echo "fn=\"181024/out.1.t163.f5.no_fossil\";make_states=false;source(\"make_mcc.Rev\");" | rb &
PIDS+=($!)
echo "fn=\"181024/out.1.t163.f5.no_fossil_unsequenced\";make_states=false;source(\"make_mcc.Rev\");" | rb &
PIDS+=($!)
echo "fn=\"181024/out.1.t163.f5.radseq_only\";make_states=false;source(\"make_mcc.Rev\");" | rb &
PIDS+=($!)
echo "fn=\"181024/out.1.t163.f5.no_epoch\";make_states=false;source(\"make_mcc.Rev\");" | rb &
PIDS+=($!)
echo "fn=\"181024/out.2.t163.f5.no_epoch\";make_states=false;source(\"make_mcc.Rev\");" | rb &
PIDS+=($!)

for pid in "${PIDS[@]}"
do
    echo "Waiting for ${pid}"
    wait $pid
done


wait
echo "...MC trees done!"

echo "Processing ancestral state trees..."
Rscript --vanilla plot_mcc.R 181024/out.1.t163.f5.no_epoch.mcc.tre
Rscript --vanilla plot_mcc.R 181024/out.2.t163.f5.no_epoch.mcc.tre

Rscript --vanilla plot_mcc.R 181024/out.1.t163.f5.mcc.tre
Rscript --vanilla plot_mcc.R 181024/out.1.t163.f5.map.tre
Rscript --vanilla plot_mcc.R 181024/out.1.t163.f5.mask_fossil_states.mcc.tre
Rscript --vanilla plot_mcc.R 181024/out.1.t163.f5.mask_fossil_states.map.tre
Rscript --vanilla plot_mcc.R 181024/out.2.t163.f5.mcc.tre
Rscript --vanilla plot_mcc.R 181024/out.2.t163.f5.map.tre
Rscript --vanilla plot_mcc.R 181024/out.2.t163.f5.mask_fossil_states.mcc.tre
Rscript --vanilla plot_mcc.R 181024/out.2.t163.f5.mask_fossil_states.map.tre

Rscript --vanilla plot_mcc.R 181024/out.1.t163.f5.no_fossil.mcc.tre
Rscript --vanilla plot_mcc.R 181024/out.1.t163.f5.no_fossil.map.tre
Rscript --vanilla plot_mcc.R 181024/out.1.t163.f5.no_fossil_unsequenced.mcc.tre
Rscript --vanilla plot_mcc.R 181024/out.1.t163.f5.no_fossil_unsequenced.map.tre
Rscript --vanilla plot_mcc.R 181024/out.1.t163.f5.radseq_only.mcc.tre
Rscript --vanilla plot_mcc.R 181024/out.1.t163.f5.radseq_only.map.tre
echo "...done!"
