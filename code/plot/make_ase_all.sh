PIDS=()

echo "Processing MCC/ASE trees..."
sleep 0.5
echo "fn=\"out.1.t163.f5\";source(\"make_mcc.Rev\");" | rb &
PIDS+=($!)
echo "fn=\"out.1.t163.f5.mask_fossil_states\";source(\"make_mcc.Rev\");" | rb &
PIDS+=($!)
echo "fn=\"out.2.t163.f5\";source(\"make_mcc.Rev\");" | rb &
PIDS+=($!)
echo "fn=\"out.2.t163.f5.mask_fossil_states\";source(\"make_mcc.Rev\");" | rb &
PIDS+=($!)

for pid in "${PIDS[@]}"
do
#    echo "Waiting for ${pid}"
    wait $pid
done

