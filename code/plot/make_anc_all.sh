FN=( "out.1.t163.f5" "out.1.t163.f5.mask_fossil_states" "out.2.t163.f5" "out.2.t163.f5.mask_fossil_states" "out.3.t163.f5.n_biomes_2" "out.3.t163.f5.mask_fossil_states.n_biomes_2" "out.3.t163.f5.n_biomes_2" "out.3.t163.f5.mask_fossil_states.n_biomes_2" )


"Processing ASE only..."

for fn in ${FN[@]}
do
    Rscript --vanilla plot_biome.R $fn &
    Rscript --vanilla --max-ppsize=500000 plot_bg.R $fn &
    #wait 
done

echo "...done!"
