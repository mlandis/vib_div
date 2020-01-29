FP="/Users/mlandis/projects/vib_div/code/plot/fig/stoch/"
FN=${FP}"out.1.t163.f5"

convert -verbose -delay 20 -loop 0 -dispose 2 -density 150 ${FN}.biome.stoch_map.pdf ${FN}.biome.stoch_map_anim.gif &
convert -verbose -delay 20 -loop 0 -dispose 2 -density 150 ${FN}.bg.stoch_map.pdf ${FN}.bg.stoch_map_anim.gif &
convert -verbose -delay 20 -loop 0 -dispose 2 -density 150 ${FN}.mask_fossil_states.biome.stoch_map.pdf ${FN}.mask_fossil_states.biome.stoch_map_anim.gif &
convert -verbose -delay 20 -loop 0 -dispose 2 -density 150 ${FN}.mask_fossil_states.bg.stoch_map.pdf ${FN}.mask_fossil_states.bg.stoch_map_anim.gif &
