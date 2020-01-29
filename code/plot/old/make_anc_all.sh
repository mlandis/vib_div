FN_PREFIX=$1
#FN_SUFFIXES=("" ".w_fossils" ".w_fossils.mask_biomes" ".w_fossils.mask_fossils" ".br_indep" ".w_fossils.br_indep" ".w_fossils.mask_biomes.br_indep" ".w_fossils.mask_fossils.br_indep")
FN_SUFFIXES=("" ".w_fossils" ".w_fossils.mask_biomes" ".w_fossils.mask_fossils") #".br_indep" ".w_fossils.br_indep" ".w_fossils.mask_biomes.br_indep" ".w_fossils.mask_fossils.br_indep")

FP="/Users/mlandis/projects/biome_range/"

#echo ${FN_SUFFIX[0]}
#echo ${FN_SUFFIX[1]}

echo "Processing ancestral states"
for FN_SUFFIX in "${FN_SUFFIXES[@]}"
do
    #echo "ok"
    echo "fn_prefix=\"$FN_PREFIX\"; fn_suffix=\"$FN_SUFFIX\"; source(\"make_mcc.Rev\");" | rb &
done
echo "...done!"

wait

echo "Generating figures"
for FN_SUFFIX in "${FN_SUFFIXES[@]}"
do
    FN=${FN_PREFIX}${FN_SUFFIX}".ase.tre"
    cp $FP"output/"${FN} $FP"code/plot/fig/"${FN}
    Rscript --vanilla plot_anc.R ${FN_PREFIX}${FN_SUFFIX} &
done
echo "...done!"

