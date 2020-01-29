
echo "
******************
* Make MCC plots * 
******************"
./make_mcc_all.sh

echo "
******************
* Make age plots *
******************"
./make_ages.sh

echo "
********************
* Make Pr(T) plots *
********************"
Rscript biome_prior.R

echo "
*********************
* Make SIMMAP plots *
*********************"
./make_stoch_all.sh

echo "
*******************
* Make LSTT plots *
*******************"
./make_lstt_all.sh

echo "
******************
* Make ASE plots *
******************"
./make_ase_all.sh
./make_anc_all.sh
