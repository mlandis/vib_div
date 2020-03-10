import dendropy as dp
import pandas as pd

fp = '/Users/mlandis/projects/gh_vib_div/'
phy_fn = fp + 'output/out.2.t163.f5.mask_fossil_states.trees'
#phy_fn = fp + 'output/out.1.t163.f5.trees'

f_burn = 0.0

clade_fn = fp + 'code/plot/fossil_clade.tre'
clade_full = dp.Tree.get_from_path(src=clade_fn, schema='newick')
#clade_no_BC = clade_full; clade_no_BC.remove('Viburnum BC')
#clade_no_CO = clade_full; clade_no_CO.remove('Porphyrotinus CO')

clade_por_CO_BC = [ t.taxon.label for t in clade_full.leaf_nodes() ]
clade_por_CO = [ t.taxon.label for t in clade_full.leaf_nodes() if t.taxon.label != 'Viburnum BC' ]
clade_por_BC = [ t.taxon.label for t in clade_full.leaf_nodes() if t.taxon.label != 'Porphyrotinus CO' ]
clade_euviburnum = [ "V mongolicum", "V burejaeticum", "V lantana", "V maculatum", "V glomeratum", "V veitchii", "V rhytidophyllum", "V buddleifolium", "V chinshanense", "Valvatotinus IS" ]
clade_lentago = [ "V cassinoides", "V nudum", "V lentago", "V elatum", "V obovatum", "V prunifolium", "V rufidulum", "Valvatotinus IS" ]
clade_nudum = [ "V cassinoides", "V nudum", "Valvatotinus IS" ]



clade_labels = { 'Por+CO+BC':clade_por_CO_BC, 
                 'Por+CO':clade_por_CO,
                 'Por+BC':clade_por_BC,
                 'Eu+IS':clade_euviburnum,
                 'Len+IS':clade_lentago,
                 'Nud+IS':clade_nudum }


print('Reading trees...')
x = pd.read_csv(phy_fn, sep='\t')
x_tree = x['tree']
n_burn = int(f_burn * len(x_tree))
x_tree = x_tree[n_burn:]
s_phy = ''.join( x_tree + '\n' )

print('Creating tree objects...')
phy = dp.TreeList.get_from_string(src=s_phy, schema='newick')

print('Searching for clade probs...')
clade_pp = {}
for k,v in clade_labels.items():
    print('    ' + k + ' = ', end='')
    clade_pp[k] = phy.frequency_of_bipartition( labels=v )
    print( str( clade_pp[k]) )

print('done!')

print(clade_pp)

