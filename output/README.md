Two independent MCMC chains were run (`out.1` and `out.2`) for each of two main analysis settings (Complete [no special label] or Masked [labeled with `mask_fossil_states`]). For the Complete setting, the tree distributions underwent three forms of taxon filtering (`no_fossil` filtered out the 5 fossil taxa; `no_fossil_unsequenced` filtered out the 5 fossil and 10 morphologically placed taxa; `radseq_only` filtered out all taxa except those with RADseq). All output is based on 2001 posterior samples that has already had burn-in removed. Several types of output were generated for each MCMC analysis, which are most easily distinguished by their file suffixes.

File types and suffixes:
- `model.log` is the posterior distribution of model parameters
- `biome.states.txt` is the ancestral biome estimates indexed by node
- `biome.stoch_map.txt` contains the stochastic mapping for ancestral biomes
- `biome.history.tsv` is the table of stochastic mappings for ancestral biomes
- `bg.states.txt` is the ancestral range estimates indexed by node
- `bg.stoch_map.txt` contains the stochastic mappings for ancestral ranges
- `bg.history.tsv` is the table of stochastic mappings for ancestral ranges
- `trees` is the posterior distribution of trees
- `mcc.tre` is the maximum clade credibility (MCC) tree
- `map.tre` is the maximum a posteriori tree
- `con.tre` is the majority rule consensus tree (p=0.5)
- `biome.ase.tre` is the MCC tree mapped with ancestral biome estimates
- `bg.ase.tre` is the MCC tree mapped with ancestral range estimates
