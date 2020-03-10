These files are needed to run the RevBayes scripts located in `../code/`. The file `../code/run_quick.Rev` controls exactly which input files are used and how they're used. A brief description of the data files:

- `viburnum.mol.nex` contains the cpDNA + ITS sequence alignment
- `viburnum.range.n6.nex` contains the 6-region range data (NOTE: these are integer-coded states, because we wanted to allow for ambiguous ranges for some fossil taxa; see RevBayes tutorials on biogeography to use presence-absence ranges)
- `viburnum.biome.n4.nex` contains the 4-biome state data
- `viburnum.backbone.tre` contains the Stage 1 RAD-seq topology constraint and the morphological constraints for the 10 unsequenced, extant taxa
- `viburnum.init.tre` contains a reasonable starting tree to speed up mixing
- `viburnum.fossil_intervals.tsv` contains the age ranges for fossil taxa
- `viburnum.taxa.tsv` contains all taxon labels and fossil taxon ages
- `viburnum.area_graph.n6.*.csv` contains the adjacencies among areas over time
- `viburnum.bg.times.txt` contains the times for applying each area-graph
