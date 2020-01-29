library(phyloch)
library(strap)
library(OutbreakTools)

out_fp = "/Users/mlandis/projects/viburnum_phylogeo/"
out_str = "output/biogeo/out.viburnum.t165.f5.w_bg_5.w_biome_4.mcc.tre"

tree_fn = paste(out_fp, out_str, sep="")

t <- read.annotated.nexus(tree_fn)
t$root.time <- max(branching.times(t))
print(t$root.time)

HPDbars(t)
geoscalePhylo(tree=ladderize(t,right=FALSE), units=c("Period", "Epoch"), boxes="Epoch",
              cex.tip=0.5, cex.age=0.7, cex.ts=0.7,
              label.offset=0, x.lim=c(-15,80), lwd=3, width=2)