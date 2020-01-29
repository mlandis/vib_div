library(ape)
library(geiger)

if (T) {
    fp = "/Users/mlandis/projects/biome_range/"
    fp2 = ""
    fn = "lawn1.w_fossils.tre"
} else {
    fp = "/Users/mlandis/projects/vib_div/"
    fp2 = "job_180524/"
    fn = "out.2.t163.f5.tre"
}
 
    
back_fn = paste(fp, "data/", "viburnum.backbone.tre", sep="")
post_fn = paste(fp, "output/", fp2, fn, sep="")
plot_fn = paste(fp, "code/plot/fig/", fn, ".tree_constraints.pdf", sep="")
post = read.tree(post_fn)
back = read.tree(back_fn)


pdf(plot_fn)
for (i in 1:length(post)) {
    phyp = post[[i]]
    phyb = back[[4]]
    
    keep_tips = phyp$tip.label[ -match(phyb$tip.label, phyp$tip.label) ]
    phyl = drop.tip( phy=phyp, tip=keep_tips )
    phyl = ladderize(phyl)
    
    plot(phyl)
}
dev.off()
