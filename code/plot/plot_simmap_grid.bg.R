library(plotrix)
library(phytools)
source("vib_div_util.R")

#base_fn = "out.1.t163.f5.mask_fossil_states.bg"
#base_fn = "out.1.t163.f5.bg"

# IO
fp      = "../../"
fn      = paste0("output/", base_fn)
phy_fn  = paste0(fp, fn, ".bg.stoch_map.txt")
col_fn  = paste0(fp, "code/plot/range_colors.n6.txt")
plot_fn = paste0(fp, "code/plot/fig/figSX_",base_fn,".stoch_grid.pdf",sep="")

# plotting settings
#pdf(plot_fn, height=17, width=3.5)
grid = c(length(iterations),1)
par(mfrow=grid)

# read data
dat_col = read.table(col_fn,sep=",",header=T)
dat_ch  = read.table(phy_fn, sep="\t", header=T)

# get vector of stochastically mapped trees
phy = as.vector(dat_ch[,ncol(dat_ch)])

# iterations to sample
#iterations = c(99000, 99050, 99100, 99150)
n_it = length(iterations)

# read/plot/append simmap trees
simphy = list()
for (j in 1:n_it) {

    # read in individual tree
    n = which(dat_ch[,1]==iterations[j])
    cat("tree",j,":",iterations[j],"->",n,"\n")
    sim2 = read.simmap(text=phy[n])
    sim2 = ladderize.simmap(sim2, right=F)
    sim2$tip.label = rep("",length(sim2$tip.label))
    
    # find relevant colors
    colors = vector()
    for (i in 1:length( sim2$maps ) ) {
        colors = c(colors, names(sim2$maps[[i]]) )
    }
    colors = sort(as.numeric(unique(colors)))

    # convert state numbers to state labels    
    colnames(sim2$mapped.edge) = as.vector(dat_col$name[as.numeric(colnames(sim2$mapped.edge))])
    for (j in 1:length(sim2$maps)) {
        if (length(sim2$maps[[j]]) > 0) {
            names(sim2$maps[[j]]) = as.vector(dat_col$name[as.numeric(names(sim2$maps[[j]]))])
        }
    }
    simphy[[n]] = sim2
    
    # assign state labels to colors
    cols = as.vector(dat_col$color)
    names(cols)=dat_col$name

    # mar : bottom, left, top, and right    
    root_age = max(branching.times(sim2))
    xlim = c(0,root_age)
    dx = max(branching.times(sim2))
    plotSimmap(sim2, cols, fsize=0.3, lwd=2, xlim=xlim,
                   split.vertical=TRUE, direction="rightwards")
}

#dev.off()
