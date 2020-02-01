library(plotrix)
library(phytools)
source("vib_div_util.R")

base_fn = "out.1.t163.f5.bg"

# IO
fp      = "../../"
fn      = paste0("output/", base_fn)
phy_fn  = paste0(fp, fn, ".stoch_map.txt")
plot_fp = paste0(fp, "code/plot/")
col_fn  = paste0(plot_fp, "range_colors.n6.txt")

# plotting settings
pdf( paste(plot_fp,"fig/fig_SX",base_fn,".stoch_map.pdf",sep="") )

# read data
dat_col = read.table(col_fn,sep=",",header=T)
dat_ch = read.table(phy_fn, sep="\t", header=T)

# get phylo n_it, n_burn
phy = as.vector(dat_ch[,ncol(dat_ch)])

# iterations to sample
iterations = dat_ch[,1]
burn = 0.99
n_it = length(iterations)
n_burn = max(1, floor(n_it * burn))

# read/plot/append simmap trees
simphy = list()
for (j in n_burn:n_it) {
    
    n = which(dat_ch[,1]==iterations[j])
    cat("tree",j,":",iterations[j],"->",n,"\n")
    sim2 = read.simmap(text=phy[n])
    sim2 = ladderize.simmap(sim2, right=F)
    sim2 = fix_vib_tip(sim2)
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
    
    cols = as.vector(dat_col$color)
    names(cols)=dat_col$name

    # mar : bottom, left, top, and right    
    max_age = 100
    root_age = max(branching.times(sim2))
    xlim = c(-3,max_age)
    dx = max(branching.times(sim2))

    plotSimmap.vib_div(sim2, cols, fsize=0.3, lwd=3, xlim=xlim,
               split.vertical=TRUE, direction="rightwards",
               x_offset=max_age-root_age-4.5) #, mar=c(0.1,dx,0.1,0.1))
    
    # lbl_names = c("Late\nCretaceous","Paleocene","Early\nEocene","Mid/Late\nEocene","Oligocene","Early\nMiocene","Mid/Late\nMiocene","Recent")
    # #x_names = max_age-c(90, 65, 56, 48, 33.9, 23, 16, 5.3)-4.5
    # x_pos = max_age-c(90, 65, 56, 48, 33.9, 23, 16, 5.3, 0)-4.5
    # x_mid = (x_pos[1:(length(x_pos)-1)] + x_pos[2:length(x_pos)])/2
    # for (k in 2:length(x_pos)) {
    #     y1=168
    #     if (k==length(x_pos)) {
    #         y1=0
    #     }
    #     segments(x0=x_pos[k], x1=x_pos[k], y0=-3, y1=y1, lty=3, col="pink2", lwd=2)
    # }
    # for (k in 1:length(x_mid)) {
    #     text(x=x_mid[k], y=-2, labels=lbl_names[k], cex=0.45 )
    # }
       
    add.simmap.legend(colors=cols[1:6],
                      vertical=TRUE,
                      shape="circle",
                      prompt=F,
                      y=166, fsize=0.8)
    
    add.simmap.legend(colors=cols[7:length(cols)],
                      vertical=TRUE,
                      shape="circle",
                      prompt=F,
                      y=166-25, fsize=0.8)
}

dev.off()
