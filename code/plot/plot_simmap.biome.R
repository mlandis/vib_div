library(plotrix)
library(phytools)
source("/Users/mlandis/projects/vib_div/code/plot/phytools.mjl.R")

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    base_fn = "out.1.t163.f5.biome"
    #base_fn = "out.1.t163.f5.mask_fossil_states.pi_cold_25.biome"
} else if (length(args)==1) {
    base_fn = args[1]
}

# IO
fp = "/Users/mlandis/projects/vib_div/"
fn = paste("output/", base_fn, sep="")
phy_fn = paste(fp, fn, ".stoch_map.txt",  sep="")
plot_fp = paste(fp, "code/plot/", sep="")
col_fn = paste(plot_fp, "biome_colors.n4.txt",sep="")

# plotting settings
write_pdf = !FALSE
prompt = FALSE
if (write_pdf) {
    pdf( paste(plot_fp,"fig/stoch/",base_fn,".stoch_map.pdf",sep="") )
}

# read data
dat_col = read.csv(col_fn)
dat_ch = read.table(phy_fn, sep="\t", header=T)

# get phylo n_it, n_burn
phy = as.vector(dat_ch[,ncol(dat_ch)])

# read/plot/append simmap trees
simphy = list()
iterations = dat_ch[,1]
n_it = length(iterations)
burn = 0.75
n_burn = max(1, floor(n_it * burn))
    
#for (n in n_burn:n_it) {
#for (n in 1:5) {
for (j in n_burn:n_it) {
    
    n = which(dat_ch[,1]==iterations[j])
    cat("tree",j,":",iterations[j],"->",n,"\n")
    sim2 = read.simmap(text=phy[n])
    sim2 = ladderize.simmap(sim2, right=F)
    
    # find relevant colors
    colors = vector()
    for (i in 1:length( sim2$maps ) ) {
        colors = c(colors, names(sim2$maps[[i]]) )
    }
    colors = sort(as.numeric(unique(colors)))
    #print(colors)

    # convert state numbers to state labels    
    # +1 correction for base-0 index
    fix=1
    colnames(sim2$mapped.edge) = as.vector(dat_col$name[as.numeric(colnames(sim2$mapped.edge))+fix])
    for (j in 1:length(sim2$maps)) {
        if (length(sim2$maps[[j]]) > 0) {
            names(sim2$maps[[j]]) = as.vector(dat_col$name[as.numeric(names(sim2$maps[[j]]))+fix])
        }
    }
    simphy[[n]] = sim2
    
    #cols = setNames( as.vector(dat_col[colors,2]), dat_col[colors,1] )
    cols = as.vector(dat_col$color)[1:4]
    names(cols)=dat_col$name[1:4]

    # mar : bottom, left, top, and right    
    max_age = 100
    root_age = max(branching.times(sim2))
    xlim = c(-3,max_age)
    dx = max(branching.times(sim2))
    plotSimmap.mjl(sim2, cols, fsize=0.3, lwd=3, xlim=xlim,
                   split.vertical=TRUE, direction="rightwards",
                   x_offset=max_age-root_age-4.5) #, mar=c(0.1,dx,0.1,0.1))
    #axis(1,side=3)

    lbl_names = c("Late\nCretaceous","Paleogene","Early\nEocene","Mid/Late\nEocene","Oligocene","Early\nMiocene","Mid/Late\nMiocene","Recent")
    #x_names = max_age-c(90, 65, 56, 48, 33.9, 23, 16, 5.3)-4.5
    x_pos = max_age-c(90, 65, 56, 48, 33.9, 23, 16, 5.3, 0)-4.5
    x_mid = (x_pos[1:(length(x_pos)-1)] + x_pos[2:length(x_pos)])/2
    for (k in 2:length(x_pos)) {
        y1=168
        if (k==length(x_pos)) {
            y1=0
        }
        segments(x0=x_pos[k], x1=x_pos[k], y0=-3, y1=y1, lty=3, col="pink2", lwd=2)
    }
    for (k in 1:length(x_mid)) {
        text(x=x_mid[k], y=-2, labels=lbl_names[k], cex=0.45 )
    }
       
    add.simmap.legend(colors=cols,
                      vertical=TRUE,
                      shape="circle",
                      prompt=F,
                      y=166, fsize=0.8)

}
if (!write_pdf) {
    invisible(readline(prompt="continue?"))
} else {
    dev.off()
}
