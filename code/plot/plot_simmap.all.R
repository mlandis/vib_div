library(plotrix)
library(phytools)
source("/Users/mlandis/projects/vib_div/code/plot/phytools.mjl.R")


# IO
fp = "/Users/mlandis/projects/vib_div/"
fn = "output/tinus/out.1.viburnum.t163.f5.bg"
phy_fn = paste(fp, fn, ".stoch_map.txt",  sep="")
col_fn = paste(fp, "code/plot/range_colors.txt",sep="")

# plotting settings
write_pdf = !FALSE
prompt = FALSE
if (write_pdf) {
    pdf( paste(fp,fn,".stoch_map.pdf",sep="") )
}

# read data
dat_col = read.csv(col_fn)
dat_ch = read.table(phy_fn, sep="\t", header=T)

# get phylo n_it, n_burn
phy = as.vector(dat_ch[,ncol(dat_ch)])

# read/plot/append simmap trees
simphy = list()
#for (n in n_burn:n_phy) {
#for (n in 10:30) {
iterations = dat_ch[,1]
#iterations = iterations[ which(iterations %in% c(8400, 8500, 8600)) ]
#iterations = seq(10000,10400,100)
#iterations = c(10300)
n_it = length(iterations)
burn = 0.50
n_burn = max(1, floor(n_it * burn))

#for (n in n_burn:n_it) {
#for (n in 1:5) {
for (j in n_burn:n_it) {
    
    n = which(dat_ch[,1]==iterations[j])
    cat("tree",j,":",iterations[j],"->",n,"\n")
    sim2 = read.simmap(text=phy[n])
    sim2 = ladderize.simmap(sim2)
    
    # find relevant colors
    colors = vector()
    for (i in 1:length( sim2$maps ) ) {
        colors = c(colors, names(sim2$maps[[i]]) )
    }
    colors = sort(as.numeric(unique(colors)))
    #print(colors)

    # convert state numbers to state labels    
    colnames(sim2$mapped.edge) = as.vector(dat_col[as.numeric(colnames(sim2$mapped.edge)),1])
    for (j in 1:length(sim2$maps)) {
        if (length(sim2$maps[[j]]) > 0) {
            names(sim2$maps[[j]]) = as.vector(dat_col[as.numeric(names(sim2$maps[[j]])),1])
        }
    }
    simphy[[n]] = sim2
    
    #cols = setNames( as.vector(dat_col[colors,2]), dat_col[colors,1] )
    cols = as.vector(dat_col$color)
    names(cols)=dat_col$range

    # mar : bottom, left, top, and right    
    max_age = 80
    root_age = max(branching.times(sim2))
    xlim = c(-3,max_age)
    #print(xlim)
    #dx = xlim[2] - max(branching.times(sim2))
    dx = max(branching.times(sim2))
    #print(dx)
    #print(dx)
    #mar = c(0.1, dx, 0.1, 0.1)
    #print(mar)
          
    #plotSimmap(sim2, cols, fsize=0.4, lwd=3, split.vertical=TRUE, direction="rightwards") #, xlim=c(0,15), mar=c(0.1,0.1,0.1,2)) #, xlim=xlim)
    #plotSimmap(sim2, cols, fsize=0.4, lwd=3, xlim=xlim, split.vertical=TRUE, direction="rightwards") #, xlim=c(0,15), mar=c(0.1,0.1,0.1,2)) #, xlim=xlim)
    plotSimmap.mjl(sim2, cols, fsize=0.4, lwd=3, xlim=xlim,
                   split.vertical=TRUE, direction="rightwards",
                   x_offset=max_age-root_age-4.5) #, mar=c(0.1,dx,0.1,0.1))
    #axis(1,side=3)
   
    n_row = 3
    n_col = floor(length(cols)/n_row)
    for (j in 0:(n_col-1)) {
        k = j * n_row
        #print(k+(1:n_row))
        #print(cols[k+(1:n_row)])
        add.simmap.legend(colors=cols[k+(1:n_row)],
                          vertical=TRUE,
                          shape="circle",
                          prompt=prompt,
                          x=-3+3*k, y=10)

    }
    #add.simmap.legend(colors=cols[(n_col+1):length(cols)],
    #                  vertical=TRUE,
    #                  shape="circle",
    #                  prompt=prompt,
    #                  x=0, y=10)
}
if (!write_pdf) {
    invisible(readline(prompt="continue?"))
} else {
    dev.off()
}
