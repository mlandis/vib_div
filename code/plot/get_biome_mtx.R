library(diagram)

make_mtx = function(fn) {
    # read trace
    dat =read.table(fn,header=T)
    r_idx = grep( "r_biome", colnames(dat) )
    r_means = colMeans(dat[,r_idx])
    rates = c(          NA,  r_means[1],  r_means[2], r_means[3],
                r_means[4],          NA,  r_means[5], r_means[6],
                r_means[7],  r_means[8],          NA, r_means[9],
               r_means[10], r_means[11], r_means[12],         NA)
    
    Q = matrix(  rates, byrow=T, nrow=4, ncol=4)
    #Q = Q / min(Q,na.rm=T)
    colnames(Q)<-rownames(Q)<-c("H","L","C","T")
    return(Q)
}

plot_mtx = function(fn, Q, round_digits=1, title_str="") {
    
    colors = c("red","darkgreen","dodgerblue","darkblue")
    color_mat = matrix( rep(colors,4), nrow=4)
    coords = matrix( c(-0.00,  0.50,
                        0.50, -0.15,
                        0.50,  1.15,
                        1.00,  0.50), byrow=T, nrow=4, ncol=2)
    colnames(coords)<-c("x","y")
    rownames(coords)<-c("H","L","C","T")
    
    tQ = t(Q)
    tQ = round(tQ, round_digits)

    pdf( fn, height=5, width=5)
    plotmat( tQ, curve=0.2, name=c("H","L","C","T"), lwd=tQ, relsize=0.65, pos=coords,
             box.lcol=colors, box.lwd=1, shadow.size=0,
             arr.lcol=color_mat, arr.col=color_mat,
             arr.type="triangle", arr.pos=0.6, arr.length=0.3, arr.width=0.3, endhead=F,
             main=title_str
             )
    
    dev.off()
}

# files
fp = "/Users/mlandis/projects/vib_div/"
out_fp = paste(fp, "output/", sep="")
fn1 = "out.1.t163.f5"
fn2 = "out.1.t163.f5.mask_fossil_states"
fnd = "out.diff.without_vs_with_ratio"
out_fn1 = paste(out_fp, fn1, ".model.log", sep="")
out_fn2 = paste(out_fp, fn2, ".model.log", sep="")
plot_fn1 = paste(fp, "code/plot/fig/matrix/", fn1, ".biome_matrix.pdf", sep="")
plot_fn2 = paste(fp, "code/plot/fig/matrix/", fn2, ".biome_matrix.pdf", sep="")
plot_fnd = paste(fp, "code/plot/fig/matrix/", fnd, ".biome_matrix.pdf", sep="")

Q1 = make_mtx(fn=out_fn1)
Q2 = make_mtx(fn=out_fn2)

Q1 = Q1 / Q1[1,4]
Q2 = Q2 / Q2[1,4]
Qd = Q2 / Q1

plot_mtx( plot_fn1, Q1, round_digits=2, title_str="Rates with fossil states" )
plot_mtx( plot_fn2, Q2, round_digits=2, title_str="Rates without fossil states" )
plot_mtx( plot_fnd, Qd, round_digits=2, title_str="Ratio of rates without/with fossil states")


