library(ggplot2)
library(Cairo)

fp = "/Users/mlandis/projects/vib_div/"
out_fp = paste(fp, "output/", sep="")
plot_fp = paste(fp, "code/plot/fig/biome_prior/", sep="")
plot_fn = paste(plot_fp, "fig6_cold_root_sensitivity.pdf", sep="")
#plot_fn2 = paste(plot_fp, "fig_temperate_entropy.pdf", sep="")

files = list.files(out_fp, full.names=T)
files = files[ grep("biome.states.txt", files) ]
files = files[ grep("pi_cold_", files) ]
files_nf = files[ grep("mask", files) ]
#files_wf = files[ !grep("mask", files) ]
#list.files(fp, pattern="mask", full.names=T)
#files_nf = files_nf[ grep("pi_cold_", files_nf) ]
files_wf = files[ !(files %in% files_nf) ]
n_files = length(files_nf)
f_burn = 0.5


pi_temp_prob = function(dat) {
    x = dat[,ncol(dat)]
    #print(x)
    p = sum(x==3) / length(x)
    return(p)
}

compute_entropy = function(dat) {
    x = dat[,ncol(dat)]
    
    p = c()
    for (j in 1:4) {
        p[j] = sum(x==j)
    }
    p = p / length(x)
    p = p[ p != 0 ]
    
    return( -sum( p * log(p, 2) ) )
}

compute_entropy_exp = function(x) {
    z = 1-x
    p = c(z/3, z/3, z/3, x)
    p = p[ p != 0 ]
    
    return( -sum( p * log(p, 2) ) )
}



df = data.frame( )
for (i in 1:n_files) {
    fp_toks = strsplit(files_nf[i], split="/")[[1]]
    fn_toks = strsplit(fp_toks[length(fp_toks)], split="\\.")[[1]]
    pi_toks = as.numeric(strsplit(fn_toks[6], split="_")[[1]][3])
    
    dat_nf = read.table(files_nf[i], sep="\t", header=T)
    dat_wf = read.table(files_wf[i], sep="\t", header=T)
    
    dat_nf = dat_nf[(f_burn*nrow(dat_nf)):nrow(dat_nf),]
    dat_wf = dat_wf[(f_burn*nrow(dat_wf)):nrow(dat_wf),]
    
    #cat( nrow(dat_nf), nrow(dat_wf), "\n")
    
    pi_temp_nf = pi_temp_prob(dat_nf)
    pi_temp_wf = pi_temp_prob(dat_wf)
    h_wf = compute_entropy(dat_wf)
    h_nf = compute_entropy(dat_nf)
    h_pi = compute_entropy_exp(pi_toks/100)
    df = rbind(df, c(pi_toks/100, pi_temp_nf, pi_temp_wf, h_nf, h_wf, h_pi))
    #df = rbind(df, c(
    #dat_nf[,ncol(dat_nf)]
}
colnames(df) = c("pi_freeze", "no_fossils", "with_fossils", "entropy_nf", "entropy_wf", "entropy_pi")

df = df[ order(df$pi_freeze), ]

df = rbind(
        c(0, 0, 0, NA, NA, NA),
        df,
        c(1, 1, 1, NA, NA, NA))

p = ggplot(df, aes(x=pi_freeze, y=no_fossils))


f_wf_inv = approxfun( df$with_fossils, df$pi_freeze )
f_wf = approxfun(  df$pi_freeze, df$with_fossils )
f_nf_inv = approxfun( df$no_fossils, df$pi_freeze )
f_nf = approxfun(  df$pi_freeze, df$no_fossils )


p = p + geom_segment( mapping=aes(x=0.25, xend=0.25, y=-Inf, yend=f_wf(0.25)), lty=2, col='gray')
p = p + geom_segment( mapping=aes(x=-Inf, xend=0.25, y=f_wf(0.25), yend=f_wf(0.25)), lty=2, col='gray')
p_flip = f_wf_inv(0.5)
p = p + geom_segment( mapping=aes(x=p_flip, xend=p_flip, y=-Inf, yend=f_wf(p_flip)), lty=2, col='gray')
p = p + geom_segment( mapping=aes(x=-Inf, xend=p_flip, y=f_wf(p_flip), yend=f_wf(p_flip)), lty=2, col='gray')
p_sig = f_wf_inv(0.95)
p = p + geom_segment( mapping=aes(x=p_sig, xend=p_sig, y=-Inf, yend=f_wf(p_sig)), lty=2, col='gray')
p = p + geom_segment( mapping=aes(x=-Inf, xend=p_sig, y=f_wf(p_sig), yend=f_wf(p_sig)), lty=2, col='gray')


p = p + annotate(geom="text", x=0.25-0.03, y=f_wf(0.25)+0.05, label=expression( paste(p[fair])))
p = p + annotate(geom="text", x=p_flip-0.03, y=f_wf(p_flip)+0.05, label=expression( paste(p[flip])))
p = p + annotate(geom="text", x=p_sig+0.05, y=f_wf(p_sig)-0.04, label=expression( paste(p[sig])))
p = p + geom_point( mapping=aes(x=0.25, y=f_wf(0.25)), col="firebrick3")
p = p + geom_point( mapping=aes(x=p_flip, y=f_wf(p_flip)), col="firebrick3")
p = p + geom_point( mapping=aes(x=p_sig, y=f_wf(p_sig)), col="firebrick3")

p = p + geom_abline(  intercept=0, slope=1, colour="lightgray", lwd=0.25)
p = p + geom_line(data=df, mapping=aes(x=pi_freeze, y=no_fossils, colour="black"), alpha=1)
p = p + geom_line(data=df, mapping=aes(x=pi_freeze, y=with_fossils, colour="red"), alpha=1)
#p = p + xlab("Prob (X_root = Cold | X_unk)")
p = p + xlab( expression( paste(Prob,"( ",X[root]," = ",Cold," | ",X[aux]," )") ) )
p = p + ylab( expression( paste(Prob,"( ",X[root]," = ",Cold," | ",X[obs],", ",X[aux]," )") ) )
#p = p + ylab("Prob (X_root = Cold | X_obs, X_unk)")

p = p + theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"))
p = p + scale_color_manual( name="Dataset", labels=c("Fossil-masked", "Complete"), values=c("dodgerblue", "firebrick3") )
p
CairoPDF(file=plot_fn, height=4, width=6)
print(p)
dev.off()


#p = p + xlab("") + ylab("")
#p = p + scale_y_log10()
#p = p + geom_hline( yintercept=0.50, linetype="dashed", colour="gray")
#p = p + geom_hline( yintercept=0.05, linetype="dashed", colour="gray")
#p = p + geom_hline( yintercept=0.95, linetype="dashed", colour="gray")
#p = p + geom_vline( xintercept=0.25, linetype="dashed", colour="gray")

# 
# 
# p2 = ggplot(df, aes(x=pi_freeze, y=entropy_pi))
# p2 = p2 + geom_line(data=df, mapping=aes(x=pi_freeze, y=entropy_nf, colour="dodgerblue"))
# p2 = p2 + geom_line(data=df, mapping=aes(x=pi_freeze, y=entropy_wf, colour="firebrick3"))
# p2 = p2 + geom_line(data=df, mapping=aes(x=pi_freeze, y=entropy_pi, colour="black"))
# p2 = p2 + xlab("Prior (MRCA = Fr.Temp.)")
# p2 = p2 + ylab("Entropy")
# p2 = p2 + theme(panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.background = element_blank(),
#               axis.line = element_line(colour = "black"))
# p2 = p2 + scale_color_manual( name="Entropy", labels=c("Prior",  "Ignore biomes", "With biomes"), values=c("black", "dodgerblue", "firebrick3") )
# 
# CairoPDF(file=plot_fn2, height=4, width=6)
# print(p2)
# dev.off()
# 
# print(p)
# #print(p2)
