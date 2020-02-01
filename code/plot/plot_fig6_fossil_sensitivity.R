library(ggplot2)
library(Cairo)

source("vib_div_util.R")

fp = "../../"
out_fp  = paste0(fp, "output/sensitivity_test", sep="")
plot_fp = paste0(fp, "code/plot/fig/", sep="")
plot_fn = paste0(plot_fp, "fig6_cold_root_sensitivity.pdf", sep="")

# collect all files
files = list.files(out_fp, full.names=T)
files = files[ grep("biome.states.txt", files) ]
files = files[ grep("pi_cold_", files) ]
files_nf = files[ grep("mask", files) ]
files_wf = files[ !(files %in% files_nf) ]
n_files = length(files_nf)
f_burn = 0.5

# collect results
df = data.frame( )
for (i in 1:n_files) {
    fp_toks = strsplit(files_nf[i], split="/")[[1]]
    fn_toks = strsplit(fp_toks[length(fp_toks)], split="\\.")[[1]]
    pi_toks = as.numeric(strsplit(fn_toks[6], split="_")[[1]][3])
    
    dat_nf = read.table(files_nf[i], sep="\t", header=T)
    dat_wf = read.table(files_wf[i], sep="\t", header=T)
    
    dat_nf = dat_nf[(f_burn*nrow(dat_nf)):nrow(dat_nf),]
    dat_wf = dat_wf[(f_burn*nrow(dat_wf)):nrow(dat_wf),]
    
    pi_temp_nf = pi_temp_prob(dat_nf)
    pi_temp_wf = pi_temp_prob(dat_wf)
    h_wf = compute_entropy(dat_wf)
    h_nf = compute_entropy(dat_nf)
    h_pi = compute_entropy_exp(pi_toks/100)
    df = rbind(df, c(pi_toks/100, pi_temp_nf, pi_temp_wf, h_nf, h_wf, h_pi))
}
colnames(df) = c("pi_freeze", "no_fossils", "with_fossils", "entropy_nf", "entropy_wf", "entropy_pi")
df = df[ order(df$pi_freeze), ]

# bookend df with Prob==0 and Prob==1
df = rbind(
        c(0, 0, 0, NA, NA, NA),
        df,
        c(1, 1, 1, NA, NA, NA))

# plotting
p = ggplot(df, aes(x=pi_freeze, y=no_fossils))

# define functions and inverse functions to get critical probs
f_wf_inv = approxfun( df$with_fossils, df$pi_freeze )
f_wf     = approxfun( df$pi_freeze, df$with_fossils )
f_nf_inv = approxfun( df$no_fossils, df$pi_freeze )
f_nf     = approxfun( df$pi_freeze, df$no_fossils )

# plot lines for critical probabilities
p = p + geom_segment( mapping=aes(x=0.25, xend=0.25, y=-Inf, yend=f_wf(0.25)), lty=2, col='gray')
p = p + geom_segment( mapping=aes(x=-Inf, xend=0.25, y=f_wf(0.25), yend=f_wf(0.25)), lty=2, col='gray')
p_flip = f_wf_inv(0.5)
p = p + geom_segment( mapping=aes(x=p_flip, xend=p_flip, y=-Inf, yend=f_wf(p_flip)), lty=2, col='gray')
p = p + geom_segment( mapping=aes(x=-Inf, xend=p_flip, y=f_wf(p_flip), yend=f_wf(p_flip)), lty=2, col='gray')
p_sig = f_wf_inv(0.95)
p = p + geom_segment( mapping=aes(x=p_sig, xend=p_sig, y=-Inf, yend=f_wf(p_sig)), lty=2, col='gray')
p = p + geom_segment( mapping=aes(x=-Inf, xend=p_sig, y=f_wf(p_sig), yend=f_wf(p_sig)), lty=2, col='gray')

# annotate critical probabilities
p = p + annotate(geom="text", x=0.25-0.03, y=f_wf(0.25)+0.05, label=expression( paste(p[fair])))
p = p + annotate(geom="text", x=p_flip-0.03, y=f_wf(p_flip)+0.05, label=expression( paste(p[flip])))
p = p + annotate(geom="text", x=p_sig+0.05, y=f_wf(p_sig)-0.04, label=expression( paste(p[sig])))
p = p + geom_point( mapping=aes(x=0.25, y=f_wf(0.25)), col="firebrick3")
p = p + geom_point( mapping=aes(x=p_flip, y=f_wf(p_flip)), col="firebrick3")
p = p + geom_point( mapping=aes(x=p_sig, y=f_wf(p_sig)), col="firebrick3")

# add main results
p = p + geom_abline(  intercept=0, slope=1, colour="lightgray", lwd=0.25)
p = p + geom_line(data=df, mapping=aes(x=pi_freeze, y=no_fossils, colour="black"), alpha=1)
p = p + geom_line(data=df, mapping=aes(x=pi_freeze, y=with_fossils, colour="red"), alpha=1)
p = p + xlab( expression( paste(Prob,"( ",X[root]," = ",Cold," | ",X[miss]," )") ) )
p = p + ylab( expression( paste(Prob,"( ",X[root]," = ",Cold," | ",X[obs],", ",X[miss]," )") ) )

# set theme and legend
p = p + theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"))
p = p + scale_color_manual( name="Dataset", labels=c("Masked", "Complete"), values=c("dodgerblue", "firebrick3") )

# print
CairoPDF(file=plot_fn, height=4, width=6)
print(p)
dev.off()
