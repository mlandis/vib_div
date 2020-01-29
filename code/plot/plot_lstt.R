#
library(cowplot)
library(ggplot2)
library(grid)
library(scales)
library(Cairo)
library(data.table)



args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    base_fn = "out.1.t163.f5"
   # base_fn = "out.1.t163.f5.mask_fossil_states"
} else if (length(args)==1) {
    base_fn = args[1]
}

# do you find that transitions into a particular biome state to precede or
# follow transitions into a bg state

# What is the value of d for the number of samples?
# SK Ernest 1978 -- compute smallest sample size, n, that satisfies inequality
#   \prod_i Pr( |f_i - \pi_i| < d ) > 1 - \alpha
# where d determines error in obs f_i from true parameter \pi_i
#
# We can also solve for alpha given D and N, and D given N and alpha.

# We compute the posterior mean for each bin/region/biome thing, S_ijt,
# where i is region, j is biome, and t is the time bin.
# We then compute S_it = sum_j S_ijt and S_jt = sum_i S_ijt.
# Some time bins have very few samples. We measure support in two ways using
# a frequentist significance test for the multinomial distribution.
#
# 1) SK Ernest computes the minimum number of samples, $N'$, needed to 
#    guarantee that all 95% of replicates simulated under the true values
#    of the approximated multinomial distribution would produce sampled
#    frequencies within +/- 0.05 units of the true simulating parameters.
#    where N' = 510 regardless of the number of categories.
#
# 2) We use this same theory to solve for what is the alpha significance
#    for the number of actual samples we have -- i.e. how often
#    do we expect to generate sampled proportions +/- d of the true values??
#

# We report values of alpha for each time bin. Alpha values correspond to
# the probability that we would expect proportions with values, f_i, that
# similar to multinomial probabilities, p_i, such that |f_i-p_i| <= 0.05
# provided the sample size is equals what we sampled for each LSTT time bin
# by MCMC.

# We report our support metric as the proportion of samples N_it = 
# and the proportion of th
# 
# 
# The support statistic is defined for the $S_{i,t} = S_{i,j,t}$ samples
# drawn for time bin $t$ for biome/region $i$. The statistic
# reports the probability, $\alpha$, that all category frequencies are
# within the true underlying posterior proportions $\pi_{i,j,t}$
# 

calculate_ske = function(s, k, alpha=0.05, D=0.05) {
    # s     : number of samples
    # k     : number of bins
    # alpha : significance level
    # D     : estimate difference must satisfy | f_i - pi_i | > D for all i in k
    # N     : min num samples needed for
    n = c()
    d = c()
    a = c()

    if (s<=1) {
        return( list(n=Inf, d=1, a=1, alpha=alpha, s=s, D=D) )
    }
    for (m in 1:k) {
        #print("m")
        #print(m)
        # significance value for alpha/m
        z = qnorm(p=(1-alpha/(2*m)),mean=0,sd=1)
        #print("z")
        #print(z)
        
        # how many samples, n, needed for z-values tolerating error D
        n[m] = ceiling(z^2*(1/m) * (1 - 1/m) / D^2)
        #print("n")
        #print(n)
        
        # how accurately can you estimate pi_i for sample size s w/ z-values
        d[m] = sqrt(z^2*(1/m) * (1 - 1/m) / s )
        #print("d")
        #print(d)
        
        # what alpha level can you obtain with s samples for tolerance D
        aval = seq(0,1,by=0.001)
        aval = aval[2:(length(aval)-1)]
        
        # get z-values for possible alpha values
        zz = qnorm(p=(1-aval/(2*m)),mean=0,sd=1)
        #print("zz")
        #print(zz)
        
        # compute n values for possible alpha values under D
        nn = (zz^2*(1/m) * (1 - 1/m) / D^2)
        #print("nn")
        #print(nn)
        
        # for values of n greater than the actual sample, which is the smallest?
        a_idx = min( sum(s < nn) + 1, length(nn) )
        #print("a_idx")
        #print(a_idx)
        
        # choose that alpha value
        a[m] = aval[a_idx]
        #print("a[m]")
        #print(a[m])
        
        #cat("\n\n")
        
    }
    return( list(n=max(n), d=max(d), a=max(a, na.rm=T), alpha=alpha, s=s, D=D) )
}

# functions
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


get_bg_state = function(s) {
    if (s==1)       return(c(1))
    else if (s==2)  return(c(2))
    else if (s==3)  return(c(3))
    else if (s==4)  return(c(4))
    else if (s==5)  return(c(5))
    else if (s==6)  return(c(6))
    else if (s==7)  return(c(1,2))
    else if (s==8)  return(c(1,3))
    else if (s==9)  return(c(2,3))
    else if (s==10) return(c(1,4))
    else if (s==11) return(c(2,4))
    else if (s==12) return(c(3,4))
    else if (s==13) return(c(1,5))
    else if (s==14) return(c(2,5))
    else if (s==15) return(c(3,5))
    else if (s==16) return(c(4,5))
    else if (s==17) return(c(1,6))
    else if (s==18) return(c(2,6))
    else if (s==19) return(c(3,6))
    else if (s==20) return(c(4,6))
    else if (s==21) return(c(5,6))
}

get_biome_state = function(s) {
    return(s+1)
}

make_joint_state = function(b1, b2) {
    
    #cat("\n\n")
    # create joint events from independent events
    b1_tmp = cbind(b1, start_state_2=NA, end_state_2=NA, transition_char=1)
    names( b1_tmp )[ match(c("start_state", "end_state"), names(b1_tmp) ) ] = c("start_state_1","end_state_1")
    b2_tmp = cbind(b2, start_state_1=NA, end_state_1=NA, transition_char=2)
    names( b2_tmp )[ match(c("start_state", "end_state"), names(b2_tmp) ) ] = c("start_state_2","end_state_2")
    b = rbind( as.data.frame(b1_tmp, stringsAsFactors=F),
               as.data.frame(b2_tmp, stringsAsFactors=F))
    
    # sort joint events by transition time
    b = b[order(b$transition_time, decreasing=T),]

    # get the first valid start state for each character
    start_state_1 = na.omit(b$start_state_1)[1]
    start_state_2 = na.omit(b$start_state_2)[1]
    
    # there must be at least two events; drop first event from list if it's no_change
    if (b$transition_type[1]=="no_change") {
        b = b[-c(1),]    
    }
    
    #print(b)
    
    # update missing states for the first event
    if (b$transition_char[1]==1) {
        b[1, c("start_state_2", "end_state_2")] = c( start_state_2, start_state_2 )            
    } else if (b$transition_char[1]==2) {
        b[1, c("start_state_1", "end_state_1")] = c( start_state_1, start_state_1 )
    }
    
    #cat("\n")
    #print(b)
    #cat("\n")
    
    # if the next remaining event is the last event, return
    if (nrow(b) > 1) {
        
        # if there's more than one event, propagate joint states
        for (i in 2:nrow(b)) {
            
            # update end states for next iteration
            end_state_1 = b$end_state_1[i-1]
            end_state_2 = b$end_state_2[i-1]
            
            # does the transition occur on character 1 or 2?
            char_idx = b$transition_char[i]
            if (char_idx==1) {
                b[i,c("start_state_2","end_state_2")] = c( end_state_2, end_state_2)            
            } else if (char_idx==2) {
                b[i,c("start_state_1","end_state_1")] = c( end_state_1, end_state_1)
            }
    
            
        }
    }
    
    if (nrow(b)==1) {
        x1 = b$branch_start_time
        x2 = b$branch_end_time
        #b = cbind(b, x1=b$branch_start_time, x2=b$branch_end_time, s1$start_state)
    } else {
        x1 = c( b$branch_start_time[1], b$transition_time )
        x2 = c( b$transition_time, b$branch_end_time[1] )
        b = rbind(as.data.frame(b, stringsAsFactors=F), 
                  as.data.frame(b[nrow(b),], stringsAsFactors=F))
        b[nrow(b),c("start_state_1","start_state_2","transition_time")] = b[nrow(b),c("end_state_1","end_state_2","branch_end_time")]
    }
    b = cbind(b, x1=x1, x2=x2, s1=b$start_state_1, s2=b$start_state_2)
    #print("after")
    #print(b)
    return(b)
    #return(list(b=b,b_orig=b_orig))
    
}


fp = "/Users/mlandis/projects/vib_div/"
out_fp = paste(fp, "output/", sep="")
plot_fp = paste(fp, "code/plot/fig/", sep="")
#out_str = paste( base_fn, ".biome", sep="" )
#tree_fn = paste(out_fp, out_str, ".ase.tre", sep="")

col_bg_fn = paste(fp, "code/plot/range_colors.n6.txt",sep="")
col_biome_fn = paste(fp, "code/plot/biome_colors.n4.txt",sep="")
bg_colors = read.csv(col_bg_fn, header=T)
biome_colors = read.csv(col_biome_fn, header=T)

# files
bg_fn = paste(out_fp, base_fn, ".bg.history.tsv", sep="")
biome_fn = paste(out_fp, base_fn, ".biome.history.tsv", sep="")

stoch_bg = read.csv(bg_fn, sep="\t", stringsAsFactors=F)
stoch_bg = stoch_bg[ stoch_bg$transition_type != "cladogenetic", ]
stoch_biome = read.csv(biome_fn, sep="\t", stringsAsFactors=F)

stoch_bg$transition_time[ stoch_bg$transition_type=="no_change" ] = stoch_bg$branch_start_time[ stoch_bg$transition_type=="no_change" ]
stoch_biome$transition_time[ stoch_biome$transition_type=="no_change" ] = stoch_biome$branch_start_time[ stoch_biome$transition_type=="no_change" ]

# data dimensions
n_areas = 6
n_biomes = 4

# iterations
f_burn = 0.0
iterations = unique(stoch_bg$iteration)
n_burn = max(1, f_burn*length(iterations))
thinby = 1
iterations = iterations[n_burn:length(iterations)]

iterations = iterations[ seq(1, length(iterations), by=thinby) ]
#iterations = iterations[1:10]
branches = 1:max(unique(stoch_bg$parent_index), na.rm=T)
#branches = c(260)

# loop over iterations
stoch_list = list()
#stoch_bg_biome = data.frame(stringsAsFactors = F)
for (i in 1:length(iterations)) {
    
    # get biome and biogeography stochastic mappings per iteration
    it = iterations[i]
    cat("Stage 1, processing iteration ",it," / ", max(iterations), "\n", sep="")
    sample_bg = stoch_bg[ stoch_bg$iteration==it, ]
    sample_biome = stoch_biome[ stoch_biome$iteration==it, ]
    
    # loop over branches
    tmp_branch_list = list()
    for (j in 1:length(branches)) {
        
        # get biome and biogeography stochastic mappings per branch
        nd_idx = branches[j]
        branch_bg = sample_bg[ sample_bg$node_index==nd_idx, ]
        branch_biome = sample_biome[ sample_biome$node_index==nd_idx, ]
        
        # interleave biome and biogeography stochastic mappings
        tmp_branch_list[[j]] = as.data.frame( make_joint_state(branch_bg, branch_biome), stringsAsFactors=F )
        #stoch_bg_biome = rbind(as.data.frame(stoch_bg_biome, stringsAsFactors=F),
        #                       as.data.frame(branch_bg_biome,stringsAsFactors=F))
    }
    stoch_list[[i]] = rbindlist(tmp_branch_list) 
}
stoch_bg_biome = rbindlist(stoch_list)
#########

# bins
# index ( bg x biome x time )
bin_width = 1
max_time = 90 #ceiling(max(stoch_bg_biome$x1)) + 1
n_bins = max_time / bin_width
state_bins = array(0, dim=c(n_areas, n_biomes, n_bins))

# 0.0 to 0.5, 0.5 to 1.0, etc
ages = seq(0.0, max_time, by=bin_width)

dat_plot_colnames = c( names(stoch_bg_biome[c(),]), "age", "joint_state" )
dat_plot = data.frame(matrix(ncol=length(dat_plot_colnames), nrow=0), stringsAsFactors=F)
colnames(dat_plot) = dat_plot_colnames

dat_tmp = data.frame(matrix(ncol=length(dat_plot_colnames), nrow=1e3), stringsAsFactors=F)
colnames(dat_tmp) = dat_plot_colnames

idx_tmp = 1
curr_it = -1
for (i in 1:nrow(stoch_bg_biome)) {
    
#for (i in 1:2500) {
    if (curr_it != stoch_bg_biome$iteration[i]) {
        curr_it = stoch_bg_biome$iteration[i]
        cat("Stage 2, processing iteration ",curr_it," / ", max(stoch_bg_biome$iteration), "\n", sep="")
    }
    
    bg_idx = get_bg_state( stoch_bg_biome$start_state_1[i] )
    biome_idx = get_biome_state( stoch_bg_biome$start_state_2[i] )
    
    #age_bins = seq(floor(stoch_bg_biome$x2[i]), ceiling(stoch_bg_biome$x1[i]), by=bin_width ) / bin_width
    start_age = floor(stoch_bg_biome$x2[i])
    end_age = ceiling(stoch_bg_biome$x1[i])
    age_bins = start_age:end_age * bin_width
    time_idx = age_bins + 1
    
    #print(i)
    for (j in 1:length(age_bins)) {
        for (k in 1:length(bg_idx)) {
            joint_state = paste(bg_idx[k],biome_idx,sep="_")
            #print(joint_state)
            #dat_plot = rbind(dat_plot, c( stoch_bg[i,], age=age_bins[j], joint_state=joint_state))
            dat_tmp[idx_tmp,] = c( stoch_bg[i,], age=age_bins[j], joint_state=joint_state)
            
            if (idx_tmp == nrow(dat_tmp)) {
                dat_plot = rbind(dat_plot, dat_tmp)
                idx_tmp = 1
            } else if (idx_tmp < nrow(stoch_bg_biome)) {
                idx_tmp = idx_tmp + 1
            }
        }
    }
    #print(time_idx)
    #time_idx = n_bins - seq(floor(stoch_bg_biome$x2[i]), ceiling(stoch_bg_biome$x1[i]), by=bin_width ) / bin_width + 1
    state_bins[ bg_idx, biome_idx, time_idx ] = state_bins[ bg_idx, biome_idx, time_idx ] + 1
}

#stop()

dat_plot = rbind(dat_plot, dat_tmp[1:idx_tmp,])

bg_names = c("SEAs", "EAs", "Eur", "NAm", "CAm", "SAm")
biome_names = c("Trop.","Warm","Cloud","Cold")
biome_cols = as.vector(biome_colors$color[1:4] )# c("brown2", "darkgreen", "cornflowerblue", "darkblue")
names(biome_cols) = c("1","2","3","4")
area_cols = as.vector( bg_colors$color[1:6] ) #c("magenta", "red", "green", "gold", "blue", "cyan")
names(area_cols) = c("1","2","3","4","5","6")

bg_label_col = c("white","white","white","white","black","black")
biome_label_col = c("white","white","black","white")
dat_plot_2 = matrix(nrow=0, ncol=6)
colnames(dat_plot_2) = c("age","count","Area","Biome","sj","Support")

for (i in 1:dim(state_bins)[1]) {
    for (j in 1:dim(state_bins)[2]) {
        for (k in 1:dim(state_bins)[3]) {
            dat_plot_2 = rbind(dat_plot_2, c( ages[k], state_bins[i,j,k], i, j, paste(i, j, sep="_"), 0))
        }
    }
}

d2 = data.frame(dat_plot_2,  stringsAsFactors=FALSE)
d2$age = as.numeric(d2$age)
d2$count = as.numeric(d2$count)
d2$Support = as.numeric(d2$Support)


biome_conf = t(apply( state_bins, 3, colSums))
bg_conf = t(apply( state_bins, 3, rowSums))
for (i in 1:n_bins) {
    for (j in 1:n_biomes) {
        if (biome_conf[i,j] > 510) {
            biome_conf[i,j] = 1
        } else {
            biome_conf[i,j] = 0
        }
        #biome_conf[i,j] = 1 - calculate_ske( biome_conf[i,j], 4 )$a
    }
    for (j in 1:n_areas) {
        #bg_conf[i,j] = 1 - calculate_ske( bg_conf[i,j], 6 )$a
        if (bg_conf[i,j] > 510) {
            bg_conf[i,j] = 1
        } else {
            bg_conf[i,j] = 0
        }
    }
}
# n_min = 510 # from SK Ernest
# biome_conf = biome_conf / n_min
# biome_conf[ biome_conf > 1 ] = 1
# biome_conf[ biome_conf < 0.01 ] = 0
# bg_conf = bg_conf / n_min
# bg_conf[ bg_conf > 1 ] = 1
# bg_conf[ bg_conf < 0.01 ] = 0

d2_ages = unique(d2$age)

#trunc_val = floor( length(iterations) / 20 )
d2_bg_trunc = d2
d2_biome_trunc = d2

for (i in 1:length(d2_ages)) {
    for (j in 1:n_areas) {
        c_ij = d2_bg_trunc[ d2_bg_trunc$age==i & d2_bg_trunc$Area==j, ]$count
        if (length(c_ij) == 0) { 
            # do nothing
        } else {
            d2_bg_trunc[ d2_bg_trunc$age==i & d2_bg_trunc$Area==j, ]$Support = bg_conf[i,j]
        }
    }
    
    for (j in 1:n_biomes) {
        
        c_ij = d2_biome_trunc[ d2_biome_trunc$age==i & d2_biome_trunc$Biome==j, ]$count
        # #cat(i,j,sum(c_ij),"\n")
        if (length(c_ij) == 0) { 
            # do nothing
        } else {
            d2_biome_trunc[ d2_biome_trunc$age==i & d2_biome_trunc$Biome==j, ]$Support = biome_conf[i,j]
        }
    }
}


sz=4
pbg=list()
for (i in 1:6) {
    p2 = ggplot(data = d2_bg_trunc[d2_bg_trunc$Area==i&d2_bg_trunc$Support>0,], aes(fill=Biome, x=age, y=count)) #, alpha=Support)) 
    title_i = bg_names[i]
    p2 = p2 + geom_bar(stat="identity",position="fill")
    #p2 = p2 + scale_y_continuous( limits=c(0,1.0) ) #labels=percent_format())
    p2 = p2 + scale_x_continuous("", trans="reverse", limits=c(max_time,0)) #, sec.axis = sec_axis(~ ., name=title_i, breaks=c() ))
    #p2 = p2 + ggtitle(title_i)
    p2 = p2 + scale_fill_manual(values=biome_cols, labels = c("Trop.", "Warm", "Cloud", "Cold"), guide="legend")
    p2 = p2 + scale_alpha_continuous(range=c(0,1), guide="legend")
    p2 = p2 + annotate("text", x=max_time/2, y=0.8, label=title_i, colour=bg_label_col[i], size=sz)
    #p2 = p2 + guides(colour = guide_legend(order = 99)) #, alpha = guide_legend(order = 4))
    
    if (i==1) {
        p2_bg = p2 + guides(alpha = F, fill=guide_legend( title.theme=element_text(size=14),
                                                          label.theme = element_text(size = 12)))
        p_bg_legend = get_legend(p2_bg)
    }
    p2 = p2 + theme(axis.title.x=element_blank(),
                    axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.text.x=element_text(size=12),
                    #axis.ticks.y=element_blank(),
                    legend.text=element_text(size=12),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(colour = "black"),
                    legend.position = "none")   
    #p2 = p2 + annotate("text", x=max_time/2, y=Inf, label = title_i, hjust = 0.5, vjust = 1, size=3)
    
    pbg[[i]] = p2
}

pbg_1 = plot_grid( pbg[[ 1]], pbg[[ 2]], pbg[[ 3]], pbg[[ 4]], pbg[[ 5]], pbg[[6]], ncol=1) 
pbg_2 = add_sub(pbg_1, "Age (Ma)", vpadding=grid::unit(0,"lines"),y=5, x=0.5, vjust=4.5)
#ggdraw(pbg_2)



pbiome=list()
for (i in 1:4) {
    p2 = ggplot(data = d2_biome_trunc[d2_biome_trunc$Biome==i&d2_biome_trunc$Support>0,], aes(fill=Area, x=age, y=count)) #, alpha=Support)) 
    title_i = biome_names[i]
    p2 = p2 + geom_bar(stat="identity",position="fill")
    #p2 = p2 + scale_y_continuous( limits=c(0,1.25) )#labels=percent_format())
    p2 = p2 + scale_x_continuous("", trans="reverse", limits=c(max_time,0) )
    #p2 = p2 + ggtitle(title_i)
    p2 = p2 + scale_fill_manual(values=area_cols, labels = c("SEAs", "EAs", "Eur", "NAm", "CAm", "SAm"), guide="legend")
    p2 = p2 + scale_alpha_continuous(range = c(0, 1), guide="legend") #, breaks=seq(0,1,by=.5), limits=seq(0,1,by=0.5))
    p2 = p2 + annotate("text", x=max_time/2, y=0.8, label=title_i, colour=biome_label_col[i], size=sz)
   # p2 = p2 + guides(colour = guide_legend(order = 1), guide_legend(order = 99))
    
    if (i == 1) {
        px = p2
        
        #p2_biome = p2 + guides(alpha = F)
        p2_biome = p2 + guides(alpha = F, fill=guide_legend( title.theme=element_text(size=14),
                                                             label.theme = element_text(size=12)))
        p_biome_legend = get_legend(p2_biome)
        
        p2_freq = px + scale_fill_manual(values=area_cols, labels = c("SEAs", "EAs", "Eur", "NAm", "LAm"), guide=FALSE)
        p_freq_legend = get_legend(p2_freq)
        
        #p2_freq = p2 + guides(colour = F)
        
    }
    
    p2 = p2 + theme(axis.title.x=element_blank(),
                    axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.text.x=element_text(size=12),
                    #axis.ticks.y=element_blank(),
                    legend.text=element_text(size=18),
                    legend.title=element_text(size=18),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(colour = "black"),
                    legend.position = "none")      
    #p2 = p2 + annotate("text", x=max_time/2, y=Inf, label = title_i, hjust = 0.5, vjust = 1, size=3)
    
    pbiome[[i]] = p2
}



plegend_left = plot_grid(NULL, p_bg_legend, NULL, ncol=1)
plegend_right = plot_grid(NULL, p_biome_legend, NULL, ncol=1)

    
pbiome_1 = plot_grid( NULL, pbiome[[ 1]], pbiome[[ 2]], pbiome[[ 3]], pbiome[[ 4]], NULL,  ncol=1)
pbiome_2 = add_sub(pbiome_1, "Age (Ma)", vpadding=grid::unit(0,"lines"),y=5, x=0.5, vjust=4.5)
#ggdraw(pbiome_2)



#pboth = plot_grid( pbg_2, p_bg_legend, pbiome_2, p_biome_legend, p_freq_legend, ncol=5, rel_heights=c(5,1,4,1,1), rel_widths=c(3,1,3,1,1) )
pboth = plot_grid( plegend_left, pbg_2, pbiome_2, plegend_right, ncol=4, rel_heights=c(1,6,4,2), rel_widths=c(1,3,3,1) )
ggdraw(pboth)

#ggsave( filename=paste(plot_fp, "lstt/", base_fn, ".lstt.pdf", sep=""), plot=pboth, device="pdf", width=7, height=7)

CairoPDF(  file=paste(plot_fp, "lstt/", base_fn, ".lstt.pdf", sep=""), width=7, height=7)
print(pboth)
dev.off()
