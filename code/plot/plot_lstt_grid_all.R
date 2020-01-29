#
library(cowplot)
library(ggplot2)
library(grid)
library(scales)
library(Cairo)



args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    base_fn = "out.1.t163.f5.mask_fossil_states"
} else if (length(args)==1) {
    base_fn = args[1]
}


# do you find that transitions into a particular biome state to precede or
# follow transitions into a bg state

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
stoch_bg_biome = data.frame(stringsAsFactors = F)
for (i in 1:length(iterations)) {
    
    # get biome and biogeography stochastic mappings per iteration
    it = iterations[i]
    cat("Stage 1, processing iteration ",it," / ", max(iterations), "\n", sep="")
    sample_bg = stoch_bg[ stoch_bg$iteration==it, ]
    sample_biome = stoch_biome[ stoch_biome$iteration==it, ]
    
    # loop over branches
    for (j in 1:length(branches)) {
        
        # get biome and biogeography stochastic mappings per branch
        nd_idx = branches[j]
        branch_bg = sample_bg[ sample_bg$node_index==nd_idx, ]
        branch_biome = sample_biome[ sample_biome$node_index==nd_idx, ]
        
        # interleave biome and biogeography stochastic mappings
        branch_bg_biome = make_joint_state(branch_bg, branch_biome)
        stoch_bg_biome = rbind(as.data.frame(stoch_bg_biome, stringsAsFactors=F),
                               as.data.frame(branch_bg_biome,stringsAsFactors=F))
    }
}

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

bg_names = c("SE Asia", "E Asia", "Europe", "N Amer", "C Amer", "S Amer")
biome_names = c("Sub/tropical","Lucidophyllous","Cloud","Temperate")
biome_cols = as.vector(biome_colors$color[1:4] )# c("brown2", "darkgreen", "cornflowerblue", "darkblue")
names(biome_cols) = c("1","2","3","4")
area_cols = as.vector( bg_colors$color[1:6] ) #c("magenta", "red", "green", "gold", "blue", "cyan")
names(area_cols) = c("1","2","3","4","5","6")

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
d2$count_norm = d2$count / length(iterations) -  1
d2$Support = as.numeric(d2$Support)

#pbb = list()
ylim = c(0, max(d2$count_norm))
p_all = list()
for (i in 1:6) {
    p_tmp = list()
    for (j in 1:4) {
        titleij = paste( bg_names[i], ",", biome_names[j], sep="")
        dij = d2[ d2$Area==i&d2$Biome==j, ]
        pij = ggplot(data = dij, aes(x=age, y=count_norm ))
        pij = pij + geom_line(colour=bg_colors$color[i], size=2)
        pij = pij + geom_line(colour=biome_colors$color[j], linetype=2, size=2)
        pij = pij + ylim(ylim)
        pij = pij + scale_x_reverse()
        pij = pij + ggtitle(titleij)
        pij = pij + xlab("") + ylab("")
        pij = pij + theme(axis.title.x=element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank(),
                          legend.position = "none")   
        p_tmp[[j]] = pij
    }
    p_all[[i]] = plot_grid( p_tmp[[1]], p_tmp[[2]], p_tmp[[3]], p_tmp[[4]], nrow=1 )
}
pg = plot_grid(p_all[[1]], p_all[[2]], p_all[[3]], p_all[[4]], p_all[[5]], p_all[[6]], ncol=1)
#pg = add_sub(pg, "Age (Ma)", vpadding=grid::unit(0,"lines"),y=5, x=0.5, vjust=4.5)

#pg2 = add_sub(pg, "Age (Ma)", vpadding=grid::unit(0,"lines"),y=5, x=0.5, vjust=4.5)

CairoPDF(  file=paste(plot_fp, "lstt/", base_fn, ".lstt_grid.pdf", sep=""), width=8, height=12)
print(pg)
dev.off()

