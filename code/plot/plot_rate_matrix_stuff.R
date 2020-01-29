#
library(cowplot)
library(ggplot2)
library(grid)
library(scales)
library(Cairo)
library(data.table)
library(igraph)

# do you find that transitions into a particular biome state to precede or
# follow transitions into a bg state

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    base_fn = "out.1.t163.f5"
} else if (length(args)==1) {
    base_fn = args[1]
}

# functions
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

get_compound_state = function(area, biome, n_areas=6, n_biomes=4) {
    return( (n_areas)*(biome-1) + area )
}

get_area_biome = function(s, n_areas=6, n_biomes=4) {
    biome = as.integer( (s-1) / n_areas)
    area = s - biome*n_areas
    return( c(area, biome+1) )
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
plot_fp = paste(fp, "code/plot/fig/matrix/", sep="")
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
n_states = n_areas * n_biomes

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
# index ( (bgxbiome) x (bgxbiome) x time )
bin_width = 1
max_time = 90 #ceiling(max(stoch_bg_biome$x1)) + 1
n_bins = max_time / bin_width
state_bins = array(0, dim=c(n_states, n_states, n_bins))

# 0.0 to 0.5, 0.5 to 1.0, etc
ages = seq(0.0, max_time, by=bin_width)

dat_tx_plot_colnames = c( names(stoch_bg_biome[c(),]), "age", "from_joint_state", "to_joint_state" )
dat_tx_plot = data.frame(matrix(ncol=length(dat_tx_plot_colnames), nrow=0), stringsAsFactors=F)
colnames(dat_tx_plot) = dat_tx_plot_colnames

dat_tx_tmp = data.frame(matrix(ncol=length(dat_tx_plot_colnames), nrow=1e3), stringsAsFactors=F)
colnames(dat_tx_tmp) = dat_tx_plot_colnames

idx_tmp = 1
curr_it = -1

stoch_bg_biome = stoch_bg_biome[ stoch_bg_biome$transition_type == "anagenetic", ]



for (i in 1:nrow(stoch_bg_biome)) {
    
    if (curr_it != stoch_bg_biome$iteration[i]) {
        curr_it = stoch_bg_biome$iteration[i]
        cat("Stage 2, processing iteration ",curr_it," / ", max(stoch_bg_biome$iteration), "\n", sep="")
    }
    
    
    if (stoch_bg_biome$transition_type == "anagenetic") {

        from_bg_idx = get_bg_state( stoch_bg_biome$start_state_1[i] )
        from_biome_idx = get_biome_state( stoch_bg_biome$start_state_2[i] )
        to_bg_idx = get_bg_state( stoch_bg_biome$end_state_1[i] )
        to_biome_idx = get_biome_state( stoch_bg_biome$end_state_2[i] )
        from_state_idx = get_compound_state( area=from_bg_idx, biome=from_biome_idx )
        to_state_idx = get_compound_state( area=to_bg_idx, biome=to_biome_idx )

        time_idx = as.integer( stoch_bg_biome$transition_time[i] )
        
        state_bins[ from_state_idx, to_state_idx, time_idx ] = state_bins[ from_state_idx, to_state_idx, time_idx ] + 1
    }
}

bg_lbl = c("SEAs","EAs","Eur","NAm","CAm","SAm")
biome_lbl = c("Tr","Wm","Cl","Fr")
state_lbl = c()
for (i in 1:n_biomes) {
    for (j in 1:n_areas) {
        state_lbl = c(state_lbl, paste(biome_lbl[i],bg_lbl[j],sep="+"))
    }
}
from_state_lbl = paste("from",state_lbl,sep="_")
to_state_lbl = paste("to",state_lbl,sep="_")

for (i in 1:n_states) { 
    state_bins[i,i,] = 0
}
sum_state_bins = rowSums(state_bins, dims=2)
rownames(sum_state_bins)=from_state_lbl
colnames(sum_state_bins)=to_state_lbl
rownames(sum_state_bins)=state_lbl
colnames(sum_state_bins)=state_lbl


s_ratio = rep(0, n_states)
names(s_ratio)=state_lbl

# posterior number of events into area-biome
r_in = apply( sum_state_bins*(1/length(iterations)), 2, sum )
# posterior number of events leaving area-biome
r_out = apply( sum_state_bins*(1/length(iterations)), 1, sum )
# posterior number of events involving each area-biome
r_total = r_in + r_out

plot(r_total,r_in/r_out, col="blue", border="red")
text(r_total, r_in/r_out, labels = state_lbl)

s_bg_lbl = rep(bg_lbl, 4)
s_biome_lbl = as.vector( sapply( biome_lbl, function(x) { rep(x, 6) }) )

dat = data.frame( count=r_total, flux=r_in-r_out, ratio=r_in/r_out, flux_ratio=(r_in-r_out)/r_total,state=state_lbl, area=s_bg_lbl, biome=s_biome_lbl ) 

s_biome_colors = as.vector(biome_colors$color)[1:4]
names(s_biome_colors) = as.vector(unique(dat$biome))
s_bg_colors = as.vector(bg_colors$color)[1:6]
names(s_bg_colors) = as.vector(unique(dat$area))

pp = ggplot(dat, aes(x=count, y=flux_ratio))
pp = pp + geom_point( aes(colour=biome), shape=21, size=3.5)
pp = pp + geom_point( aes(colour=biome), shape=21, size=4)
pp = pp + geom_point( aes(colour=biome), shape=21, size=4.5)
pp = pp + geom_point( aes(colour=biome), shape=21, size=5)
pp = pp + geom_point( aes(colour=biome), shape=21, size=5.5)
pp = pp + geom_point( aes(colour=biome), shape=21, size=6)
pp = pp + geom_point( aes(fill=area), shape=21, size=3)
pp = pp + scale_colour_manual( values=s_biome_colors )
pp = pp + scale_fill_manual( values=s_bg_colors )
pp = pp + ylim(-1,1)
pp = pp + xlab( "Mean posterior # total events" )
pp = pp + ylab( "Mean posterior relative flux \n(# incoming - # outgoing) / # total" )
#pp = pp + geom_point( aes(colour=area), shape=21, size=4)


plot_fn = paste(plot_fp, base_fn, ".flux.pdf", sep="")
CairoPDF( plot_fn, height=6, width=6 )
print(pp)
dev.off()

