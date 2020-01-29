#
library(cowplot)
library(ggplot2)
library(grid)
library(scales)
library(Cairo)
library(data.table)
library(igraph)
library(ape)

source("/Users/mlandis/projects/reml_levy/code/ggbiplot_custom.r")
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

bins_to_probs = function(x, lbl) {
    n = dim(x)[3]
    y = rowSums(x, dims=2)
    rownames(y) = lbl
    colnames(y) = lbl
    y = y * (1/n)
    return(y)
}

max_bins_to_probs = function(x, lbl) {
    n = dim(x)[3]
    y = rep(0, length(lbl))
    names(y) = lbl
    for (i in 1:n) {
        y = y + apply(x[,,i], 2, max)
    }
    y = y * (1/n)
    return(y)
}


get_shift_idx = function(x1, x2) {
    if (length(x1)==0 || length(x2)==0) {
        return( list(NA,NA) )
    } else if (length(x1)==1 && length(x2)==1) {
        return( list(x1, x2) )
    } else if (length(x1)==1 && length(x2)==2) {
        return( list(x1, x2[ !x1==x2 ] ) )
    } else if (length(x1)==2 && length(x2)==2) {
        return( list(x1, x2) )
    } else if (length(x1)==2 && length(x2)==1) {
        return( list(NA,NA) )
    }
}

max_age = function(x, type) {
    x = x[ x$transition_type=="anagenetic", ]
    if (type=="biome") {
        y = matrix(0, nrow=4, ncol=4)
    } else if (type=="area") {
        y = matrix(0, nrow=6, ncol=6)
    } else if (type=="state") {
        y = matrix(0, nrow=24, ncol=24)
    }
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
   
    # update missing states for the first event
    if (b$transition_char[1]==1) {
        b[1, c("start_state_2", "end_state_2")] = c( start_state_2, start_state_2 )            
    } else if (b$transition_char[1]==2) {
        b[1, c("start_state_1", "end_state_1")] = c( start_state_1, start_state_1 )
    }
    
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
    b = cbind(b, x1=x1, x2=x2, dx=x1-x2, s1=b$start_state_1, s2=b$start_state_2)

    return(b)
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
phy_fn = paste(out_fp, base_fn, ".tre", sep="")

stoch_bg = read.csv(bg_fn, sep="\t", stringsAsFactors=F)
stoch_bg = stoch_bg[ stoch_bg$transition_type != "cladogenetic", ]
#stoch_bg = stoch_bg[ stoch_bg$start_state!=0, ]
stoch_biome = read.csv(biome_fn, sep="\t", stringsAsFactors=F)
phy_dat=read.table(phy_fn, sep="\t", stringsAsFactors=F, header=T)

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
        if (length(branch_bg$start_state)==1 && !any(branch_bg$start_state[1]==0)) {
            tmp_branch_list[[j]] = as.data.frame( make_joint_state(branch_bg, branch_biome), stringsAsFactors=F )
        }
        #stoch_bg_biome = rbind(as.data.frame(stoch_bg_biome, stringsAsFactors=F),
        #                       as.data.frame(branch_bg_biome,stringsAsFactors=F))
    }
    stoch_list[[i]] = rbindlist(tmp_branch_list) 
}
stoch_bg_biome = rbindlist(stoch_list)


#########
#########
#########


# bins
# index ( (bgxbiome) x (bgxbiome) x time )
bin_width = 1
max_time = 90 #ceiling(max(stoch_bg_biome$x1)) + 1
n_bins = max_time / bin_width
n_iter = length(iterations)

state_freq = array(0, dim=c(n_states, n_states, n_iter))

state_bins = array(0, dim=c(n_states, n_states, n_iter))
state_biome_bins = array(0, dim=c(n_biomes, n_biomes, n_iter))
state_area_bins = array(0, dim=c(n_areas, n_areas, n_iter))

time_state_bins = array(0, dim=c(n_states, n_states, n_iter))
time_biome_bins = array(0, dim=c(n_biomes, n_biomes, n_iter))
time_area_bins = array(0, dim=c(n_areas, n_areas, n_iter))

age_state_bins = array(0, dim=c(n_states, n_states, n_iter))
age_biome_bins = array(0, dim=c(n_biomes, n_biomes, n_iter))
age_area_bins = array(0, dim=c(n_areas, n_areas, n_iter))

max_age_state_bins = array(0, dim=c(n_states, n_states, n_iter))
max_age_biome_bins = array(0, dim=c(n_biomes, n_biomes, n_iter))
max_age_area_bins = array(0, dim=c(n_areas, n_areas, n_iter))


# 0.0 to 0.5, 0.5 to 1.0, etc
#ages = seq(0.0, max_time, by=bin_width)

#dat_tx_plot_colnames = c( names(stoch_bg_biome[c(),]), "age", "from_joint_state", "to_joint_state" )
#dat_tx_plot = data.frame(matrix(ncol=length(dat_tx_plot_colnames), nrow=0), stringsAsFactors=F)
#colnames(dat_tx_plot) = dat_tx_plot_colnames

#dat_tx_tmp = data.frame(matrix(ncol=length(dat_tx_plot_colnames), nrow=1e3), stringsAsFactors=F)
#colnames(dat_tx_tmp) = dat_tx_plot_colnames

idx_tmp = 1
curr_it = -1

#stoch_bg_biome = stoch_bg_biome[ stoch_bg_biome$transition_type == "anagenetic", ]

#for (i in 1:500) { #nrow(stoch_bg_biome)) {
#n_iter = 20
for (k in 1:n_iter) {

    # get iteration
    it = iterations[k]
    cat("Stage 2, state times/counts ", it," / ", max(stoch_bg_biome$iteration), "\n", sep="")
    
    # get data for iteration
    dat_it = stoch_bg_biome[ stoch_bg_biome$iteration == it, ]
    
    # total tree length
    tree_length = sum(dat_it$dx)
      
    max_age_biome = matrix(0, nrow=4, ncol=4)
    max_age_area  = matrix(0, nrow=6, ncol=6)
    max_age_state = matrix(0, nrow=24, ncol=24)
    for (j in 1:nrow(dat_it)) {
    
        # get info for event row in iteration
        from_bg_idx     = get_bg_state( dat_it$start_state_1[j] )
        from_biome_idx  = get_biome_state( dat_it$start_state_2[j] )
        to_bg_idx       = get_bg_state( dat_it$end_state_1[j] )
        to_biome_idx    = get_biome_state( dat_it$end_state_2[j] )
        from_state_idx  = get_compound_state( area=from_bg_idx, biome=from_biome_idx )
        to_state_idx    = get_compound_state( area=to_bg_idx, biome=to_biome_idx )
        age             = dat_it$transition_time[j]
        time_idx        = as.integer( age )
        dt              = dat_it$dx[j]
        
        # anagenetic event counts
        event_type     = dat_it$transition_type[j]
        is_root        = is.na(dat_it$parent_index)
        
        # get shift index        
        bg_shift_idx    = get_shift_idx( from_bg_idx, to_bg_idx )
        state_shift_idx = get_shift_idx( from_state_idx, to_state_idx )
        biome_shift_idx = get_shift_idx( from_biome_idx, to_biome_idx )
        shift_from_bg_idx    = bg_shift_idx[[1]];    shift_to_bg_idx    = bg_shift_idx[[2]]
        shift_from_state_idx = state_shift_idx[[1]]; shift_to_state_idx = state_shift_idx[[2]]
        shift_from_biome_idx = biome_shift_idx[[1]]; shift_to_biome_idx = biome_shift_idx[[2]]
        
        # time spent in state/biome/area
        time_biome_bins[ from_biome_idx, to_biome_idx, k ] = time_biome_bins[ from_biome_idx, to_biome_idx, k ]
        time_area_bins[ from_bg_idx, to_bg_idx, k ] = time_area_bins[ from_bg_idx, to_bg_idx, k ] + dt
        time_state_bins[ from_state_idx, to_state_idx, k ] = time_state_bins[ from_state_idx, to_state_idx, k ] + dt

        
        # count of shifts into state/biome/areas
        if (event_type == "anagenetic") {

            # did a particular type of shift happen in this sample?
            state_freq[ from_state_idx, to_state_idx, k ] = 1
            
            if (shift_from_state_idx != shift_to_state_idx && !is.na(shift_from_state_idx) && !is.na(shift_to_state_idx)) {
                state_bins[ shift_from_state_idx, shift_to_state_idx, k ] = state_bins[ shift_from_state_idx, shift_to_state_idx, k ] + 1
                age_state_bins[ shift_from_state_idx, shift_to_state_idx, k ] = age_state_bins[ shift_from_state_idx, shift_to_state_idx, k ] + age
                max_age_state[ shift_from_state_idx, shift_to_state_idx ] = max( c(age, max_age_state[ shift_from_state_idx, shift_to_state_idx ]) )            
            }
            if (shift_from_biome_idx != shift_to_biome_idx && !is.na(shift_from_biome_idx) && !is.na(shift_to_biome_idx)) {
                state_biome_bins[ shift_from_biome_idx, shift_to_biome_idx, k ] = state_biome_bins[ shift_from_biome_idx, shift_to_biome_idx, k] + 1
                age_biome_bins[ shift_from_biome_idx, shift_to_biome_idx, k ] = age_biome_bins[ shift_from_biome_idx, shift_to_biome_idx, k ] + age
                max_age_biome[ shift_from_biome_idx, shift_to_biome_idx ] = max( c(age, max_age_biome[ shift_from_biome_idx, shift_to_biome_idx ]) )
            }
            if (shift_from_bg_idx != shift_to_bg_idx && !is.na(shift_from_bg_idx) && !is.na(shift_to_bg_idx)) {
                state_area_bins[ shift_from_bg_idx, shift_to_bg_idx, k ] = state_area_bins[ shift_from_bg_idx, shift_to_bg_idx, k ] + 1
                age_area_bins[ shift_from_bg_idx, shift_to_bg_idx, k ] = age_area_bins[ shift_from_bg_idx, shift_to_bg_idx, k ] + age
                max_age_area[ shift_from_bg_idx, shift_to_bg_idx ] = max( c(age, max_age_area[ shift_from_bg_idx, shift_to_bg_idx ]) )
            }
            
        }
    }
    
    #readline(prompt="Press [enter] to continue")
    max_age_biome_bins[,,k] = max_age_biome
    max_age_area_bins[,,k] = max_age_area
    max_age_state_bins[,,k] = max_age_state
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

# ignore self-transitions for counts (but include for times) in history tsv
for (i in 1:n_states) { 
    state_bins[i,i,] = 0
}

sum_state_bins      = bins_to_probs( state_bins, state_lbl )
sum_area_bins       = bins_to_probs( state_area_bins, bg_lbl )
sum_biome_bins      = bins_to_probs( state_biome_bins, biome_lbl )
sum_time_state_bins = bins_to_probs( time_state_bins, state_lbl )
sum_time_biome_bins = bins_to_probs( time_biome_bins, biome_lbl )
sum_time_area_bins  = bins_to_probs( time_area_bins, bg_lbl )
sum_age_state_bins  = bins_to_probs( age_state_bins, state_lbl )
sum_age_biome_bins  = bins_to_probs( age_biome_bins, biome_lbl )
sum_age_area_bins   = bins_to_probs( age_area_bins, bg_lbl )
sum_max_age_state_bins = max_bins_to_probs( max_age_state_bins, state_lbl )
sum_max_age_biome_bins = max_bins_to_probs( max_age_biome_bins, biome_lbl )
sum_max_age_area_bins  = max_bins_to_probs( max_age_area_bins, bg_lbl )

sum_freq = bins_to_probs( state_freq, state_lbl )
for (i in 1:n_states) {
    for (j in 1:n_states) {
        if (sum_freq[i,j] > 0.95) {
            sum_freq[i,j] = 1
        } else {
            sum_freq[i,j] = 0
        }
    }
}
sum_in_freq = colSums(sum_freq)
sum_in_freq[ sum_in_freq > 1 ] = 1

sum_time_state_bins = sum_time_state_bins # * (sum_freq)
sum_state_bins = sum_state_bins# * sum_freq


s_ratio = rep(0, n_states)
names(s_ratio)=state_lbl

# posterior number of events into area-biome

n_state_in = colSums(sum_state_bins, na.rm=T)
n_state_out = rowSums(sum_state_bins, na.rm=T)
l_state_in = colSums(sum_time_state_bins, na.rm=T)
l_state_out = rowSums(sum_time_state_bins, na.rm=T)
a_state_in = colSums(sum_age_state_bins, na.rm=T)
a_state_out = rowSums(sum_age_state_bins, na.rm=T)
#m_state_in = colSums(sum_max_age_state_bins, na.rm=T)
m_state_in = sum_max_age_state_bins

n_biome_in = colSums(sum_biome_bins, na.rm=T)
l_biome_in = colSums(sum_time_biome_bins, na.rm=T)
a_biome_in = colSums(sum_age_biome_bins, na.rm=T)
#m_biome_in = colSums(sum_max_age_biome_bins, na.rm=T)
m_biome_in = sum_max_age_biome_bins

n_area_in = colSums(sum_area_bins, na.rm=T)
l_area_in = colSums(sum_time_area_bins, na.rm=T)
a_area_in = colSums(sum_age_area_bins, na.rm=T)
#m_area_in = colSums(sum_max_age_area_bins, na.rm=T)
m_area_in = sum_max_age_area_bins

#n_in = colSums(sum_state_bins, na.rm=T) #apply( sum_state_bins, 2, sum )
#l_in = colSums(sum_time_state_bins, na.rm=T) # apply( sum_time_state_bins, 2, sum )
#a_in = colSums(sum_age_state_bins, na.rm=T)
# posterior number of events leaving area-biome
#n_out = rowSums(sum_state_bins, na.rm=T) #apply( sum_state_bins, 1, sum )
#l_out = rowSums(sum_time_state_bins-diag(sum_time_state_bins), na.rm=T)#apply( sum_time_state_bins, 1, sum)
#a_out = rowSums(sum_age_state_bins, na.rm=T)

# filter out low freq histories
#n_in[n_in==0] = NA
#n_out[n_out==0] = NA
#l_in[l_in==0] = NA
#l_out[l_out==0] = NA

n_state_in = n_state_in * sum_in_freq
n_state_out = n_state_out * sum_in_freq
l_state_in = l_state_in * sum_in_freq
a_state_in = a_state_in * sum_in_freq
#a_out = a_out * sum_in_freq

# posterior number of events involving each area-biome
n_state_total = n_state_in + n_state_out
l_state_total = l_state_in + l_state_out
a_state_total = a_state_in + a_state_out

### first shift time vs. total time/n_events?


s_bg_lbl = rep(bg_lbl, 4)
s_biome_lbl = as.vector( sapply( biome_lbl, function(x) { rep(x, 6) }) )

dat = data.frame( n_state_total=n_state_total,
                  n_state_in=n_state_in,
                  n_state_flux=n_state_in-n_state_out,
                  n_state_ratio=n_state_in/n_state_out,
                  n_flux_state_ratio=(n_state_in-n_state_out)/n_state_total,
                  l_state_total=l_state_total,
                  l_state_in=l_state_in,
                  l_state_ratio=l_state_in/l_state_total,
                  a_state_total=a_state_total,
                  a_state_in=a_state_in,
                  a_state_out=a_state_out,
                  a_ratio=a_state_in/a_state_total,
                  da_in=a_state_in/n_state_in,
                  dl_in=l_state_in/n_state_in,
                  #da_biome_in=a_biome_in/n_biome_in,
                  #dl_biome_in=l_biome_in/n_biome_in,
                  #da_area_in=a_area_in/n_area_in,
                  #dl_area_in=l_area_in/n_area_in,
                  state=state_lbl, area=s_bg_lbl, biome=s_biome_lbl ) 
dat$biome = factor(x=dat$biome, levels=c("Tr","Wm","Cl","Fr"), ordered=T)
dat$area = factor(x=dat$area, levels=c("SEAs","EAs","Eur","NAm","CAm","SAm"), ordered=T)
#dat = dat[dat$biome=="Cl",]

s_biome_colors = as.vector(biome_colors$color)[1:4]
names(s_biome_colors) = levels(dat$biome)
s_bg_colors = as.vector(bg_colors$color)[1:6]
names(s_bg_colors) = levels(dat$area)

pp = ggplot(dat, aes(y=dl_in, x=da_in))
#pp = pp + geom_hline(yintercept=0, linetype=2, colour="gray")
pp = pp + geom_point( aes(colour=biome), shape=21, size=3.5)
pp = pp + geom_point( aes(colour=biome), shape=21, size=4)
pp = pp + geom_point( aes(colour=biome), shape=21, size=4.5)
pp = pp + geom_point( aes(colour=biome), shape=21, size=5)
pp = pp + geom_point( aes(colour=biome), shape=21, size=5.5)
pp = pp + geom_point( aes(colour=biome), shape=21, size=6)
pp = pp + geom_point( aes(fill=area), shape=21, size=3)
pp = pp + scale_colour_manual( values=s_biome_colors, name="Biome" )
pp = pp + scale_fill_manual( values=s_bg_colors, name="Area" )
#pp = pp + ylim(-1,1)
#pp = pp + xlab("Area-biome state duration\n(l_k / n_j)")
xlab = bquote('Mean area-biome shift length ('*l[k]/n[k]^'in'*')')
ylab = bquote('Mean area-biome shift age ('*a[k]/n[k]^'in'*')')
#ylab = bquote('Area-biome state flux ('*frac(n[k]^'in' - n[k]^'out', n[k]^'total')*')')#bquote('Area-biome state flux (('*n[k]^'in' - n[k]^'out'*')/'*n[k]*')')
pp = pp + ylab(xlab) + xlab(ylab)
#xlab = bquote('Assimilation ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')
#xlab = bquote('Area-biome state duration ('l[k] / n[k]')')
#pp = pp + ylab("Area-biome flux\n(n_in,j - n_out,j)")
#pp = pp + xlab( "Mean posterior # total events" )
#pp = pp + ylab( "Mean posterior relative flux \n(# incoming - # outgoing) / # total" )
#pp = pp + geom_point( aes(colour=area), shape=21, size=4)
pp



pp = ggplot(dat, aes(y=dl_in, x=da_in))
#pp = pp + geom_hline(yintercept=0, linetype=2, colour="gray")
pp = pp + geom_point( aes(fill=area), shape=21, size=3)
pp = pp + scale_colour_manual( values=s_biome_colors, name="Biome" )
pp = pp + scale_fill_manual( values=s_bg_colors, name="Area" )
#pp = pp + ylim(-1,1)
#pp = pp + xlab("Area-biome state duration\n(l_k / n_j)")
xlab = bquote('Mean area-biome shift length ('*l[k]/n[k]^'in'*')')
ylab = bquote('Mean area-biome shift age ('*a[k]/n[k]^'in'*')')
#ylab = bquote('Area-biome state flux ('*frac(n[k]^'in' - n[k]^'out', n[k]^'total')*')')#bquote('Area-biome state flux (('*n[k]^'in' - n[k]^'out'*')/'*n[k]*')')
pp = pp + ylab(xlab) + xlab(ylab)
#xlab = bquote('Assimilation ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')
#xlab = bquote('Area-biome state duration ('l[k] / n[k]')')
#pp = pp + ylab("Area-biome flux\n(n_in,j - n_out,j)")
#pp = pp + xlab( "Mean posterior # total events" )
#pp = pp + ylab( "Mean posterior relative flux \n(# incoming - # outgoing) / # total" )
#pp = pp + geom_point( aes(colour=area), shape=21, size=4)
pp



plot( a_biome_in/n_biome_in, l_biome_in/n_biome_in, col=s_biome_colors, pch=16)
legend( "topright", legend=names(s_biome_colors), col=s_biome_colors, pch=16 )

plot( a_area_in/n_area_in, l_area_in/n_area_in, col=s_bg_colors, pch=16)
legend( "topright", legend=names(s_bg_colors), col=s_bg_colors, pch=16 )

pdf_fn = paste(plot_fp, base_fn, ".flux.pdf", sep="")
CairoPDF( pdf_fn, height=6, width=6 )
print(pp)
dev.off()

png_fn = paste(plot_fp, base_fn, ".flux.png", sep="")
CairoPNG( png_fn, width = 6.5, height=6, units="in", dpi=300)
print(pp)
dev.off()

if (F) {

pca_dat = data.frame( n_in=n_in, l_in=l_in, a_in=a_in, n_out=n_out )
pca_dat = pca_dat[ pca_dat[,1] > 0, ]
#pca_dat = log(pca_dat)

pc = prcomp(pca_dat, scale.=T) #, center=T, scale=T)
comp = data.frame(pc$x[,1:3])

# k-means
n_clust = 4
k = kmeans(pca_dat, centers=n_clust, nstart=25, iter.max=1000)

#k_clust_orig = c("Punctuated","Radiation","Gradual")[as.vector(k$cluster)]
k_cluster = as.vector(k$cluster)
k_order = unique(k_cluster)
k_cluster[k_cluster==k_order[1]] = "A"
k_cluster[k_cluster==k_order[2]] = "B"
k_cluster[k_cluster==k_order[3]] = "C"
k_cluster[k_cluster==k_order[4]] = "D"
#print(k_cluster)

pcs = matrix( c(2,1,3,1,3,2), nrow=3, byrow=T )
g = list()
#sz = 8*c(-1,1)
clust_colors = c("orchid3","goldenrod","#7777DD","#DD7777")[c(1,4,3,2)]
for (i in 1:3) {
   #   g[[i]] <- ggbiplot(pc, groups=k_clust, obs.scale = 1, var.scale = 1, ellipse = TRUE, choices=pcs[i,], labels=taxon_labels, xlim=c(-5,5), ylim=c(-5,5), alpha=0.6)
   g[[i]] <- ggbiplot(pc, groups=k_cluster, 
                      var.scale=1,
                      circle=FALSE,
                      ellipse = TRUE, 
                      choices=pcs[i,], 
                      labels=rownames(pca_dat),
                      #xlim=sz, ylim=sz,
                      #alpha=alpha,
                      #show_legend=F,
                      loadings.label.size=40, label.size=40,size=40,varname.size=4) +
             scale_color_manual(name="Mode", values=clust_colors)
}

p = plot_grid(g[[1]], g[[2]], g[[3]], labels=c("A","B","C"), align="v", ncol=3) +
    theme(plot.margin = unit( c(0,0,0,0), "cm"))

}



#### GRAVEYARD

