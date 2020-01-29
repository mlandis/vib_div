####################
## Load libraries ##
####################

library(cowplot)
library(ggplot2)
library(grid)
library(scales)
library(Cairo)
library(data.table)
library(igraph)
library(ggnetwork)
library(ggnewscale)


######################
## Define functions ##
######################

# functions
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

radian.rescale <- function(x, start=0, direction=1, PI=3.141592653589793) {
  c.rotate <- function(x) (x + start) %% (2 * PI) * direction
  c.rotate(scales::rescale(x, from=c(0, 2 * PI), to=range(x)))
}

get_compound_state = function(area, biome, n_areas=6, n_biomes=4) {
    return( (n_areas)*(biome-1) + area )
}

# Returns a matrix of transition events indexed by compound state identities.
# Each row of the matrix
make_compound_tx = function( from_bg_idx, from_biome_idx, to_bg_idx, to_biome_idx ) {
    ret = matrix(NA, nrow=0, ncol=2)
    
    is_biome_event = !identical(from_biome_idx, to_biome_idx)
    is_bg_event = !identical(from_bg_idx, to_bg_idx)
    if (is_biome_event && is_bg_event) {
        cat("from_bg_idx    =", from_bg_idx, "\n")
        cat("to_bg_idx      =", to_bg_idx, "\n")
        cat("from_biome_idx =", from_biome_idx, "\n")
        cat("to_biome_idx   =", to_biome_idx, "\n")
        error("Event identified as both a biome and bg event!")
    }
    is_dispersal = FALSE
    if (is_bg_event && (length(from_bg_idx) < length(to_bg_idx))) {
        is_dispersal = TRUE
    }
    
    for (from_bg in from_bg_idx) {
        for (to_bg in to_bg_idx) {
            for (from_biome in from_biome_idx) {
                for (to_biome in to_biome_idx) {
                    if (is_biome_event && from_bg == to_bg) {
                        # only identify biome shifts within each region
                        from_state = get_compound_state(from_bg, from_biome)
                        to_state = get_compound_state(to_bg, to_biome)
                        ret = rbind(ret, c(from_state, to_state))
                    } else if (is_dispersal && from_bg != to_bg) {
                        from_state = get_compound_state(from_bg, from_biome)
                        to_state = get_compound_state(to_bg, to_biome)
                        ret = rbind(ret, c(from_state, to_state))
                    }
                }
            }
        }
    }
    colnames(ret) = c("from","to")
    return(ret)
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


check_bin_diag_zero = function(x) {
    n_states = dim(x)[1]
    n_bins = dim(x)[3]
    for (i in 1:n_states) {
        for (j in 1:n_bins) {
            if (x[i,i,j] != 0) {
                error("bin",j,"+","state",i,"=",x[i,i,j],"> 0\n")
            }
        }
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
    
    #cat("......................\n")
    #cat("Step 1\n")
    #print(b);cat("\n")
    #print(b1_tmp);cat("\n")
    #print(b2_tmp);cat("\n")
    
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
    
    #cat("......................\n")
    #cat("Step 2\n")
    #print(b);cat("\n")
    
    
    if (nrow(b)==1) {
        x1 = b$branch_start_time
        x2 = b$branch_end_time
        #b = cbind(b, x1=b$branch_start_time, x2=b$branch_end_time, s1$start_state)
        #cat("......................\n")
        #cat("Step 3a\n")
        #print(b);cat("\n")
        
    } else {
        x1 = c( b$branch_start_time[1], b$transition_time )
        x2 = c( b$transition_time, b$branch_end_time[1] )
        b = rbind(as.data.frame(b, stringsAsFactors=F), 
                  as.data.frame(b[nrow(b),], stringsAsFactors=F))
        b[nrow(b),c("start_state_1","start_state_2","transition_time")] = b[nrow(b),c("end_state_1","end_state_2","branch_end_time")]
        #cat("......................\n")
        #cat("Step 3b\n")
        #print(b);cat("\n")
    }
    b = cbind(b, x1=x1, x2=x2, s1=b$start_state_1, s2=b$start_state_2)
    #cat("Step 4\n")
    #print(b);cat("\n\n")
    #print("after")
    #print(b)
    return(b)
    #return(list(b=b,b_orig=b_orig))
    
}


# default arguments
args = commandArgs(trailingOnly=TRUE)
if        (length(args)==0) {
    base_fn = "out.1.t163.f5"
} else if (length(args)==1) {
    base_fn = args[1]
}

# file system
fp           = "/Users/mlandis/projects/vib_div/"
out_fp       = paste(fp, "output/", sep="")
plot_fp      = paste(fp, "code/plot/fig/matrix/", sep="")
bg_fn        = paste(out_fp, base_fn, ".bg.history.tsv", sep="")
biome_fn     = paste(out_fp, base_fn, ".biome.history.tsv", sep="")
col_bg_fn    = paste(fp, "code/plot/range_colors.n6.txt",sep="")
col_biome_fn = paste(fp, "code/plot/biome_colors.n4.txt",sep="")

# read files
bg_colors    = read.csv(col_bg_fn, header=T)
biome_colors = read.csv(col_biome_fn, header=T)
stoch_bg     = read.csv(bg_fn, sep="\t", stringsAsFactors=F)
stoch_biome = read.csv(biome_fn, sep="\t", stringsAsFactors=F)

# define state space
n_areas  = 6
n_biomes = 4
n_states = n_areas * n_biomes
branches = 1:max(unique(stoch_bg$parent_index), na.rm=T)

# labels
bg_lbl = c("SEAs","EAs","Eur","NAm","CAm","SAm")
biome_lbl = c("Trop","Warm","Cloud","Cold")
state_lbl = c()
for (i in 1:n_biomes) {
    for (j in 1:n_areas) {
        state_lbl = c(state_lbl, paste(biome_lbl[i],bg_lbl[j],sep="+"))
    }
}
from_state_lbl = paste("from",state_lbl,sep="_")
to_state_lbl = paste("to",state_lbl,sep="_")



# filter files for events
stoch_bg                                                                = stoch_bg[ stoch_bg$transition_type != "cladogenetic", ]
stoch_bg$transition_time[ stoch_bg$transition_type=="no_change" ]       = stoch_bg$branch_start_time[ stoch_bg$transition_type=="no_change" ]
stoch_biome$transition_time[ stoch_biome$transition_type=="no_change" ] = stoch_biome$branch_start_time[ stoch_biome$transition_type=="no_change" ]

# filter files for number of samples
f_burn     = 0.0
iterations = unique(stoch_bg$iteration)
n_burn     = max(1, f_burn*length(iterations))
thinby     = 1
iterations = iterations[n_burn:length(iterations)]
iterations = iterations[ seq(1, length(iterations), by=thinby) ]
#iterations = iterations[ length(iterations) ]
n_iter     = length(iterations)

#######################################################################
## Stage 1: merge list of biome and biogeo events into one structure ##
#######################################################################

stoch_list = list()
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
        nd_idx       = branches[j]
        branch_bg    = sample_bg[ sample_bg$node_index==nd_idx, ]
        branch_biome = sample_biome[ sample_biome$node_index==nd_idx, ]
        
        # interleave biome and biogeography stochastic mappings
        tmp_branch_list[[j]] = as.data.frame( make_joint_state(branch_bg, branch_biome), stringsAsFactors=F )
    }
    stoch_list[[i]] = rbindlist(tmp_branch_list) 
}
# merge all events
stoch_bg_biome = rbindlist(stoch_list)
# filter for anagenetic events
stoch_bg_biome = stoch_bg_biome[ stoch_bg_biome$transition_type == "anagenetic", ]


#################################################
## Stage 2: reformat event list as LSTT matrix ##
#################################################

# define state bin deimensions
# index bins as [from-area-biome] x [to-area-biome] x time
bin_width  = 1
#max_time   = 90
#n_bins     = max_time / bin_width
state_bins = array(0, dim=c(n_states, n_states, n_iter))
hit_bins   = array(0, dim=c(n_states, n_states, n_iter))
#ages       = seq(0.0, max_time, by=bin_width)

# populate the LSTT matrix
idx_tmp = 1
curr_it = -1
for (i in 1:nrow(stoch_bg_biome)) {
    
    if (curr_it != stoch_bg_biome$iteration[i]) {
        curr_it = stoch_bg_biome$iteration[i]
        cat("Stage 2, processing iteration ",curr_it," / ", max(stoch_bg_biome$iteration), "\n", sep="")
    }
    it_idx = match(curr_it, iterations)
    
    if (stoch_bg_biome$transition_type[i] == "anagenetic") {

        from_bg_idx      = get_bg_state( stoch_bg_biome$start_state_1[i] )
        from_biome_idx   = get_biome_state( stoch_bg_biome$start_state_2[i] )
        to_bg_idx        = get_bg_state( stoch_bg_biome$end_state_1[i] )
        to_biome_idx     = get_biome_state( stoch_bg_biome$end_state_2[i] )
        #time_idx         = as.integer( stoch_bg_biome$transition_time[i] )
       
        # cat("sample\n")
        # print(stoch_bg_biome[i,])
        # cat("iteration\n")
        # print(curr_it)
        # cat("time_idx\n")
        # print(time_idx)
        # cat("from_biome_idx\n")
        # print(from_biome_idx)
        # cat("to_biome_idx\n")
        # print(to_biome_idx)
        # cat("from_bg_idx\n")
        # print(from_bg_idx)
        # cat("to_bg_idx\n")
        # print(to_bg_idx)
        # #cat("from_state_idx\n")
        # #print(from_state_idx)
        # #cat("to_state_idx\n")
        # #print(to_state_idx)
        # cat("\n")

        # get the matrix of state transitions
        state_tx = make_compound_tx( from_bg_idx, from_biome_idx, to_bg_idx, to_biome_idx )
        
        # store each event's info as count or hit
        if (nrow(state_tx) > 0) {
            for (i in 1:nrow(state_tx)) {
                state_bins[ state_tx[i,1], state_tx[i,2], it_idx ] = state_bins[ state_tx[i,1], state_tx[i,2], it_idx ] + 1
                hit_bins[ state_tx[i,1], state_tx[i,2], it_idx ] = 1
            }
        }
    }
}

# validate bins don't include self transitions
check_bin_diag_zero(state_bins)

# convert from posterior samples
mean_event = rowSums(state_bins,dims=2) * (1/n_iter)
prob_event = rowSums(hit_bins,dims=2) * (1/n_iter)

# format & filter
rownames(mean_event)=state_lbl; colnames(mean_event)=state_lbl
rownames(prob_event)=state_lbl; colnames(prob_event)=state_lbl




# threshold values, below which values set to zero
th_mean = 1
th_prob = 0.05
for (i in 1:n_states) {
    mean_event[ mean_event<th_mean ] = 0
    prob_event[ prob_event<th_prob ] = 0
}
#mean_event = mean_event / max(mean_event)
#mean_event = mean_event * 5
mean_event = round(mean_event, digits=3)
prob_event = round(prob_event, digits=3)





#################
# another style #
#################

make_count_class = function(x, p) {
    n = length(p)
    y = rep(0, length(x))
    for (i in 1:length(x)) {
        z = sum(x[i] > p)
        if (z == 0) {
            y[i] = NA #paste("<", p[1], sep="")
        } else {
            y[i] = paste(">", p[z], sep="")
        }
    }
   return(y)
}

m_mean = melt(mean_event)
m_prob = melt(prob_event)
dat_plot = m_mean # cbind( melt(mean_event), melt(prob_event)[3] )
colnames(dat_plot) = c("From_State","To_State", "Mean") #, "Prob")
#dat_plot$From_State  = factor(dat_plot$From_State, ordered=T, rev(levels(dat_plot$From_State)))
#dat_plot$To_State    = factor(dat_plot$To_State, ordered=T, rev(levels(dat_plot$To_State)))
dat_plot$Log_Mean    = log(dat_plot$Mean)
dat_plot$From_Biome  = sapply( as.vector(dat_plot$From_State), function(x) { strsplit(x, "\\+")[[1]][1] })
dat_plot$From_Area   = sapply( as.vector(dat_plot$From_State), function(x) { strsplit(x, "\\+")[[1]][2] })
dat_plot$To_Biome    = sapply( as.vector(dat_plot$To_State),   function(x) { strsplit(x, "\\+")[[1]][1] })
dat_plot$To_Area     = sapply( as.vector(dat_plot$To_State),   function(x) { strsplit(x, "\\+")[[1]][2] })
cc = make_count_class(x=dat_plot$Mean, p=c(2,5,10,15,30))
dat_plot$Count_Class = factor( cc, ordered=T, levels=c(">2",">5",">10",">15",">30") )
dat_plot$Present = dat_plot$Mean > 1


biome_colors =  c( rep("red",6), rep("darkgreen",6), rep("dodgerblue",6), rep("darkblue", 6))
names(biome_colors) = c( rep("Trop",6), rep("Warm",6), rep("Cloud",6), rep("Cold", 6))
area_colors = rep( c("magenta","red","green","gold","cyan","blue"), 4)
names(area_colors) = rep( c("SEAs","EAs","Eur","NAm","CAm","SAm"), 4)

#dat_plot = dat_plot[dat_plot$Present, ]
#dat_plot = dat_plot[ dat_plot$Mean > 2, ]
#n_from = unique(dat_plot$From_State)
#n_to = unique(dat_plot$To_State)
#dat_plot$From_Idx = 


rangeTransform = function(x) (x - min(x)) / (max(x) - min(x))

pp = ggplot()
#pp = pp + scale_y_continuous(limits = c(-1,25) )
#pp = pp + annotate(geom="rect", xmin= -0.2, xmax=0.2, ymin=-0.2, ymax=6.2, alpha=0.2)
pp = pp + geom_point( data=dat_plot, mapping=aes( x=0, y=From_State, col=From_Biome), size=9)
pp = pp + geom_point( data=dat_plot, mapping=aes( x=1, y=To_State,   col=To_Biome ), size=9)
pp = pp + scale_color_manual( values=biome_colors) #, labels=unique(names(biome_colors)) )
pp = pp + new_scale_color()

pp = pp + geom_point( data=dat_plot, mapping=aes( x=0, y=From_State), size=6, col="white")
pp = pp + geom_point( data=dat_plot, mapping=aes( x=1, y=To_State), size=6, col="white")
#pp = pp + scale_color_manual( values=biome_colors, labels=unique(names(biome_colors)) )
pp = pp + new_scale_color()

pp = pp + geom_point( data=dat_plot, mapping=aes( x=0, y=From_State, col=From_Area), size=5)
pp = pp + geom_point( data=dat_plot, mapping=aes( x=1, y=To_State,   col=To_Area ), size=5)
pp = pp + scale_color_manual( values=area_colors) #, labels=unique(names(area_colors)) )
pp = pp + new_scale_color()

pp = pp + geom_text( data=dat_plot,  mapping=aes( x=0, y=From_State, label=From_State), col="white", size=0.7 )
pp = pp + geom_text( data=dat_plot,  mapping=aes( x=1, y=To_State,   label=To_State), col="white", size=0.7 )
#pp = pp + scale_color_manual(

my_arrow = arrow(length=unit(0.15,"cm"), ends="last", type = "open")
arrow_cols = c("#DDDDDD","#BBBBBB","#999999","#777777","#000000")
names(arrow_cols) = c(">2",">5",">10",">15",">30")

n_cc_lvl = names(arrow_cols)
for (i in 1:length(n_cc_lvl)) {
    d_tmp = dat_plot[ dat_plot$Count_Class==n_cc_lvl[i],]
    if (nrow(d_tmp) > 0) {
        pp = pp + geom_segment( data=d_tmp,
                            mapping=aes( x=0.05, xend=0.95, y=From_State, yend=To_State, col=Count_Class, size=Count_Class),
                            arrow=my_arrow, size=1)
    }
}

pp = pp + scale_colour_manual(values=arrow_cols)
#pp = pp + scale_color_gradient( low="#DDDDDD", high="#000000")
pp = pp + theme_void()
pp = pp + theme( axis.text = element_blank(),
                 axis.ticks = element_blank() )
pp = pp + coord_cartesian(expand=T)
pp

plot_fn = paste( fp, "code/plot/fig/tx/fig8_tx_graph.pdf", sep="")
pdf( plot_fn, height=6, width=5 )
print(pp)
dev.off()

#pp = pp + geom_segment( data=dat_plot[ dat_plot$Prob > 0.50, ], mapping=aes( x=0, xend=1, y=From_State, yend=To_State, lwd=Mean), col="gray", arrow=arrow())
#pp = pp + geom_segment( data=dat_plot[ dat_plot$Prob > 0.95, ], mapping=aes( x=0, xend=1, y=From_State, yend=To_State, lwd=Mean), col="black", arrow=arrow())
#pp = pp + annotate( "text", x=0, y=24:1, label=state_lbl)
#pp = pp + annotate( "text", x=1, y=24:1, label=state_lbl)
  
#, aes(x=To_State, y=From_State, fill=Mean) )
#pp + geom_tile()










#mean_event = mean_event / sum(mean_event)


###################
## Stage 4: plot ##
###################

class_event = mean_event
class_event[ mean_event < 01 ] = 0
class_event[ mean_event > 01 ] = 1
class_event[ mean_event > 05 ] = 2
class_event[ mean_event > 10 ] = 3
class_event[ mean_event > 15 ] = 4
class_event[ mean_event > 30 ] = 5





# build graph
g <- graph.adjacency(class_event, weighted=T, mode="directed")
# assign edget weights, order, etc.
weights = E(g)$weight
nm_weight = c(">1",">5",">10",">15",">30")
edge_order = order(E(g)$weight)
new_edges = data.frame(from=NA, to=NA, weight=NA)
new_weights = rep(0, length(weights))
for (i in 1:length(edge_order)) {
    j = edge_order[i]
    new_edges[i,] = c( get.edgelist(g)[j,], E(g)$weight[j] )
}
new_edges$weight = as.numeric(new_edges$weight)

# reformat graph
g <- graph_from_data_frame(d = new_edges, vertices = state_lbl )

# coordinates
coords = layout_in_circle(g)
if (T) {
    # Tropical
    coords[01:06,1] = coords[01:06,1] + 1/6 # + c(0.5, 0.5)
    coords[01:06,2] = coords[01:06,2] + 1/6 # + c(0.5, 0.5)
    # Warm Temp
    coords[07:12,1] = coords[07:12,1] - 1/6 # + c(0.5, 0.5)
    coords[07:12,2] = coords[07:12,2] + 1/6 # + c(0.5, 0.5)
    # Cloud
    coords[13:18,1] = coords[13:18,1] - 1/6 # + c(0.5, -0.5)
    coords[13:18,2] = coords[13:18,2] - 1/6 # + c(0.5, -0.5)
    # Cold Temp
    coords[19:24,1] = coords[19:24,1] + 1/6 # + c(0.5, 0.5)
    coords[19:24,2] = coords[19:24,2] - 1/6 # + c(0.5, 0.5)
}
coords_rotate = function(x, angle=0, origin=c(0,0)) {
    for (i in 1:nrow(x)) {
        xy = x[i,]
        xy_new = xy
        xy_new[1] = origin[1] + cos(angle) * (xy[1] - origin[1]) - sin(angle) * (xy[2] - origin[2])
        xy_new[2] = origin[2] + sin(angle) * (xy[1] - origin[1]) + cos(angle) * (xy[2] - origin[2])
        x[i,] = xy_new
    }
    return(x)
}
coords = coords_rotate(coords, angle=1/24/2*2*pi)
coords[,1]-0.3
scale_factor = 9
coord_labels = (scale_factor*1.5)*coords
coords = scale_factor*coords
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}
lab.locs = radian.rescale(x=1:24, start=0, direction=-1)
xlim<-c(-1.3*scale_factor-3,1.3*scale_factor+2)
ylim<-c(-1.3*scale_factor,1.3*scale_factor)

# colors
biome_colors =  c( rep("red",6), rep("darkgreen",6), rep("dodgerblue",6), rep("darkblue", 6))
mark_biome_colors = c("red","darkgreen","dodgerblue","darkblue")
mark_biome_colors = adjustcolor( mark_biome_colors, alpha.f=0.2 )
area_colors = rep( c("magenta","red","green","gold","cyan","blue"), 4)
mark_groups = list( 1:6, 7:12, 13:18, 19:24 )

# apply/rescale colors with graph
V(g)$color = biome_colors
V(g)$label.cex = 0.6
E(g)$color = "black"

# determine edge color (darkness) scale
f_col = 2
n_col = length(nm_weight)
edge_colors = rev(gray.colors(n=n_col, start=0.0, end=0.9)) #[ E(g)$weight ]
edge_sizes = c(6,10,18,25,35)/10
edge_curve = 0.3
arrow_heads = c(2,5,10,15,30)/10
#edge_curves = seq(-1*edge_curve, 1*edge_curve, length = ecount(g))
edge_curves = sample( seq(-0.3, -0.1, length.out = ecount(g)) )
edge_lty = c(3,2,5,1,1)
    
# plot system
plot_pdf_fn = paste(plot_fp, base_fn, ".transition_graph.pdf", sep="")
CairoPDF( plot_pdf_fn, height=6, width=7 )
plot.igraph(g,
            asp=0,xlim=xlim,ylim=ylim,
            layout=coords,
            vertex.size=80,
            rescale=F,
            vertex.color=area_colors,
            vertex.frame.color=biome_colors,
            vertex.frame.width=30,
            vertex.label.dist=35,
            vertex.label.degree=lab.locs,
            mark.groups=mark_groups,
            mark.col=mark_biome_colors,
            mark.border=mark_biome_colors,
            mark.expand=100,
            mark.shape=-1/4,
            edge.width=edge_sizes[ E(g)$weight ], #E(g)$weight,
            edge.arrow.size=0.5, #0.5, #(E(g)$weight/max(E(g)$weight))^2,
            edge.arrow.width=1.75, #arrow_heads[ E(g)$weight ], #0.5, #(E(g)$weight/max(E(g)$weight))^2,
            edge.color="black", #edge_colors[  E(g)$weight ],
            edge.curved=edge_curves,
            edge.lty=edge_lty[ E(g)$weight ])

legend(x=xlim[1]-5,y=ylim[1]+2,
       title="Num. changes",
       legend=nm_weight,
       box.lwd=0,
       lwd=edge_sizes,
       lty=edge_lty,
       cex=0.65,
       col="black") #,pt.cex=scaled,col='black',pch=21, pt.bg='orange')


dev.off()

# 
# 
# plot_png_fn = paste(plot_fp, base_fn,".transition_graph.png", sep="")
# CairoPNG( plot_png_fn, height=6, width=6, res=300, units="in")
# #dev.off()
# 
# 
# #par(mar=c(0,0,0,0)+6, oma=c(0,0,0,0))
# plot.igraph(g,
#             asp=0,xlim=xlim,ylim=ylim,
#             layout=coords,
#             vertex.size=80,
#             rescale=F,
#             vertex.color=area_colors,
#             vertex.frame.color=biome_colors,
#             vertex.frame.width=30,
#             vertex.label.dist=30,
#             vertex.label.degree=lab.locs,
#             mark.groups=mark_groups,
#             mark.col=mark_biome_colors,
#             mark.border=NA,
#             mark.expand=100,
#             mark.shape=1/2,
#             edge.labels=E(g)$weight,
#             edge.width=2, #E(g)$weight,
#             edge.arrow.size=0.5, #(E(g)$weight/max(E(g)$weight))^2,
#             edge.color=edge_colors,
#             edge.curved=seq(-1*curve, 1*curve, length = ecount(g)))
# 
# #legend('left',legend=seq(0.2,1.0,by=0.1), lwd=2, col=rev(gray.colors(n=6, start=0.0, end=0.8))) #,pt.cex=scaled,col='black',pch=21, pt.bg='orange')
# 
# 
# dev.off()
# 



# 







# GRAVEYARD

# option_str = c("all","into_LatAm","into_Cloud")[1]
# if (option_str=="all"){
#     # do nothing
#     curve=0.2
# } else if (option_str=="into_LatAm") {
#     # into Latin America from anywhere
#     mean_event[ ,c(1:4, 7:10, 13:16, 19:22) ] = 0
#     mean_event[ c(5:6, 11:12, 17:18, 23:24), ] = 0
#     curve=0.2
# } else if (option_str=="into_Cloud") {
#     mean_event[ , c(1:12, 19:24) ] = 0
#     mean_event[ 13:18, ] = 0
#     curve=0.2
# }
# q_total_ij = n_states^2-n_states #(n_states) * (n_areas-1 + n_biomes-1)
# q_total = sum(sum_state_bins)
# q_mean = q_total / (n_states^2-n_states)
