####################
## Load libraries ##
####################

library(Cairo)
library(ggplot2)
library(igraph)

source("vib_div_util.R")

# default arguments
base_fn = "out.1.t163.f5"

# file system
fp           = "../../"
out_fp       = paste(fp, "output/", sep="")
plot_fp      = paste(fp, "code/plot/fig/", sep="")
bg_fn        = paste(out_fp, base_fn, ".bg.history.tsv", sep="")
biome_fn     = paste(out_fp, base_fn, ".biome.history.tsv", sep="")
col_bg_fn    = paste(fp, "code/plot/range_colors.n6.txt",sep="")
col_biome_fn = paste(fp, "code/plot/biome_colors.n4.txt",sep="")

# read files
bg_colors    = read.csv(col_bg_fn, header=T)
biome_colors = read.csv(col_biome_fn, header=T)
stoch_bg     = read.csv(bg_fn, sep="\t", stringsAsFactors=F)
stoch_biome  = read.csv(biome_fn, sep="\t", stringsAsFactors=F)

# define state space
n_areas  = 6
n_biomes = 4
n_states = n_areas * n_biomes

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
f_burn     = 0
thinby     = 1
iterations = unique(stoch_bg$iteration)
n_burn     = max(1, f_burn*length(iterations))
iterations = iterations[ n_burn:length(iterations) ]
iterations = iterations[ seq(1, length(iterations), by=thinby) ]
n_iter     = length(iterations)

#### Stage 1: merge list of biome and biogeo events into one structure ##
stoch_bg_biome = make_stoch_bg_biome(sample_bg, sample_biome, iterations)

#### Stage 2: reformat event list as LSTT matrix ##
mean_event = make_bg_biome_tx(stoch_bg_biome, 1, iterations, state_lbl)

#### Define and collect classes of event counts
class_vals = c(1, 2, 5, 15)
class_lbls = paste0(">", class_vals)
mean_event = round(mean_event, digits=1)
class_event = mean_event
class_event[ mean_event < class_vals[1] ] = 0
for (i in 1:length(class_vals)) {
    class_event[ mean_event >= class_vals[i] ] = i
}

### Reformat data
m_mean = melt(mean_event)
dat_plot = m_mean
colnames(dat_plot) = c("From_State","To_State", "Mean") #, "Prob")
dat_plot$From_Biome  = sapply( as.vector(dat_plot$From_State), function(x) { strsplit(x, "\\+")[[1]][1] })
dat_plot$From_Area   = sapply( as.vector(dat_plot$From_State), function(x) { strsplit(x, "\\+")[[1]][2] })
dat_plot$To_Biome    = sapply( as.vector(dat_plot$To_State),   function(x) { strsplit(x, "\\+")[[1]][1] })
dat_plot$To_Area     = sapply( as.vector(dat_plot$To_State),   function(x) { strsplit(x, "\\+")[[1]][2] })
cc                   = make_count_class(x=dat_plot$Mean, p=class_vals)
dat_plot$Count_Class = factor( cc, ordered=T, levels=class_lbls )
dat_plot$Present     = dat_plot$Mean > 0.0


### Pre-process graph
g <- graph.adjacency(class_event, weighted=T, mode="directed")

# assign edget weights, order, etc.
weights = E(g)$weight
edge_order = order(E(g)$weight)
new_edges = data.frame(from=NA, to=NA, weight=NA)
new_weights = rep(0, length(weights))
for (i in 1:length(edge_order)) {
    j = edge_order[i]
    new_edges[i,] = c( get.edgelist(g)[j,], E(g)$weight[j] )
}
new_edges$weight = as.numeric(new_edges$weight)

### Create the main graph object
state_lbl_alt = c( state_lbl[1:12], state_lbl[19:24], state_lbl[13:18] )
g <- graph_from_data_frame(d = new_edges, vertices = state_lbl_alt )

### Tweaking the coordinate system
coords = layout_in_circle(g)
coords = coords_rotate(coords, angle=3.5*(1/24*2*pi))

# create a little space among biomes
# Tropical
coords[01:06,1] = coords[01:06,1] + 0
coords[01:06,2] = coords[01:06,2] + 1/6
# Warm Temp
coords[07:12,1] = coords[07:12,1] - 1/6
coords[07:12,2] = coords[07:12,2] + 0
# Cold Temp
coords[13:18,1] = coords[13:18,1] - 0
coords[13:18,2] = coords[13:18,2] - 1/6
# Cloud
coords[19:24,1] = coords[19:24,1] + 1/6
coords[19:24,2] = coords[19:24,2] - 0

scale_factor = 9
coord_labels = (scale_factor*1.5)*coords
coords = scale_factor*coords

lab.locs = -apply( coords, 1, get_radian )
xlim<-c(-1.3*scale_factor-3,1.3*scale_factor+2)
ylim<-c(-1.3*scale_factor,1.3*scale_factor)

# colors
biome_colors      =  c( rep("red",6), rep("darkgreen",6), rep("darkblue", 6), rep("dodgerblue",6))
mark_biome_colors = c("red","darkgreen","dodgerblue","darkblue")
mark_biome_colors = adjustcolor( mark_biome_colors, alpha.f=0.2 )
area_colors       = rep( c("magenta","red","green","gold","cyan","blue"), 4)
mark_groups       = list( 1:6, 7:12, 19:24, 13:18 )

# apply/rescale colors with graph
V(g)$color     = biome_colors
V(g)$label.cex = 0.6
E(g)$color     = "black"

# determine edge color (darkness) scale
f_col       = 2
n_col       = length(class_lbls)
edge_colors = rev(gray.colors(n=n_col, start=0.0, end=0.7)) #[ E(g)$weight ]
edge_sizes  = c(8,16,25,35)/10
edge_curve  = 0.3
edge_curves = sample( seq(-0.3, -0.1, length.out = ecount(g)) )
edge_lty    = c(1,1,1,1) #c(3,2,5,1,1)

# plot figure
plot_pdf_fn = paste(plot_fp, base_fn, ".transition_graph.pdf", sep="")
CairoPDF( plot_pdf_fn, height=6, width=7 )
plot.igraph(g,
            asp=0,xlim=xlim,ylim=ylim,
            layout=coords,
            rescale=F,
            vertex.shape="fcircle", 
            vertex.color=area_colors,
            vertex.frame.color=biome_colors,
            vertex.frame.width=4,
            vertex.size=100,
            vertex.label.dist=50,
            vertex.label.degree=lab.locs,
            vertex.label.color="black",
            mark.groups=mark_groups,
            mark.col=mark_biome_colors,
            mark.border=mark_biome_colors,
            mark.expand=200,
            mark.shape=1/4,
            edge.width=edge_sizes[ E(g)$weight ], #E(g)$weight,
            edge.arrow.size=0.5, #0.5, #(E(g)$weight/max(E(g)$weight))^2,
            edge.arrow.width=1.5, #arrow_heads[ E(g)$weight ], #0.5, #(E(g)$weight/max(E(g)$weight))^2,
            edge.color=edge_colors[  E(g)$weight ],
            edge.curved=edge_curves,
            edge.lty=edge_lty[ E(g)$weight ])

legend(x=xlim[1]-3,y=ylim[1]+1,
       title="Num. events\n(posterior mean)",
       legend=class_lbls,
       box.lwd=0,
       lwd=edge_sizes,
       lty=edge_lty,
       cex=0.65,
       col=edge_colors)


dev.off()
