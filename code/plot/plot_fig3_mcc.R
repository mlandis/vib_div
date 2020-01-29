library(ape)
library(ggplot2)
library(deeptime)
library(ggtree)
library(dplyr)

source("vib_div_util.R")

# default args
fp = "../../"
fn = "out.2.t163.f5.mask_fossil_states.mcc.tre"

# filepaths
phy_fp = paste(fp, "output/", sep="")
plot_fp = paste(fp, "code/plot/fig/", sep="")

# read tree
phy_fn = paste(phy_fp, fn, sep="")
#csv_fn = paste(fp, "data/emp/alpha/alphavirus_table.csv", sep="")
plot_fn = paste(plot_fp, "fig3_",fn,".pdf", sep="")

# read the tree in various formats
phy     = read.beast(phy_fn)
phy_rad = read.beast( paste(fp, "data/radseq_stage1.tre", sep=""))@phylo
phy_seq = read.beast( paste(fp, "data/radseq_cpdna_stage2.tre", sep=""))@phylo

# clean up tips
phy@phylo = fix_vib_tip(phy@phylo)
phy_rad   = fix_vib_tip(phy_rad)
phy_seq   = fix_vib_tip(phy_seq)

# gather information about clades
n_rad = phy_rad$tip.label
n_seq = setdiff( phy_seq$tip.label, n_rad )
idx_seq = match( n_seq, phy@phylo$tip.label )
n_unseq = setdiff( phy@phylo$tip.label, union(n_seq, n_rad))
n_unseq = n_unseq[ -match(c("CO","BC","IS","NWT","PB"), n_unseq) ]
idx_unseq = match( n_unseq, phy@phylo$tip.label )

# color tip labels based on whether taxon is sequenced
tip_colors = rep("black", length(phy@phylo$tip.label))
tip_colors[ idx_unseq ] = "#AAAAAA" 

# format posterior data
phy@data$posterior[ phy@data$posterior == 1 ] = NA

# Encountered problems with using geom_range to plot age HPDs in ggtree. It
# appears that geom_range incorrectly rotates the HPD relative to the height
# of the node unnecessarily. My guess for this would be because older version
# of ggtree primarily supported length measurements, and not height measurements
# so the new capability to handle height might contain a "reflection" bug.
#
# For example, suppose a node has height 3 with HPD [2, 7]. You can think of
# this offset as h + [2-h, 7-h]. ggtree seems to "rotate" this
# causing the HPD to appear as [-1, 4]. Figtree displays this correctly.
# 
# See this excellent trick by Tauana: https://groups.google.com/forum/#!msg/bioc-ggtree/wuAlY9phL9Q/L7efezPgDAAJ
# Adapted this code to also plot fossil tip uncertainty in red
#

phy@data$age_0.95_HPD = lapply(phy@data$age_0.95_HPD, function(z) { if (is.na(z)) { return(c(NA,NA)) } else { return(z) } })
minmax = t(matrix(unlist(phy@data$age_0.95_HPD), nrow=2))
bar_df = data.frame(node_id=as.integer(phy@data$node), as.data.frame(minmax))
names(bar_df) = c("node_id", "min", "max")
fossil_df = bar_df %>% filter(node_id %in% match(c("IS","PB","NWT","CO","BC"),phy@phylo$tip.label))
node_df = bar_df %>% filter(node_id > Ntip(phy@phylo))

# get dimensions
n_nodes = Nnode(phy)
max_age = 80
dx = max_age %% 10

# plot
pp = ggtree(phy, right=F )

# root edge
pp = pp + geom_rootedge(rootedge=5)

# plot age densities
node_df = node_df %>% left_join( pp$data, by=c("node_id"="node")) %>% select( node_id, min, max, y)
fossil_df = fossil_df %>% left_join( pp$data, by=c("node_id"="node")) %>% select( node_id, min, max, y)
pp = pp + geom_segment(aes(x=-min, y=y, xend=-max, yend=y), data=node_df, color="blue", size=1.5, alpha=0.3)
pp = pp + geom_segment(aes(x=-min, y=y, xend=-max, yend=y), data=fossil_df, color="red", size=1.5, alpha=0.4)

# set coordinates
pp = pp + coord_cartesian(xlim = c(-max_age,30), ylim=c(-7, n_nodes+1.5), expand=F)
pp = pp + scale_x_continuous(breaks = seq(-max_age-dx,0,10), labels = rev(seq(0,max_age+dx,10)))
pp = pp + theme_tree2()
pp = pp + labs(x="Age (Ma)")
pp = pp + theme(legend.position=c(.05, .955), axis.line = element_line(colour = "black"))
pp = revts(pp)
pp = add_epoch_times(pp, max_age, dy_bars=-7, dy_text=-3)

# plot clade support
pp$data$posterior_class = NA
pp$data$posterior_class[ which(pp$data$posterior>=0.99) ] = ">0.99"
pp$data$posterior_class[ which(pp$data$posterior<0.99&pp$data$posterior>0.95) ] = ">0.95"
pp$data$posterior_class[ which(pp$data$posterior<=0.95&pp$data$posterior>0.75) ] = ">0.75"
pp$data$posterior_class[ which(pp$data$posterior<=0.75&pp$data$posterior>0.5) ] = ">0.50"
pp$data$posterior_class[ which(pp$data$posterior<=0.5) ] = "<0.50"
pp$data$posterior_class = factor(pp$data$posterior_class, levels=c(">0.99",">0.95",">0.75", ">0.50", "<0.50"))
pp$data$label = sapply(pp$data$label, function(x) { gsub("_"," ",x) })
pp$data$label = sapply(pp$data$label, function(x) { gsub("subsp. ","",x) })

pp = pp + geom_nodepoint( data=pp$data[ !is.na(pp$data$posterior), ], size=2, color="black")
pp = pp + geom_nodepoint( data=pp$data[ !is.na(pp$data$posterior), ], aes(color=posterior_class), size=1.5)
col_post = c("#000000","#666666","#999999","#BBBBBB","#EEEEEE")
names(col_post) = levels(pp$data$posterior_class)
pp = pp + scale_color_manual(values=col_post, name="Posterior")

# set tip & clade names
pp = pp + geom_tiplab(size=2.5, offset=0.2, color=tip_colors)
clade_df = make_vib_clade_mtx(phy@phylo)
for (i in 1:nrow(clade_df)) {
    hjust = 0
    angle = 0
    offset.text = 1.1
    if (clade_df$level[i] > 2) {
        angle = 90
        offset.text = 1.5
        hjust = 0.5
    }
    pp = pp + geom_cladelabel( node=clade_df$node[i],
                               label=clade_df$name[i],
                               offset= -10 + 11*clade_df$level[i], 
                               offset.text=offset.text,
                               barsize=1,
                               extend=0.2,
                               hjust=hjust,
                               align=T, angle=angle, fontsize = 3.2 )
}

# print pdf
pdf(height=18, width=10, file=plot_fn)
print(pp)
dev.off()

#print(pp)


