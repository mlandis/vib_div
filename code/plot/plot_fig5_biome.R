library(RevGadgets)
library(cowplot)
library(ggtree)
library(ggimage)

source("vib_div_util.R")

# Files
fp      = "../../"
base_fn = "out.2.t163.f5.mask_fossil_states"; iterations = 60000 + (0:5)*50
#base_fn = "out.1.t163.f5"; iterations = 44000 + (0:5)*50
out_fp  = paste(fp, "output/", sep="")
plot_fp = paste(fp, "code/plot/fig/", sep="")
out_str = paste( base_fn, ".biome", sep="" )
tree_fn = paste(out_fp, out_str, ".ase.tre", sep="")

# add simmap to margin?
plot_simmap = TRUE
#iterations = 44000 + (0:5)*50
#iterations = c(99000, 99050, 99100, 99150, 99200, 99250)

# set up state labels
n_biomes = 4
if (n_biomes == 4) {
    st_lbl = list("0"="Tropical",
                  "1"="Warm Temp.",
                  "2"="Cloud",
                  "3"="Cold Temp.",
                  "..."="...")
    st_colors = c("red","forestgreen","deepskyblue","darkblue","gray") #,"magenta","orange","darkgray","gray","gray")
} else if (n_biomes == 2) {
    st_lbl = list("0"="H",
                  "1"="F",
                  "..."="...")
    st_colors = c("red","darkblue","gray") #,"magenta","orange","darkgray","gray","gray")
}
names(st_colors) = st_lbl

# Build figure
cat("Processing...\n")

phy = read.nexus(tree_fn)
phy = fix_vib_tip(phy)
    
summary_statistic = "PieState"
zz=plot_ancestral_states(tree_file=tree_fn,
                      summary_statistic=summary_statistic,
                      state_labels=st_lbl,
                      state_colors=st_colors,
                      node_label_size=2,
                      node_size_range=c(0.5,2.0),
                      node_label_nudge_x=1.0,
                      node_label_nudge_y=3.5, #2.5
                      tip_node_size=2,
                      tip_label_size=0,
                      tip_label_offset=0.5,
                      xlim_visible=c(0,70),
                      show_posterior_legend=T,
                      node_pie_diameter=4.0,
                      tip_pie_diameter=3.0,
                      pie_nudge_x=0.2,
                      pie_nudge_y=0.2,
                      alpha=1)


p1 = zz
x_height = max(p1$data$x)
x_new_root = 75
x_offset = x_new_root - x_height
x_label = 30
dx = x_new_root %% 10
n_node = length(phy$tip.label)

# set coordinates
p1 = p1 + geom_rootedge(rootedge=5)
p1 = p1 + theme_tree2()
p1 = p1 + coord_cartesian(xlim = c(-(x_offset+8),x_height+x_label), ylim=c(-7,n_node+1), expand=F)
p1 = p1 + scale_x_continuous(breaks = seq(0,x_new_root,10)-x_offset+dx, labels = rev(seq(0,x_new_root,10)))
p1 = p1 + labs(x="Age (Ma)")

# set geological ages
p1 = add_epoch_times(p1, x_new_root,  dy_bars=-7, dy_text=-3)

# correct tip labels
tip_lbl = p1$data$label
tip_idx = which(!is.na(tip_lbl))
tip_lbl = tip_lbl[tip_idx]
tip_lbl = sapply( tip_lbl, function(x) { strsplit(x, split=" ")[[1]][2] })
tip_lbl[ tip_lbl=="stellato" ] = "stellato-tomentosum"
tip_lbl[ tip_lbl=="rigidum" ] = "rugosum"
p1$data$label[tip_idx] = tip_lbl

# assign tip colors
tip_colors = make_tip_colors()
p1 = p1 + geom_tiplab(size=2.5, offset=0.5, color=tip_colors[ match(tip_lbl,names(tip_colors)) ])

# Add clade labels
clade_df = make_vib_clade_mtx(phy)
for (i in 1:nrow(clade_df)) {
    hjust = 0
    angle = 0
    offset.text = 1.1
    if (clade_df$level[i] > 2) {
        angle = 90
        offset.text = 1.5
        hjust = 0.5
    }
    p1 = p1 + geom_cladelabel( node=clade_df$node[i],
                               label=clade_df$name[i],
                               offset= -10 + 11*clade_df$level[i], 
                               offset.text=offset.text,
                               barsize=1, #*clade_df$barsize[i],
                               extend=0.2,
                               hjust=hjust,
                               align=T, angle=angle, fontsize = 3.2 )
}


# correct legend fill
p1 = p1 + theme(plot.title = element_text(size=18, face="bold"),
                legend.position=c(0.02, 0.98),
                legend.justification=c("left", "top"),
                legend.key = element_blank(),
                axis.line = element_line(colour = "black"))

p1 = p1 + guides(colour=guide_legend(title="Biome", override.aes=list(size=5)))

# fix run-on x-axis??
if (plot_simmap) {
    source("plot_simmap_grid.biome.R")
    pr1 = recordPlot() 
    dev.off()
    ps1 = plot_grid(NULL, pr1, NULL, NULL, ncol=2, nrow=2, rel_heights=c(17,1), rel_widths=c(1.5,5))
    #ps3 = plot_grid(NULL, ps2, ncol=2, nrow=2, rel_widths=c(1.5,5))
    #print(ps2)
    #psd = ggdraw(ps)
    pg = plot_grid(p1, ps1, ncol=2, axis="tr", align="h", rel_widths = c(12,3), labels=c("A)","B)"))
    plot_combine_fn = paste0(plot_fp, "fig5_", out_str,".ase_stoch.pdf")
    pdf(file=plot_combine_fn, height=18, width=18*(8.5/11))
    print(pg)
    dev.off()
} else {
    pdf(file=plot_fn, height=18, width=12)
    print(p2)
    dev.off()
}