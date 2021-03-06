library(RevGadgets)
library(cowplot)
library(ggtree)
library(ggimage)

source("vib_div_util.R")

# Files
fp = "../../" #/Users/mlandis/projects/vib_div/"
base_fn = "out.2.t163.f5.mask_fossil_states"; iterations = 60000 + (0:5)*50
#base_fn = "out.1.t163.f5"; iterations = 44000 + (0:5)*50
out_fp = paste(fp, "output/", sep="")
plot_fp = paste(fp, "code/plot/fig/", sep="")
out_str = paste( base_fn, ".bg", sep="" )
tree_fn = paste(out_fp, out_str, ".ase.tre", sep="")
plot_fn = paste(plot_fp, "fig4_", out_str,".ase.pdf",sep="")

# add simmap to margin?
plot_simmap = TRUE
#iterations = 44000 + (0:5)*50

# State labels
st_lbl = list(
           "1"="SE As",
           "2"="E As",
           "3"="Eur",
           "4"="N Am",
           "5"="C Am",
           "6"="S Am",
           "A"="SE As + E As",
           "B"="SE As + Eur",
           "C"="E As + Eur",
           "D"="SE As + N Am",
           "E"="E As + N Am",
           "F"="Eur + N Am",
           "G"="SE As + C Am",
           "H"="E As + C Am",
           "I"="Eur + C Am",
           "J"="N Am + C Am",
           "K"="SE As + S Am",
           "L"="E As + S Am",
           "M"="Eur + S Am",
           "N"="N Am + S Am",
           "O"="C Am + S Am",
           "..."="...")

# State colors
bg_colors = read.csv(paste(fp, "code/plot/range_colors.n6.txt", sep=""), header=T)
st_colors = c( as.vector(bg_colors$color), "gray" )
names(st_colors) = st_lbl

# Build figure
cat("Processing...\n")
#source("/Users/mlandis/projects/RevGadgets/R/plot_ancestral_states.R")

# get tree
phy = read.nexus(tree_fn)
phy = fix_vib_tip(phy)

# plot
zz=plot_ancestral_states(tree_file=tree_fn,
                      include_start_states=T,
                      shoulder_label_size=0,
                      summary_statistic="PieRange",
                      state_labels=st_lbl,
                      state_colors=st_colors,
                      node_label_size=0,
                      node_size_range=c(0.5,2.0),
                      node_label_nudge_x=0.5,
                      tip_node_size=0.75,
                      tip_label_size=0.0,
                      tip_label_offset=0.25,
                      xlim_visible=c(0,70),
                      shoulder_label_nudge_x=-0.1,
                      show_posterior_legend=T,
                      node_pie_diameter=4.0,
                      tip_pie_diameter=3.0,
                      pie_nudge_x=0.2,
                      pie_nudge_y=0.2,
                      alpha=1)

p2 = zz
x_height = max(p2$data$x)
x_new_root = 75
x_offset = x_new_root - x_height
x_label = 30
dx = x_new_root %% 10
n_node = length(phy$tip.label)

# set coordinates
p2 = p2 + geom_rootedge(rootedge=5)
p2 = p2 + theme_tree2()
p2 = p2 + coord_cartesian(xlim = c(-(x_offset+8),x_height+x_label),
                          ylim=c(-7,n_node+1), expand=F)
p2 = p2 + scale_x_continuous(breaks = seq(0,x_new_root,10)-x_offset+dx,
                             labels = rev(seq(0,x_new_root,10)))
p2 = p2 + labs(x="Age (Ma)")

# set geological ages
p2 = add_epoch_times(p2, x_new_root,  dy_bars=-7, dy_text=-3)

# correct tip labels
tip_lbl = p2$data$label
tip_idx = which(!is.na(tip_lbl))
tip_lbl = tip_lbl[tip_idx]
tip_lbl = sapply( tip_lbl, function(x) { strsplit(x, split=" ")[[1]][2] })
tip_lbl[ tip_lbl=="stellato" ] = "stellato-tomentosum"
tip_lbl[ tip_lbl=="rigidum" ] = "rugosum"
p2$data$label[tip_idx] = tip_lbl

# assign tip colors
tip_colors = make_tip_colors()
p2 = p2 %<+% data.frame(cbind(tip_colors))
p2 = p2 + geom_tiplab(size=2.5, offset=0.5, color=tip_colors[ match(tip_lbl,names(tip_colors)) ])

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
    p2 = p2 + geom_cladelabel( node=clade_df$node[i],
                               label=clade_df$name[i],
                               offset= -10 + 11*clade_df$level[i], 
                               offset.text=offset.text,
                               barsize=1, #*clade_df$barsize[i],
                               extend=0.2,
                               hjust=hjust,
                               align=T, angle=angle, fontsize = 3.2 )
}

# correct legend fill
p2 = p2 + theme(plot.title = element_text(size=18, face="bold"),
                legend.position=c(0.02, 0.98),
                legend.justification=c("left", "top"),
                legend.key = element_blank(),
                axis.line = element_line(colour = "black"))

p2 = p2 + guides(colour=guide_legend(title="Range", override.aes=list(size=5)))


if (plot_simmap) {
    source("plot_simmap_grid.bg.R")
    pr2 = recordPlot() 
    dev.off()
    ps2 = plot_grid(NULL, pr2, NULL, NULL, ncol=2, nrow=2, rel_heights=c(17,1), rel_widths=c(1.5,5))
    #ps3 = plot_grid(NULL, ps2, ncol=2, nrow=2, rel_widths=c(1.5,5))
    #print(ps2)
    #psd = ggdraw(ps)
    pg = plot_grid(p2, ps2, ncol=2, axis="tr", align="h", rel_widths = c(12,3), labels=c("A)","B)"))
    plot_combine_fn = paste0(plot_fp, "fig4_", out_str,".ase_stoch.pdf")
    pdf(file=plot_combine_fn, height=18, width=18*(8.5/11))
    print(pg)
    dev.off()
} else {
    pdf(file=plot_fn, height=18, width=12)
    print(p2)
    dev.off()
}