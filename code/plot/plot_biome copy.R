stable = !T
if (stable) {
    library(RevGadgets)
} else {
    source("/Users/mlandis/projects/RevGadgets/R/plot_ancestral_states.R")
}
library(ggtree)



add_epoch_times <- function( p, max_age, x_offset=0, dy=4 ) {
    
    dy2 = dy+2
    max_x = max(p$data$x)
    max_y = max(p$data$y)
    epoch_names = c("Late\nCretaceous","Paleogene","Early\nEocene","Mid/Late\nEocene","Oligocene","Early\nMiocene","Mid/Late\nMiocene","Recent")
  

    x_pos = max_x-c(max_age, 65, 56, 48, 33.9, 23, 16, 5.3, 0)
    y_pos = rep(max_y, length(x_pos))
    y_pos[ length(x_pos) ] = 0
    print(y_pos)
    x_pos_mid = ( x_pos[1:(length(x_pos)-1)] + x_pos[2:length(x_pos)] ) / 2 

    t = "dotted"
    for (k in 2:(length(x_pos))) {
        #p = p + geom_vline(xintercept=x_pos[k], color="pink2", linetype=t, alpha=0.5)
        print(c(x_pos[k], y_pos[k]))
        p = p + geom_segment( x=x_pos[k], xend=x_pos[k], y=0-dy2, yend=y_pos[k], color="pink2", linetype=t, alpha=0.5)
    }
    for (k in 1:length(epoch_names)) {
        p = p + geom_text( label=epoch_names[k], x=x_pos_mid[k], y=0-dy, hjust=0.5, size=2)
        ##text(x=x_names[k]+1, y=-2, labels=epoch_names[k], cex=0.35, adj=0)
    }
    return(p)

}



args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    base_fn = "out.1.t163.f5.mask_fossil_states"
} else if (length(args)==1) {
    base_fn = args[1]
}

fp = "/Users/mlandis/projects/vib_div/"
out_fp = paste(fp, "output/", sep="")
plot_fp = paste(fp, "code/plot/fig/ase/", sep="")
out_str = paste( base_fn, ".biome", sep="" )
tree_fn = paste(out_fp, out_str, ".ase.tre", sep="")

n_biomes = 4

if (n_biomes == 4) {
    st_lbl = list("0"="Tropical",
                  "1"="W.Temp",
                  "2"="Cloud",
                  "3"="Fr.Temp",
                  "..."="...")
    #biome_colors = read.csv("biome_colors.n4.txt",header=T)
    st_colors = c("red","forestgreen","deepskyblue","darkblue","gray") #,"magenta","orange","darkgray","gray","gray")
} else if (n_biomes == 2) {
    st_lbl = list("0"="H",
                  "1"="F",
                  "..."="...")
    #biome_colors = read.csv("biome_colors.n4.txt",header=T)
    st_colors = c("red","darkblue","gray") #,"magenta","orange","darkgray","gray","gray")

}
names(st_colors) = st_lbl

    
summary_statistic = "PieState"
zz=plot_ancestral_states(tree_file=tree_fn,
                      summary_statistic=summary_statistic,
                      state_labels=st_lbl,
                      state_colors=st_colors,
                      node_label_size=0,
                      node_size_range=c(0.5,2.0),
                      node_label_nudge_x=0.5,
                      tip_node_size=0.75,
                      tip_label_size=1.0,
                      tip_label_offset=0.5,
                      xlim_visible=c(0,70),
                      show_posterior_legend=T,
                      node_pie_diameter=0.055,
                      tip_pie_diameter=0.042,
                      pie_nudge_x=0.3,
                      pie_nudge_y=0.45,
                      alpha=1)


p2 = zz
x_height = max(p2$data$x)
x_new_root = 80
x_offset = x_new_root - x_height
dy = 4
p2 = p2 + theme_tree2()
p2 = p2 + coord_cartesian(xlim = c(-(x_offset+5),x_height+3), expand=TRUE)
p2 = p2 + scale_x_continuous(breaks = seq(0,x_new_root,10)-x_offset, labels = rev(seq(0,x_new_root,10)))
p2 = p2 + labs(x="Age (Ma)")
p2 = p2 + theme(legend.position=c(0.07, 0.85), axis.line = element_line(colour = "black"))
p2 = add_epoch_times(p2, x_new_root, x_offset, dy)

ggsave(paste(plot_fp,out_str,".pdf",sep=""), height=8, width=8)
