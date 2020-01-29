stable = !T
if (stable) {
    library(RevGadgets)
} else {
    source("/Users/mlandis/projects/RevGadgets/R/plot_ancestral_states.R")
}
library(ggtree)
library(ggimage)


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


# Arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    base_fn = "out.1.t163.f5"
} else if (length(args)==1) {
    base_fn = args[1]
}

# Files
fp = "/Users/mlandis/projects/vib_div/"
out_fp = paste(fp, "output/", sep="")
plot_fp = paste(fp, "code/plot/fig/ase/", sep="")
out_str = paste( base_fn, ".bg", sep="" )
col_fn = paste(plot_fp, "range_colors.n6.txt",sep="")
tree_fn = paste(out_fp, out_str, ".ase.tre", sep="")

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
source("/Users/mlandis/projects/RevGadgets/R/plot_ancestral_states.R")
n_expr = options()$expressions
options(expressions=50000)

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
                      tip_label_size=1.0,
                      tip_label_offset=0.25,
                      xlim_visible=c(0,70),
                      shoulder_label_nudge_x=-0.1,
                      show_posterior_legend=T,
                      node_pie_diameter=0.055,
                      tip_pie_diameter=0.04,
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
p2 = p2 + theme(legend.position=c(0.07, 0.7), axis.line = element_line(colour = "black"))
p2 = add_epoch_times(p2, x_new_root, x_offset, dy)
print(p2)
#options(expressions=n_expr)

#print(vv)
ggsave(paste(plot_fp,out_str,".pdf",sep=""), height=8, width=8)
options(expressions=n_expr)
