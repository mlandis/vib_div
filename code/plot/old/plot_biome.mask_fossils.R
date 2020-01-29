stable = !T
if (stable) {
    library(RevGadgets)
} else {
    source("/Users/mlandis/projects/RevGadgets/R/plot_ancestral_states.R")
}
library(ggtree)

out_fp = "/Users/mlandis/projects/vib_div/"
out_str = "output/tinus/out.1.viburnum.t163.f5.biome"
out_str = "output/tinus/out.1.viburnum.t163.f5.mask_fossil_states.biome"
st_lbl = list("0"="H",
              "1"="L",
              "2"="C",
              "3"="T",
              "(1 3)"="LT?",
              #"(0 1)"="HL?",
              "?"="?",
              "..."="...")
st_colors = c("red","forestgreen","deepskyblue","darkblue","magenta","black","gray")
names(st_colors) = st_lbl
    

tree_fn = paste(out_fp, out_str, ".ase.tre", sep="")

#st_lbl = NULL

#summary_statistic = "MAP"
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
                      pie_diameter=0.065,
                      pie_nudge_x=0.2,
                      pie_nudge_y=0.2,
                      alpha=1)
