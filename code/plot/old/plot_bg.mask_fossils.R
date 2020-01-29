stable = !T
if (stable) {
    library(RevGadgets)
} else {
    source("/Users/mlandis/projects/RevGadgets/R/plot_ancestral_states.R")
}
library(ggtree)

out_fp = "/Users/mlandis/projects/vib_div/"

#out_str = "output/tinus/out.1.viburnum.t163.f5.bg"
out_str = "output/tinus/out.1.viburnum.t163.f5.mask_fossil_states.bg"

st_lbl = list("1"="SE Asia",
           "2"="E Asia",
           "3"="Eur",
           "4"="N Am",
           "5"="Lat Am",
           "6"="SE Asia + E Asia",
           "7"="SE Asia + Eur",
           "8"="E Asia + Eur",
           "9"="SE Asia + N Am",
           "10"="E Asia + N Am",
           "11"="Eur + N Am",
           "12"="SE Asia + Lat Am",
           "13"="E Asia + Lat Am",
           "14"="Eur + Lat Am",
           "15"="N Am + Lat Am",
           "..."="...")

st_colors = c("magenta",
              "red",
              "green",
              "gold",
              "blue",
              "coral",
              "pink",
              "darkgoldenrod3",
              "tan",
              "brown",
              "limegreen",
              "cadetblue",
              "purple",
              "deepskyblue",
              "seagreen",
              "gray")
    
names(st_colors) = st_lbl

tree_fn = paste(out_fp, out_str, ".ase.tre", sep="")
summary_statistic = "PieRange"
#summary_statistic = "MAPRange"
#summary_statistic = "PieState"
#summary_statistic = "MAP"
vv=plot_ancestral_states(tree_file=tree_fn,
                      include_start_states=T,
                      shoulder_label_size=0,
                      summary_statistic=summary_statistic,
                      state_labels=st_lbl,
                      state_colors=st_colors,
                      node_label_size=0,
                      node_size_range=c(0.5,2),
                      node_label_nudge_x=0.5,
                      tip_node_size=0.75,
                      tip_label_size=1.0,
                      tip_label_offset=0.5,
                      xlim_visible=c(0,70),
                      shoulder_label_nudge_x=-0.1,
                      show_posterior_legend=T,
                      pie_diameter=0.07,
                      pie_nudge_x=0.2,
                      pie_nudge_y=0.2,
                      alpha=1)
ggsave( paste(out_str,".pdf",sep="") )
