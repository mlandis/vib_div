# load libraries
library(Cairo)
library(cowplot)
library(ggplot2)
#library(grid)
#library(scales)
library(data.table)

# load various helper functions
source("vib_div_util.R")

# filename prefix
base_fn = "out.1.t163.f5"

# directories
fp           = "../../"
out_fp       = paste0(fp, "output/", sep="")
plot_fp      = paste0(fp, "code/plot/fig/", sep="")

# filenames
col_bg_fn    = paste0(fp, "code/plot/range_colors.n6.txt")
col_biome_fn = paste0(fp, "code/plot/biome_colors.n4.txt")
plot_fn      = paste0(plot_fp, base_fn, ".lstt.pdf")
bg_fn        = paste0(out_fp, base_fn, ".bg.history.tsv")
biome_fn     = paste0(out_fp, base_fn, ".biome.history.tsv")

# read in color files
bg_colors    = read.csv(col_bg_fn, header=T)
biome_colors = read.csv(col_biome_fn, header=T)

# biogeography plot settings
bg_names = c("SEAs", "EAs", "Eur", "NAm", "CAm", "SAm")
area_cols = as.vector( bg_colors$color[1:6] )
names(area_cols) = c("1","2","3","4","5","6")
bg_label_col = c("white","white","white","white","black","black")

# biome plot settings
biome_names = c("Trop.","Warm","Cloud","Cold")
biome_cols = as.vector(biome_colors$color[1:4] )
names(biome_cols) = c("1","2","3","4")
biome_label_col = c("white","white","black","white")

# data dimensions
n_areas   = length(bg_names)
n_biomes  = length(biome_names)
n_states  = n_areas * n_biomes
max_time  = 90
bin_width = 1
n_bins    = max_time / bin_width
ages      = seq(0.0, max_time, by=bin_width)

# settings to control LSTT accuracy
D_tol     = 0.05 # how close all sampled LSTT frequencies in a bin must be to the true LSTT frequencies
alpha_tol = 0.05 # false-positive rate of failing to satisfy the D_tol level

# apply burnin/thinning if desired
f_burn    = 0.998
thinby    = 1

###################

# read in stochastic mappings
stoch_bg     = read.csv(bg_fn, sep="\t", stringsAsFactors=F)
stoch_bg     = stoch_bg[ stoch_bg$transition_type != "cladogenetic", ]
stoch_biome  = read.csv(biome_fn, sep="\t", stringsAsFactors=F)

# filter out non-events
stoch_bg$transition_time[ stoch_bg$transition_type=="no_change" ] = stoch_bg$branch_start_time[ stoch_bg$transition_type=="no_change" ]
stoch_biome$transition_time[ stoch_biome$transition_type=="no_change" ] = stoch_biome$branch_start_time[ stoch_biome$transition_type=="no_change" ]

# iterations
iterations = unique(stoch_bg$iteration)
n_it       = length(iterations)
n_burn     = max(1, f_burn*length(iterations))
iterations = iterations[n_burn:length(iterations)]
iterations = iterations[ seq(1, length(iterations), by=thinby) ]

# Stage 1: construct combined bg + biome stochastic mappings
stoch_bg_biome = make_stoch_bg_biome(stoch_bg, stoch_biome, iterations)

# Stage 2: create time-binned bg + biome occupancies (main obj. needed for LSTT)
state_bins= make_bg_biome_bins(stoch_bg_biome, bin_width, iterations, n_areas, n_biomes)

# Stage 3: create plotting table
min_sample = calculate_ske(s=Inf,k=n_states,alpha=alpha_tol,D=D_tol)$n
plot_dat = make_lstt_dat(state_bins, ages, min_sample)
plot_dat_bg = plot_dat$bg
plot_dat_biome = plot_dat$biome

sz=4
pbg=list()
for (i in 1:6) {
    title_i = bg_names[i]
    di = plot_dat_bg[ plot_dat_bg$Area==i & plot_dat_bg$Support > 0, ]
    p2 = ggplot(data = di, aes(fill=Biome, x=age, y=count))
    p2 = p2 + geom_bar(stat="identity",position="fill")
    p2 = p2 + scale_x_continuous("", trans="reverse", limits=c(max_time,0)) #, sec.axis = sec_axis(~ ., name=title_i, breaks=c() ))
    p2 = p2 + scale_fill_manual(values=biome_cols, labels = c("Trop.", "Warm", "Cloud", "Cold"), guide="legend")
    p2 = p2 + scale_alpha_continuous(range=c(0,1), guide="legend")
    p2 = p2 + annotate("text", x=max_time/2, y=0.8, label=title_i, colour=bg_label_col[i], size=sz)

    if (i==1) {
        p2_bg = p2 + guides(alpha = F, fill=guide_legend( title.theme=element_text(size=14),
                                                          label.theme = element_text(size = 12)))
        p_bg_legend = get_legend(p2_bg)
    }
    p2 = p2 + theme(axis.title.x=element_blank(),
                    axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.text.x=element_text(size=12),
                    legend.text=element_text(size=12),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(colour = "black"),
                    legend.position = "none")   
    pbg[[i]] = p2
}

pbg_1 = plot_grid( pbg[[ 1]], pbg[[ 2]], pbg[[ 3]], pbg[[ 4]], pbg[[ 5]], pbg[[6]], ncol=1) 
pbg_2 = add_sub(pbg_1, "Age (Ma)", vpadding=grid::unit(0,"lines"),y=5, x=0.5, vjust=4.5)

pbiome=list()
for (i in 1:4) {
    title_i = biome_names[i]
    di = plot_dat_biome[ plot_dat_biome$Biome==i & plot_dat_biome$Support > 0, ]
    p2 = ggplot(data = di, aes(fill=Area, x=age, y=count))
    p2 = p2 + geom_bar(stat="identity",position="fill")
    p2 = p2 + scale_x_continuous("", trans="reverse", limits=c(max_time,0) )
    p2 = p2 + scale_fill_manual(values=area_cols, labels = c("SEAs", "EAs", "Eur", "NAm", "CAm", "SAm"), guide="legend")
    p2 = p2 + scale_alpha_continuous(range = c(0, 1), guide="legend") #, breaks=seq(0,1,by=.5), limits=seq(0,1,by=0.5))
    p2 = p2 + annotate("text", x=max_time/2, y=0.8, label=title_i, colour=biome_label_col[i], size=sz)
    
    if (i == 1) {
        p2_biome = p2 + guides(alpha = F, fill=guide_legend( title.theme=element_text(size=14),
                                                             label.theme = element_text(size=12)))
        p_biome_legend = get_legend(p2_biome)
        p2_freq = p2 + scale_fill_manual(values=area_cols, labels = c("SEAs", "EAs", "Eur", "NAm", "LAm"), guide=FALSE)
        p_freq_legend = get_legend(p2_freq)
    }
    
    p2 = p2 + theme(axis.title.x=element_blank(),
                    axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.text.x=element_text(size=12),
                    legend.text=element_text(size=18),
                    legend.title=element_text(size=18),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(colour = "black"),
                    legend.position = "none")      
    pbiome[[i]] = p2
}

plegend_left = plot_grid(NULL, p_bg_legend, NULL, ncol=1)
plegend_right = plot_grid(NULL, p_biome_legend, NULL, ncol=1)

pbiome_1 = plot_grid( NULL, pbiome[[ 1]], pbiome[[ 2]], pbiome[[ 3]], pbiome[[ 4]], NULL,  ncol=1)
pbiome_2 = add_sub(pbiome_1, "Age (Ma)", vpadding=grid::unit(0,"lines"),y=5, x=0.5, vjust=4.5)

pboth = plot_grid( plegend_left, pbg_2, pbiome_2, plegend_right, ncol=4, rel_heights=c(1,6,4,2), rel_widths=c(1,3,3,1) )

CairoPDF(  file=plot_fn, width=7, height=7)
print(pboth)
dev.off()

