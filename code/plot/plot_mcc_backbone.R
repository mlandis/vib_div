library(ggplot2)
library(ggtree)
library(ape)

# written by Sebastian Duchene
allnode.times <- function(phylo, tipsonly = FALSE){
    di.phylo <- dist.nodes(phylo)
    root.phylo <- phylo$edge[, 1][!(phylo$edge[, 1] %in% phylo$edge[, 2])][1]
    phylo.depth <- max(di.phylo[as.numeric(colnames(di.phylo)) == root.phylo, ])
    node.times <- phylo.depth - di.phylo[as.numeric(rownames(di.phylo)) == root.phylo, ]
    if(tipsonly == TRUE){
    	node.times <- phylo.depth - di.phylo[as.numeric(rownames(di.phylo)) == root.phylo, 1:length(phylo$tip.label)]
    }
    return(node.times)
}

add_epoch_times <- function( p, dy=4 ) {
    max_x = max(p$data$x)
    max_y = max(p$data$y)
    epoch_names = c("Late\nCretaceous","Paleogene","Early\nEocene","Mid/Late\nEocene","Oligocene","Early\nMiocene","Mid/Late\nMiocene","Recent")
    x_pos = max_x-c(90, 65, 56, 48, 33.9, 23, 16, 5.3, 0)
    x_pos_mid = (x_pos[1:(length(x_pos)-1)] + x_pos[2:length(x_pos)])/2 
    #x_names = max_x - (x_pos[1:(length(x_pos)-1)] + x_pos[2:length(x_pos)])/2 
    
    t = "dotted"
    for (k in 2:length(x_pos)) {
        p = p + geom_vline(xintercept=x_pos[k], color="pink2", linetype=t)
    }
    for (k in 1:length(epoch_names)) {
        p = p + geom_text( label=epoch_names[k], x=x_pos_mid[k], y=0-dy, hjust=0.5, size=2)
        #text(x=x_names[k]+1, y=-2, labels=epoch_names[k], cex=0.35, adj=0)
    }
    return(p)

}

# default args
fp = "/Users/mlandis/projects/vib_div/"
fn = "out.1.t163.f5.radseq_only.con.tre"

# update args
args = commandArgs(trailingOnly=TRUE)
if (length(args)>0) {
    fn = args[1]
}

# filepaths
phy_fp = paste(fp, "output/", sep="")
plot_fp = paste(fp, "code/plot/fig/mcc/", sep="")

# read tree
phy_fn = paste(phy_fp, fn, sep="")
plot_fn = paste(plot_fp, fn, ".pdf", sep="")

cat("Reading \"",phy_fn,"\"\n",sep="")

# read the tree in various formats
phy = read.beast(phy_fn)
phy2 = read.nexus(phy_fn)

# subtract the MCC age from the HPD range to get lengths
bt = allnode.times(phy2)
heights = bt[ as.numeric(phy@stats$node) ]
phy@stats$height = heights
heights_hpd = phy@stats$age_0.95_HPD
for (i in 1:length(heights_hpd)) {
    hpd = as.numeric(heights_hpd[[i]])
    h = heights[i]
    if (!is.na(hpd)) {
        h_tmp = c( 2*h-hpd[2], 2*h-hpd[1] ) 
        heights_hpd[[i]] = h_tmp
    }
}
phy@stats$height_0.95_HPD = heights_hpd
phy@fields = c(phy@fields, "height", "height_0.95_HPD")





# make tree object, add plottable features
pp = ggtree(phy,  ndigits=2, branch.length="none" )
pp$data$x = pp$data$x + 8 #(90 - max(pp$data$x))
pp$data$posterior_class = NA
pp$data$posterior_class[ which(pp$data$posterior > 0.95) ] = ">0.95"
pp$data$posterior_class[ which(pp$data$posterior<=0.95&pp$data$posterior>0.75) ] = ">0.75"
pp$data$posterior_class[ which(pp$data$posterior<=0.75&pp$data$posterior>0.5) ] = ">0.50"
pp$data$posterior_class[ which(pp$data$posterior<=0.5) ] = "<0.50"
pp$data$posterior_class = factor(pp$data$posterior_class, levels=c(">0.95",">0.75", ">0.50", "<0.50"))
pp$data$label = sapply(pp$data$label, function(x) { gsub("_"," ",x) })
pp$data$label = sapply(pp$data$label, function(x) { gsub("subsp. ","",x) })
dy = 4
#ylim = c(0, max(pp$data$y)+dy)

#pp = pp + ylim(ylim)
#pp = add_epoch_times(pp, dy)
#pp = pp + geom_range(range='height_0.95_HPD', color='dodgerblue', alpha=.6, size=1)
pp = pp + geom_tiplab(size=1.6, offset=0.25)
pp = pp + geom_nodepoint( size=2, color="black")
pp = pp + geom_nodepoint(aes(color=posterior_class), size=1.5)
pp = pp + scale_color_manual(values=c("#000000","#666666","#BBBBBB","#EEEEEE"), name="Posterior")
#pp = pp + scale_x_continuous(breaks = seq(0,90,10), labels = rev(seq(0,90,10)))
#pp = pp + theme_tree2()
pp = pp + coord_cartesian(xlim = c(8.5,40), expand=TRUE)
#pp = pp + labs(x="Age (Ma)")
pp = pp + theme(legend.position=c(0.05,0.9), axis.line = element_line(colour = "black"))


dx = 1
clade_mtx_1 = matrix(ncol=3, byrow=T,
                    data=c(
                           "Urceolata",    0, getMRCA(phy@phylo,tip=c("V_urceolatum","V_taiwanianum")),
                           "Lentago",      0, getMRCA(phy@phylo,tip=c("V_nudum","V_obovatum")),
                           "Euviburnum",   0, getMRCA(phy@phylo,tip=c("V_lantana","V_cotinifolium")),
                           "Pseudotinus",  0, getMRCA(phy@phylo,tip=c("V_lantanoides","V_nervosum")),
                           "Solenotinus",  0, getMRCA(phy@phylo,tip=c("V_foetens","V_sieboldii")),
                           "Lutescentia",  0, getMRCA(phy@phylo,tip=c("V_lutescens","V_plicatum")),
                           "Tinus",        0, getMRCA(phy@phylo,tip=c("V_tinus","V_cinnamomifolium")),
                           "Mollotinus",   0, getMRCA(phy@phylo,tip=c("V_molle","V_ellipticum")),
                           "Dentata",      0, getMRCA(phy@phylo,tip=c("V_dentatum","V_scabrellum")),
                           "Oreintotinus", 0, getMRCA(phy@phylo,tip=c("V_loeseneri","V_stipitatum")),
                           "Opulus",       0, getMRCA(phy@phylo,tip=c("V_opulus","V_edule")),
                           "Sambucina",    0, getMRCA(phy@phylo,tip=c("V_sambucinum","V_beccarii")),
                           "Coriacea",     0, getMRCA(phy@phylo,tip=c("V_coriaceum","V_cylindricum")),
                           "Lobata",       0, getMRCA(phy@phylo,tip=c("V_kansuense","V_acerifolium")),
                           "Succotinus",   0, getMRCA(phy@phylo,tip=c("V_foetidum","V_erosum")),
                           "Punctata",     0, which( phy@phylo$tip.label=="V_punctatum") #getMRCA(phy@phylo,tip=c("V_punctatum","V_punctatum"))
                    ))
clade_df_1 = data.frame(clade_mtx_1, stringsAsFactors=F)
colnames(clade_df_1) =c("name","depth","index")
clade_df_1$depth=as.numeric(clade_df_1$depth)
clade_df_1$index=as.numeric(clade_df_1$index)
    
for (i in 1:nrow(clade_df_1)) {
    drow = clade_df_1[i,]
    pp = pp + geom_cladelabel(node=drow$index,
                              label=drow$name,
                              align=T,
                              fontsize=3,
                              #offset.text=.2,
                              offset=6.5,
                              hjust='left',
                              barsize=1., extend=0.17)
}

clade_mtx_2 = matrix(ncol=3, byrow=T,
                    data=c(
                           "Valvatotinus",   1, getMRCA(phy@phylo,tip=c("V_lantana","V_punctatum")),
                           "Oreinodentinus", 1, getMRCA(phy@phylo,tip=c("V_dentatum","V_obtectum")),
                           "Porphyrotinus",  1, getMRCA(phy@phylo,tip=c("V_obtectum","V_molle")),
                           "Imbricotinus",   3, getMRCA(phy@phylo,tip=c("V_erosum","V_molle")),
                           "Nectarotinus",   4, getMRCA(phy@phylo,tip=c("V_tinus","V_molle")),
                           "Laminotinus",    1, getMRCA(phy@phylo,tip=c("V_sambucinum","V_erosum")),
                           "Crenotinus",     1, getMRCA(phy@phylo,tip=c("V_lutescens","V_foetens")),
                           "Amplicrenotinus",2, getMRCA(phy@phylo,tip=c("V_lutescens","V_amplificatum"))
                    ))
clade_df_2 = data.frame(clade_mtx_2, stringsAsFactors=F)
colnames(clade_df_2) =c("name","depth","index")
clade_df_2$depth=as.numeric(clade_df_2$depth)
clade_df_2$index=as.numeric(clade_df_2$index)


for (i in 1:nrow(clade_df_2)) {
    drow = clade_df_2[i,]
    dx=1.5
    #pp = pp + geom_cladelabel(node=clade_idx, label=drow[1], angle=270, align=T, offset.text=.2, offset=2.5+as.numeric(drow[2])*dx, hjust='center', barsize=1.25)
    pp = pp + geom_cladelabel(node=drow$index,
                              label=drow$name,
                              align=T,
                              fontsize=3,
                              angle=270,
                              offset.text=.4,
                              offset=12+drow$depth*dx,
                              hjust='center',
                              barsize=1., extend=0.17)
}



print(pp)




cat("Writing \"",plot_fn,"\"\n",sep="")
pdf(plot_fn, width=4, height=10)
print(pp)
dev.off()

