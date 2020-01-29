library(ggplot2)
library(ggtree)
library(ape)
library(cowplot)

source("vib_div_util.R")

# filepaths
fp       = "../../"
phy_fp   = paste(fp, "data/", sep="")
plot_fp  = paste(fp, "code/plot/fig/", sep="")
fn_a     = "radseq_stage1.tre"
fn_b     = "radseq_cpdna_stage2.tre"
phy_a_fn = paste(phy_fp, fn_a, sep="")
phy_b_fn = paste(phy_fp, fn_b, sep="")
plot_fn = paste(plot_fp, "fig2_topology.pdf", sep="")

# read tree
#cat("Reading \"",phy_b_fn,"\"\n",sep="")

# read the tree in various formats
phy_a = read.beast(phy_a_fn)
phy_a@phylo = fix_vib_tip(phy_a@phylo)
phy_a@data$posterior = phy_a@data$posterior / 100
phy_a@data$posterior[ phy_a@data$posterior == 1 ] = NA


phy_b = read.beast(phy_b_fn)
phy_b@phylo = fix_vib_tip(phy_b@phylo)
phy_b@data$posterior[ phy_b@data$posterior == 1 ] = NA

# tip colors
lbl_a = phy_a@phylo$tip.label
lbl_b = phy_b@phylo$tip.label
tip_col = rep("#AAAAAA", length(lbl_b))
tip_col[ match(lbl_a,lbl_b) ] = "black"

# make tree object, add plottable features
print("Fig 2A")
pp_a = ggtree(phy_a, branch.length="none" )
pp_a$data$x = pp_a$data$x# + 8
pp_a$data$posterior_class = NA
pp_a$data$posterior_class[ which(pp_a$data$posterior==0.99) ] = "0.99"
pp_a$data$posterior_class[ which(pp_a$data$posterior<=0.98&pp_a$data$posterior>0.95) ] = ">0.95"
pp_a$data$posterior_class[ which(pp_a$data$posterior<=0.95&pp_a$data$posterior>0.75) ] = ">0.75"
pp_a$data$posterior_class[ which(pp_a$data$posterior<=0.75&pp_a$data$posterior>0.5) ] = ">0.50"
pp_a$data$posterior_class[ which(pp_a$data$posterior<=0.5) ] = "<0.50"
pp_a$data$posterior_class = factor(pp_a$data$posterior_class, levels=c("0.99",">0.95",">0.75", ">0.50", "<0.50"))
pp_a$data$label = sapply(pp_a$data$label, function(x) { gsub("_"," ",x) })
pp_a$data$label = sapply(pp_a$data$label, function(x) { gsub("subsp. ","",x) })
dy = 4
pp_a = pp_a + geom_tiplab(size=1.6, offset=0.25)
pp_a = pp_a + geom_nodepoint( data=pp_a$data[ !is.na(pp_a$data$posterior), ], size=2, color="black")
pp_a = pp_a + geom_nodepoint( data=pp_a$data[ !is.na(pp_a$data$posterior), ], aes(color=posterior_class), size=1.5)
col_a = c("#000000","#666666","#999999","#BBBBBB","#EEEEEE")
names(col_a) = levels(pp_a$data$posterior_class)
pp_a = pp_a + scale_color_manual(values=col_a, name="Bootstrap")
pp_a = pp_a + coord_cartesian(xlim = c(-5,25), ylim=c(0,155+2)-(155-120), expand=TRUE)
pp_a = pp_a + theme(legend.position=c(0.1,0.86), axis.line = element_line(colour = "black"))

# make tree object, add plottable features
print("Fig 2B")
pp_b = ggtree(phy_b, branch.length="none")
pp_b$data$x = pp_b$data$x# + 8
pp_b$data$posterior_class = NA
pp_b$data$posterior_class[ which(pp_b$data$posterior>=0.985) ] = "0.99"
pp_b$data$posterior_class[ which(pp_b$data$posterior<0.985&pp_b$data$posterior>0.95) ] = ">0.95"
pp_b$data$posterior_class[ which(pp_b$data$posterior<=0.95&pp_b$data$posterior>0.75) ] = ">0.75"
pp_b$data$posterior_class[ which(pp_b$data$posterior<=0.75&pp_b$data$posterior>0.5) ] = ">0.50"
pp_b$data$posterior_class[ which(pp_b$data$posterior<=0.5) ] = "<0.50"
pp_b$data$posterior_class = factor(pp_b$data$posterior_class, levels=c("0.99",">0.95",">0.75", ">0.50", "<0.50"))
pp_b$data$label = sapply(pp_b$data$label, function(x) { gsub("_"," ",x) })
pp_b$data$label = sapply(pp_b$data$label, function(x) { gsub("subsp. ","",x) })
dy = 4
pp_b = pp_b + geom_tiplab( color=tip_col, size=1.6, offset=0.25)
pp_b = pp_b + geom_nodepoint( data=pp_b$data[ !is.na(pp_b$data$posterior), ], size=2, color="black")
pp_b = pp_b + geom_nodepoint( data=pp_b$data[ !is.na(pp_b$data$posterior_class), ], aes(color=posterior_class), size=1.5)
col_b = c("#000000","#666666","#999999","#BBBBBB","#EEEEEE")
names(col_b) = levels(pp_b$data$posterior_class)
pp_b = pp_b + scale_color_manual(values=col_b, name="Posterior")
pp_b = pp_b + coord_cartesian(xlim = c(0,30), ylim=c(0,155+2), expand=TRUE)
pp_b = pp_b + theme(legend.position=c(0.09,0.86), axis.line = element_line(colour = "black"))


pg = plot_grid( pp_a, pp_b,
                nrow=1, ncol=2,
                rel_widths = c(13,17),
                labels=c("(A) Stage 1 topology\n", "          (B) Stage 2 topology\n(only sequenced taxa shown)"),
                hjust=-0.25,
                label_size=12)
pg

pdf(plot_fn, height=11, width=6)
print(pg)
dev.off()



if (FALSE) {
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

}
