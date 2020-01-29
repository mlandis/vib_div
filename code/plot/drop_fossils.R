library(geiger)

fp = "../../"
data_fp = paste(fp, "data/", sep="")
out_fp = paste(fp, "output/", sep="")

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    fn = "out.1.t163.f5"
} else if (length(args)==1) {
    fn = args[1]
}

fossil_tips = c("Porphyrotinus_CO","Viburnum_BC", "Valvatotinus_NWT", "Valvatotinus_IS", "Valvatotinus_PB")
unsequenced_tips = c("V_junghunii","V_garrettii","V_shweliense","V_wardii","V_hondurense","V_longiradiatum","V_dalzielii","V_chunii","V_hainanense","V_tengyuehense")

# drop fossils
dat_phy_no_fossil=read.csv(paste(out_fp,fn,".trees",sep=""), sep="\t", stringsAsFactors=F)
for (i in 1:nrow(dat_phy_no_fossil)) {
    phy = read.tree( text=dat_phy_no_fossil$tree[i] )
    phy = drop.tip(phy, fossil_tips)
    #print(phy)
    dat_phy_no_fossil$tree[i] = write.tree(phy, file="")
}

write.table(dat_phy_no_fossil, file=paste(out_fp,fn,".no_fossil.trees",sep=""), sep="\t", quote=F, row.names=F)


# drop fossils/unsequenced
dat_phy_no_fossil_unsequenced=read.csv(paste(out_fp,fn,".trees",sep=""), sep="\t", stringsAsFactors=F)
for (i in 1:nrow(dat_phy_no_fossil_unsequenced)) {
    phy = read.tree( text=dat_phy_no_fossil_unsequenced$tree[i] )
    phy = drop.tip(phy, c(fossil_tips, unsequenced_tips))
    #print(phy)
    dat_phy_no_fossil_unsequenced$tree[i] = write.tree(phy, file="")
}

write.table(dat_phy_no_fossil_unsequenced, file=paste(out_fp,fn,".no_fossil_unsequenced.trees",sep=""), sep="\t", quote=F, row.names=F)

# drop non-RADseq
phy_radseq = read.tree( paste(data_fp, "viburnum.backbone.tre", sep="" ) )[[1]]
dat_phy_radseq_only=read.csv(paste(out_fp,fn,".trees",sep=""), sep="\t", stringsAsFactors=F)
phy_all = read.tree( text=dat_phy_radseq_only$tree[1] )
non_radseq_tips = phy_all$tip.label[ !(phy_all$tip.label %in% phy_radseq$tip.label) ]

#drop_tips = c("Porphyrotinus_CO","Viburnum_BC", "Valvatotinus_NWT", "Valvatotinus_IS", "Valvatotinus_PB")
for (i in 1:nrow(dat_phy_radseq_only)) {
    phy = read.tree( text=dat_phy_radseq_only$tree[i] )
    phy = drop.tip(phy, non_radseq_tips)
    #print(phy)
    dat_phy_radseq_only$tree[i] = write.tree(phy, file="")
}

write.table(dat_phy_radseq_only, file=paste(out_fp,fn,".radseq_only.trees",sep=""), sep="\t", quote=F, row.names=F)
