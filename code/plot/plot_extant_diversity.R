
library(stringr)
library(ggmap)
library(reshape2)
library(dplyr)
library(maps)
library(Cairo)
library(ggforce)
library(oce)

clean_str = function(s) {
    s = str_remove_all(s, pattern = "[()]")
    s = unlist( str_split(s, pattern=" ") )
    return(s)
}

get_bg_state = function(s) {
    if      (s=="1")  return(c(1))
    else if (s=="2")  return(c(2))
    else if (s=="3")  return(c(3))
    else if (s=="4")  return(c(4))
    else if (s=="5")  return(c(5))
    else if (s=="6")  return(c(6))
    else if (s=="A")  return(c(1,2))
    else if (s=="B")  return(c(1,3))
    else if (s=="C")  return(c(2,3))
    else if (s=="D") return(c(1,4))
    else if (s=="E") return(c(2,4))
    else if (s=="F") return(c(3,4))
    else if (s=="G") return(c(1,5))
    else if (s=="H") return(c(2,5))
    else if (s=="I") return(c(3,5))
    else if (s=="J") return(c(4,5))
    else if (s=="K") return(c(1,6))
    else if (s=="L") return(c(2,6))
    else if (s=="M") return(c(3,6))
    else if (s=="N") return(c(4,6))
    else if (s=="O") return(c(5,6))
}

fp= "/Users/mlandis/projects/vib_div/"
dat_fp = paste(fp, "data/", sep="")
plot_fp = paste(fp, "code/plot/", sep="")
bg_fn = paste(dat_fp, "vib_bg_raw.txt", sep="")
biome_fn = paste(dat_fp, "vib_biomes_raw.txt", sep="")
taxa_fn = paste(dat_fp, "viburnum.taxa.tsv", sep="")
plot_fn = paste(plot_fp, "vib_map.pdf", sep="")

col_biome_fn = paste(plot_fp, "biome_colors.n4.txt",sep="")
col_bg_fn = paste(plot_fp, "range_colors.n6.txt",sep="")
col_biome = read.table(col_biome_fn,sep=",",header=T, stringsAsFactors=F)
col_bg = read.table(col_bg_fn,sep=",",header=T, stringsAsFactors=F)

dat_bg = readLines(bg_fn)
dat_biome = readLines(biome_fn)
taxa = read.csv(taxa_fn, sep="\t", stringsAsFactors=F)$taxon
taxa = taxa[1:163]

n_biomes=4
n_areas=6
n_taxa = length(taxa)
m = matrix(0, nrow=n_areas, ncol=n_biomes)
for (i in 1:length(dat_bg)) {
    s_bg = clean_str( dat_bg[i] )
    s_biome = clean_str( dat_biome[i] )
    if (s_bg[1] != s_biome[1]) {
        stop("Name mismatch!")
    }
    s_name  = s_bg[1]
    x_bg    = s_bg[2:length(s_bg)]
    n_bg    = 1/length(x_bg)
    y_bg = c()
    for (i in 1:length(x_bg)) {
        y_bg = c(y_bg, as.numeric(get_bg_state(x_bg[i])))
    }
    #y_bg    = as.numeric( get_bg_state(x_bg) )
    y_biome = as.numeric( s_biome[2:length(s_biome)] ) + 1
    n_biome = 1/length(y_biome)
    
    m[y_bg,y_biome] = m[y_bg,y_biome] + n_biome*n_bg
}
rownames(m)=c("SE Asia", "E Asia", "Europe", "N. America", "C. America", "S. America")
colnames(m)=c("Tropical", "Warm Temp.", "Cloud", "Cold Temp.")
mm = melt(m)
colnames(mm) = c("area","biome","count")
mm$lat = 0
mm$long = 0

#mm$area = as.vector( mm$area )
#mm$biome = as.vector( mm$biome )
mm$count = as.vector( as.integer(mm$count) )
mm$count[mm$count==0] = NA

rel_offset = 10
biome_offset = list( "Tropical"=c(0,0)*rel_offset,
                     "Warm Temp."=c(1,0)*rel_offset,
                     "Cloud"=c(2,0)*rel_offset,
                     "Cold Temp."=c(3,0)*rel_offset)

area_coords = list( "SE Asia"=c(100,10),
                    "E Asia"=c(95,45),
                    "Europe"=c(20,50),
                    "N. America"=c(-90,45),
                    "C. America"=c(-90,15),
                    "S. America"=c(-60,-15))

circle_dat = data.frame( matrix(unlist(area_coords),ncol=2,byrow=T) )
circle_dat$area = names(area_coords)
circle_dat$radius = c( 15, 23, 18, 20, 12, 15 )
circle_dat$area = factor( circle_dat$area, levels = names(area_coords), ordered=T)
colnames(circle_dat)=c("long","lat","area", "radius")

circle_dat$long[1] = 107
circle_dat$long[2] = 104
circle_dat$lat[2] = 48
circle_dat$long[3] = 26
circle_dat$lat[4] = 50
circle_dat$long[5] = -83

for (i in 1:nrow(mm)) {
    a = mm$area[i]
    b = mm$biome[i]
    cat(a,b,"\n")
    mm$long[i] = area_coords[[ a ]][1] #+ biome_offset[[ b ]][1]
    mm$lat[i] = area_coords[[ a ]][2]# + biome_offset[[ b ]][2]
}

dx = 10
dy = 10
# SE Asia
mm$long[ c(1,7) ] = mm$long[ c(1,7) ] + c(0,1.36)*dx
# E Asia
mm$long[ c(2,8,20) ] = mm$long[ c(2,8,20) ] + c(0.5,0,2.75)*dx
mm$lat[ c(2,8,20) ] = mm$lat[ c(2,8,20) ] + c(-0.87,0,0)*dx
# Eur
mm$long[ c(9,21) ] = mm$long[ c(9,21) ] + c(0,1.45)*dx
# N Am
mm$long[ c(10,22) ] = mm$long[ c(10,22) ] + c(0,1.67)*dx
# C Am
mm$long[ c(11,17) ] = mm$long[ c(11,17) ] + c(0,1.25)*dx

m_fossil = data.frame(
    locality=c("Northwest Territories","Paris Basin","Iceland","British Columbia","Florissant"),
    clade=c("Valvatotinus","Valvatotinus","Valvatotinus","Viburnum","Viburnum"),
    age=c(48,45,15,48,35),
    lat =c( 79, 49, 65, 50, 39),
    long=c(-91,  2,-18,-120, -105),
    pollen_type=c("Ib","Ib","Ic","Ia","Ia"))
    


# plot limits
xlim = c(-150,150)
ylim = c(-63,90)


mp1 <- fortify(map(fill=TRUE, plot=FALSE))
mp2 <- mp1
mp2$long <- mp2$long + 360
mp2$group <- mp2$group + max(mp2$group) + 1
mp <- rbind(mp1, mp2)

#worldmap = map_data("world")
#setnames(worldmap, c("X","Y","PID","POS","region","subregion"))
#worldmap = clipPolys(worldmap, xlim=xlim,ylim=ylim, keepExtra=TRUE)


# good map
#worldMap <- getMap()
#world.points <- fortify(worldMap)
#world.points$region <- world.points$id
#world.df <- world.points[,c("long","lat","group", "region")]


biome_col = col_biome$color[1:4]
names(biome_col) = col_biome$name[1:4]
bg_col = col_bg$color[1:6]
names(bg_col) = names(area_coords)

fossils_only = !T
include_fossils = T

pp = ggplot() #zz, aes(x=long, y=lat, group=group))
pp = pp + geom_polygon(data=mp1, aes(x=long, y=lat, group=group), fill="lightgray")
pp = pp + coord_map(projection="mollweide",xlim=xlim, ylim=ylim)
pp = pp + theme_void()  # Remove ugly grey background

if (!fossils_only) {


    pp = pp + geom_circle(data=circle_dat, aes(x0=long, y0=lat, r=radius, fill=area), alpha=0.5, colour=NA )
    pp = pp + scale_fill_manual( values=bg_col )


    pp = pp + geom_point(data=mm, aes(x=long, y=lat, colour=biome, size=count))
    pp = pp + geom_text(data=mm, aes(x=long, y=lat, label=count), colour="white")
    pp = pp + scale_colour_manual( values=biome_col )

    pp = pp + guides(fill=guide_legend(title="Areas"), colour=guide_legend(title="Biomes", override.aes=list(size=5)), size=FALSE)
    pp = pp + scale_size(range = c(5, 20)*0.8)
}

if (include_fossils) {
    pp = pp + geom_point(data=m_fossil, aes(x=long, y=lat), size=10, shape=42)
}
#pp = pp + geom_text(data=m_fossil, aes(x=long, y=lat, label=pollen_type), nudge_x=0, nudge_y=5, family="Courier")

CairoPDF(plot_fn, height=4, width=8)
pp
dev.off()

