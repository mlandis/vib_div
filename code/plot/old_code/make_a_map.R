library(ggplot2)
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
world2 = subset(world, continent != "Antarctica")
#world
coast <- ne_coastline(scale = "medium", returnclass = "sf")
bb    <- ne_download(type = "wgs84_bounding_box", category = "physical",
                  returnclass = "sf")
#coast2 = subset(coast, continent != "Antarctica")

p = ggplot() 
#p = p + geom_sf(data = bb, col = "grey20", fill = "transparent")
p = p + geom_sf(data = world2, size=0, colour="antiquewhite",fill="antiquewhite")
p = p + geom_sf(data=coast, colour="black", size=0.3)
#p = p + coord_sf(crs = "+proj=laea +lat_0=89 +lon_0=300 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs ")
p = p + coord_sf(crs = "+proj=laea +lat_0=89 +lon_0=300 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs ")
p = p + theme_minimal()
p = p + theme(axis.text = element_blank())
p = p + theme(panel.grid.major = element_line(colour = 'transparent'))
p

fp = "/Users/mlandis/projects/vib_div/"
plot_fn = paste(fp, "code/plot/fig1_vib_map.pdf", sep="")


pdf(plot_fn, height=7, width=7)
print(p)
dev.off()
