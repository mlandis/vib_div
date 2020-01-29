library(ape)
library(cowplot)
library(deeptime)
library(gginnards)
library(ggplot2)
library(ggtree)
library(RevGadgets)

make_height_hpd = function(phy_fn) {
    
    phy = read.beast(phy_fn)
    phy2 = read.nexus(phy_fn)
    
    # subtract the MCC age from the HPD range to get lengths
    bt = allnode.times(phy2)
    heights = bt[ as.numeric(phy@data$node) ]
    phy@data$height = heights
    heights_hpd = phy@data$age_0.95_HPD
    for (i in 1:length(heights_hpd)) {
        hpd = as.numeric(heights_hpd[[i]])
        h = heights[i]
        if (!is.na(hpd)) {
            h_tmp = c( h-hpd[2], h-hpd[1] )
            heights_hpd[[i]] = h_tmp
        }
    }
    phy@data$height_0.95_HPD = heights_hpd
    
    return(heights_hpd)
}

fix_vib_tip = function(phy) {
    lbl = phy$tip.label
    lbl = as.vector(sapply(lbl, function(x) { strsplit(x, "_")[[1]][2] } ))
    lbl[ lbl=="rigidum" ] = "rugosum"
    lbl[ lbl=="stellato" ] = "stellato-tomentosum"
    phy$tip.label = lbl
    return(phy)
}

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


add_epoch_times <- function( p, max_age, dy_bars, dy_text ) {
    
    max_x = max(p$data$x)
    max_y = max(p$data$y)
    epoch_names = c("Late\nCretaceous","Paleocene","Early\nEocene","Mid/Late\nEocene","Oligocene","Early\nMiocene","Mid/Late\nMiocene","Recent")
  

    x_pos = max_x-c(max_age, 65, 56, 48, 33.9, 23, 16, 5.3, 0)
    y_pos = rep(max_y, length(x_pos))
    x_pos_mid = ( x_pos[1:(length(x_pos)-1)] + x_pos[2:length(x_pos)] ) / 2 

    for (k in 2:(length(x_pos))) {
        box_col = "gray92"
        if (k %% 2 == 0) box_col = "white"
        box = geom_rect( xmin=x_pos[k-1], xmax=x_pos[k], ymin=dy_bars, ymax=y_pos[k], fill=box_col )
        p = append_layers(p, box, position = "bottom")
    }
    for (k in 1:length(epoch_names)) {
        p = p + annotate( geom="text", label=epoch_names[k], x=x_pos_mid[k], y=dy_text, hjust=0.5, size=3.25)
    }
    return(p)

}

insertLayer <- function(P, after=0, ...) {
  #  P     : Plot object
  # after  : Position where to insert new layers, relative to existing layers
  #  ...   : additional layers, separated by commas (,) instead of plus sign (+)

      if (after < 0)
        after <- after + length(P$layers)

      if (!length(P$layers))
        P$layers <- list(...)
      else 
        P$layers <- append(P$layers, list(...), after)

      return(P)
}

make_vib_clade_mtx = function(phy) {
    
    # set up clade names
    clade_mtx = matrix(NA, ncol=4, nrow=0)
    clade_mtx = rbind( clade_mtx, c( "Succotinus",    getMRCA(phy, c("erosum","betulifolium")),       2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Lobata",        getMRCA(phy, c("orientale","kansuense")),       2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Sambucina",     getMRCA(phy, c("sambucinum","beccarii")),       2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Coriacea",      getMRCA(phy, c("cylindricum","coriaceum")),     2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Opulus",        getMRCA(phy, c("opulus","edule")),              2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Oreinotinus",   getMRCA(phy, c("undulatum","caudatum")),        2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Dentata",       getMRCA(phy, c("dentatum","scabrellum")),       2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Mollotinus",    getMRCA(phy, c("molle","australe")),            2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Tinus",         getMRCA(phy, c("tinus","davidii")),             2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Solenotinus",   getMRCA(phy, c("sieboldii","foetens")),         2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Lutescentia",   getMRCA(phy, c("plicatum","amplifolium")),      2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Euviburnum",    getMRCA(phy, c("chinshanense","cotinifolium")), 2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Lentago",       getMRCA(phy, c("rufidulum","nudum")),           2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Punctata",      getMRCA(phy, c("punctatum","lepidotulum")),     2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Pseudotinus",   getMRCA(phy, c("sympodiale","nervosum")),       2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "Urceolata",     getMRCA(phy, c("urceolatum","taiwanianum")),    2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "V. clemensiae", which(phy$tip.label=="clemensiae"),             2, 1 ))
    clade_mtx = rbind( clade_mtx, c( "V. amplificatum", which(phy$tip.label=="amplificatum"),  2, 1 ))
    
    clade_mtx = rbind( clade_mtx, c( "Laminotinus",   getMRCA(phy, c("erosum","cylindricum")), 3.2, 2 ))
    clade_mtx = rbind( clade_mtx, c( "Porphyrotinus", getMRCA(phy, c("undulatum","australe")), 3.2, 2 ))
    clade_mtx = rbind( clade_mtx, c( "Crenotinus", getMRCA(phy, c("amplifolium","foetens")),   3.2, 2 ))
    clade_mtx = rbind( clade_mtx, c( "Valvatotinus", getMRCA(phy, c("lantana","punctatum")),   3.2, 2 ))
    
    colnames(clade_mtx) = c("name","node","level","barsize")
    clade_df = data.frame(clade_mtx, stringsAsFactors=F)
    clade_df$node = as.numeric(clade_df$node)
    clade_df$level = as.numeric(clade_df$level)
    clade_df$barsize = as.numeric(clade_df$barsize)
    
    return(clade_df)
    
}


make_states = function(label_fn, color_fn, fp="./") {

    # generate colors for ranges
    range_color_list = read.csv(color_fn, header=T, sep=",", colClasses="character")
    
    # get area names
    area_names = unlist(sapply(range_color_list$range, function(y) { if (nchar(y)==1) { return(y) } }))
    
    # get state labels
    state_descriptions = read.csv(label_fn, header=T, sep=",", colClasses="character")
    
    # map presence-absence ranges to area names
    range_labels = sapply(state_descriptions$range[2:nrow(state_descriptions)],
        function(x) {
            present = as.vector(gregexpr(pattern="1", x)[[1]])
            paste( area_names[present], collapse="")
        })

    # map labels to colors 
    range_colors = range_color_list$color[ match(range_labels, range_color_list$range) ]
    
    # generate state/color labels
    idx = 1
    st_lbl = list()
    st_colors = c()
    for (j in 1:(nrow(state_descriptions)-1)) {
        st_lbl[[ as.character(j) ]] = range_labels[j]
        st_colors[j] = range_colors[j]
    }
    st_colors[ length(st_colors)+1 ] = "lightgray"
    st_lbl[["..."]] = "..."   
    
    return( list(state_labels=st_lbl, state_colors=st_colors) )
}



make_tip_colors = function(s="unseq") {
    n_rad = c("cotinifolium","lantana","veitchii","glomeratum","rhytidophyllum","buddleifolium","chinshanense","macrocephalum","congestum","bitchiuense","carlesii","utile","schensianum","punctatum","nudum","lentago","rufidulum","prunifolium","obovatum","elatum","amplificatum","plicatum","hanceanum","amplifolium","pyramidatum","lutescens","sieboldii","farreri","foetens","grandiflorum","taitoense","suspensum","brachybotryum","henryi","subalpinum","erubescens","urceolatum","taiwanianum","lantanoides","furcatum","sympodiale","nervosum","clemensiae","tinus","cinnamomifolium","propinquum","davidii","dentatum","ciliatum","stenocalyx","caudatum","microcarpum","sulcatum","microphyllum","acutifolium","jucundum","lautum","obtusatum","hartwegii","pichinchense","hallii","tinoides","jamesonii","toronis","stipitatum","seemenii","obtectum","divaricatum","opulus","sargentii","coriaceum","cylindricum","hispidulum","vernicosum","glaberrimum","beccarii","ternatum","leiocarpum","inopinatum","sambucinum","mullaha","foetidum","japonicum","wrightii","parvifolium","dilatatum","corylifolium","brachyandrum","erosum","luzonicum","tashiroi","formosanum","fordiae","sempervirens","anamensis","setigerum","phlebotrichum","integrifolium","kansuense","orientale","acerifolium")
    n_all = c("anabaptista","subsessile","tinoides","undulatum","glabratum","jamesonii","lasiophyllum","hallii","pichinchense","divaricatum","stipitatum","obtectum","triphyllum","seemenii","toronis","ayavacense","costaricanum","hartwegii","obtusatum","disjunctum","stellato-tomentosum","villosum","jucundum","lautum","hondurense","acutifolium","sulcatum","venustum","blandum","microphyllum","discolor","ciliatum","microcarpum","caudatum","stenocalyx","loeseneri","dentatum","recognitum","scabrellum","molle","rafinesquianum","bracteatum","australe","ellipticum","CO","BC","brevipes","corylifolium","dilatatum","lancifolium","parvifolium","flavescens","ichangense","adenophorum","dalzielii","wrightii","fordiae","formosanum","chunii","luzonicum","tashiroi","phlebotrichum","setigerum","anamensis","longiradiatum","sempervirens","brachyandrum","erosum","hainanense","japonicum","integrifolium","betulifolium","hupehense","lobophyllum","melanocarpum","mullaha","foetidum","acerifolium","orientale","kansuense","hispidulum","vernicosum","glaberrimum","beccarii","inopinatum","sambucinum","leiocarpum","ternatum","coriaceum","hebanthum","cylindricum","opulus","sargentii","trilobum","edule","koreanum","atrocyaneum","calvum","cinnamomifolium","propinquum","davidii","rugosum","tinus","treleasei","chingii","wardii","yunnanense","corymbiflorum","henryi","tengyuehense","brachybotryum","brevitubum","erubescens","oliganthum","shweliense","foetens","grandiflorum","farreri","suspensum","taitoense","subalpinum","awabuki","odoratissimum","sieboldii","colebrookeanum","garrettii","junghunii","hanceanum","plicatum","lutescens","pyramidatum","amplifolium","amplificatum","bitchiuense","carlesii","congestum","schensianum","utile","macrocephalum","buddleifolium","chinshanense","rhytidophyllum","glomeratum","veitchii","lantana","maculatum","burejaeticum","mongolicum","cotinifolium","elatum","obovatum","prunifolium","rufidulum","lentago","cassinoides","nudum","IS","lepidotulum","punctatum","NWT","PB","furcatum","sympodiale","lantanoides","nervosum","taiwanianum","urceolatum","clemensiae")
    n_unseq = c("hondurense","dalzielii","chunii","longiradiatum","hainanense","wardii","tengyuehense","shweliense","garrettii","junghunii")
    n_fossil = c("NWT","PB","CO","BC","IS")
    n_seq = n_all[ -match(c(n_unseq,n_fossil), n_all) ]
    
    
    col = rep("black", length(n_all))
    names(col) = n_all
    if (s=="unseq") {
        col[ n_unseq ] =  "#AAAAAA" 
    }
    return(col)
}    


plot.igraph.mjl = function (x, axes = FALSE, add = FALSE, xlim = c(-1, 1), ylim = c(-1, 
    1), mark.groups = list(), mark.shape = 1/2, mark.col = rainbow(length(mark.groups), 
    alpha = 0.3), mark.border = rainbow(length(mark.groups), 
    alpha = 1), mark.expand = 15, ...) 
{
    graph <- x
    if (!is_igraph(graph)) {
        stop("Not a graph object")
    }
    params <- igraph:::i.parse.plot.params(graph, list(...))
    vertex.size <- 1/200 * params("vertex", "size")
    label.family <- params("vertex", "label.family")
    label.font <- params("vertex", "label.font")
    label.cex <- params("vertex", "label.cex")
    label.degree <- params("vertex", "label.degree")
    label.color <- params("vertex", "label.color")
    label.dist <- params("vertex", "label.dist")
    labels <- params("vertex", "label")
    shape <- igraph:::igraph.check.shapes(params("vertex", "shape"))
    edge.color <- params("edge", "color")
    edge.width <- params("edge", "width")
    edge.lty <- params("edge", "lty")
    arrow.mode <- params("edge", "arrow.mode")
    edge.labels <- params("edge", "label")
    loop.angle <- params("edge", "loop.angle")
    edge.label.font <- params("edge", "label.font")
    edge.label.family <- params("edge", "label.family")
    edge.label.cex <- params("edge", "label.cex")
    edge.label.color <- params("edge", "label.color")
    elab.x <- params("edge", "label.x")
    elab.y <- params("edge", "label.y")
    arrow.size <- params("edge", "arrow.size")
    arrow.width <- params("edge", "arrow.width")
    curved <- params("edge", "curved")
    if (is.function(curved)) {
        curved <- curved(graph)
    }
    layout <- params("plot", "layout")
    margin <- params("plot", "margin")
    margin <- rep(margin, length = 4)
    rescale <- params("plot", "rescale")
    asp <- params("plot", "asp")
    frame <- params("plot", "frame")
    main <- params("plot", "main")
    sub <- params("plot", "sub")
    xlab <- params("plot", "xlab")
    ylab <- params("plot", "ylab")
    palette <- params("plot", "palette")
    if (!is.null(palette)) {
        old_palette <- palette(palette)
        on.exit(palette(old_palette), add = TRUE)
    }
    arrow.mode <- igraph:::i.get.arrow.mode(graph, arrow.mode)
    maxv <- max(vertex.size)
    if (rescale) {
        layout <- norm_coords(layout, -1, 1, -1, 1)
        xlim <- c(xlim[1] - margin[2] - maxv, xlim[2] + margin[4] + 
            maxv)
        ylim <- c(ylim[1] - margin[1] - maxv, ylim[2] + margin[3] + 
            maxv)
    }
    if (!add) {
        plot(0, 0, type = "n", xlab = xlab, ylab = ylab, xlim = xlim, 
            ylim = ylim, axes = axes, frame = frame, asp = asp, 
            main = main, sub = sub)
    }
    if (!is.list(mark.groups) && is.numeric(mark.groups)) {
        mark.groups <- list(mark.groups)
    }
    mark.shape <- rep(mark.shape, length = length(mark.groups))
    mark.border <- rep(mark.border, length = length(mark.groups))
    mark.col <- rep(mark.col, length = length(mark.groups))
    mark.expand <- rep(mark.expand, length = length(mark.groups))
    for (g in seq_along(mark.groups)) {
        v <- V(graph)[mark.groups[[g]]]
        if (length(vertex.size) == 1) {
            vs <- vertex.size
        }
        else {
            vs <- rep(vertex.size, length = vcount(graph))[v]
        }
        igraph:::igraph.polygon(layout[v, , drop = FALSE], vertex.size = vs, 
            expand.by = mark.expand[g]/200, shape = mark.shape[g], 
            col = mark.col[g], border = mark.border[g])
    }
    el <- as_edgelist(graph, names = FALSE)
    loops.e <- which(el[, 1] == el[, 2])
    nonloops.e <- which(el[, 1] != el[, 2])
    loops.v <- el[, 1][loops.e]
    loop.labels <- edge.labels[loops.e]
    loop.labx <- if (is.null(elab.x)) {
        rep(NA, length(loops.e))
    }
    else {
        elab.x[loops.e]
    }
    loop.laby <- if (is.null(elab.y)) {
        rep(NA, length(loops.e))
    }
    else {
        elab.y[loops.e]
    }
    edge.labels <- edge.labels[nonloops.e]
    elab.x <- if (is.null(elab.x)) 
        NULL
    else elab.x[nonloops.e]
    elab.y <- if (is.null(elab.y)) 
        NULL
    else elab.y[nonloops.e]
    el <- el[nonloops.e, , drop = FALSE]
    edge.coords <- matrix(0, nrow = nrow(el), ncol = 4)
    edge.coords[, 1] <- layout[, 1][el[, 1]]
    edge.coords[, 2] <- layout[, 2][el[, 1]]
    edge.coords[, 3] <- layout[, 1][el[, 2]]
    edge.coords[, 4] <- layout[, 2][el[, 2]]
    if (length(unique(shape)) == 1) {
        ec <- igraph:::.igraph.shapes[[shape[1]]]$clip(edge.coords, el, 
            params = params, end = "both")
    }
    else {
        shape <- rep(shape, length = vcount(graph))
        ec <- edge.coords
        ec[, 1:2] <- t(sapply(seq(length = nrow(el)), function(x) {
            igraph:::.igraph.shapes[[shape[el[x, 1]]]]$clip(edge.coords[x, 
                , drop = FALSE], el[x, , drop = FALSE], params = params, 
                end = "from")
        }))
        ec[, 3:4] <- t(sapply(seq(length = nrow(el)), function(x) {
            igraph:::.igraph.shapes[[shape[el[x, 2]]]]$clip(edge.coords[x, 
                , drop = FALSE], el[x, , drop = FALSE], params = params, 
                end = "to")
        }))
    }
    x0 <- ec[, 1]
    y0 <- ec[, 2]
    x1 <- ec[, 3]
    y1 <- ec[, 4]
    if (length(loops.e) > 0) {
        ec <- edge.color
        if (length(ec) > 1) {
            ec <- ec[loops.e]
        }
        point.on.cubic.bezier <- function(cp, t) {
            c <- 3 * (cp[2, ] - cp[1, ])
            b <- 3 * (cp[3, ] - cp[2, ]) - c
            a <- cp[4, ] - cp[1, ] - c - b
            t2 <- t * t
            t3 <- t * t * t
            a * t3 + b * t2 + c * t + cp[1, ]
        }
        compute.bezier <- function(cp, points) {
            dt <- seq(0, 1, by = 1/(points - 1))
            sapply(dt, function(t) point.on.cubic.bezier(cp, 
                t))
        }
        plot.bezier <- function(cp, points, color, width, arr, 
            lty, arrow.size, arr.w) {
            p <- compute.bezier(cp, points)
            polygon(p[1, ], p[2, ], border = color, lwd = width, 
                lty = lty)
            if (arr == 1 || arr == 3) {
                igraph:::igraph.Arrows(p[1, ncol(p) - 1], p[2, ncol(p) - 
                  1], p[1, ncol(p)], p[2, ncol(p)], sh.col = color, 
                  h.col = color, size = arrow.size, sh.lwd = width, 
                  h.lwd = width, open = FALSE, code = 2, width = arr.w)
            }
            if (arr == 2 || arr == 3) {
                igraph:::igraph.Arrows(p[1, 2], p[2, 2], p[1, 1], p[2, 
                  1], sh.col = color, h.col = color, size = arrow.size, 
                  sh.lwd = width, h.lwd = width, open = FALSE, 
                  code = 2, width = arr.w)
            }
        }
        loop <- function(x0, y0, cx = x0, cy = y0, color, angle = 0, 
            label = NA, width = 1, arr = 2, lty = 1, arrow.size = arrow.size, 
            arr.w = arr.w, lab.x, lab.y) {
            rad <- angle
            center <- c(cx, cy)
            cp <- matrix(c(x0, y0, x0 + 0.4, y0 + 0.2, x0 + 0.4, 
                y0 - 0.2, x0, y0), ncol = 2, byrow = TRUE)
            phi <- atan2(cp[, 2] - center[2], cp[, 1] - center[1])
            r <- sqrt((cp[, 1] - center[1])^2 + (cp[, 2] - center[2])^2)
            phi <- phi + rad
            cp[, 1] <- cx + r * cos(phi)
            cp[, 2] <- cy + r * sin(phi)
            plot.bezier(cp, 50, color, width, arr = arr, lty = lty, 
                arrow.size = arrow.size, arr.w = arr.w)
            if (is.language(label) || !is.na(label)) {
                lx <- x0 + 0.3
                ly <- y0
                phi <- atan2(ly - center[2], lx - center[1])
                r <- sqrt((lx - center[1])^2 + (ly - center[2])^2)
                phi <- phi + rad
                lx <- cx + r * cos(phi)
                ly <- cy + r * sin(phi)
                if (!is.na(lab.x)) {
                  lx <- lab.x
                }
                if (!is.na(lab.y)) {
                  ly <- lab.y
                }
                text(lx, ly, label, col = edge.label.color, font = edge.label.font, 
                  family = edge.label.family, cex = edge.label.cex)
            }
        }
        ec <- edge.color
        if (length(ec) > 1) {
            ec <- ec[loops.e]
        }
        vs <- vertex.size
        if (length(vertex.size) > 1) {
            vs <- vs[loops.v]
        }
        ew <- edge.width
        if (length(edge.width) > 1) {
            ew <- ew[loops.e]
        }
        la <- loop.angle
        if (length(loop.angle) > 1) {
            la <- la[loops.e]
        }
        lty <- edge.lty
        if (length(edge.lty) > 1) {
            lty <- lty[loops.e]
        }
        arr <- arrow.mode
        if (length(arrow.mode) > 1) {
            arr <- arrow.mode[loops.e]
        }
        asize <- arrow.size
        if (length(arrow.size) > 1) {
            asize <- arrow.size[loops.e]
        }
        awidth <- arrow.width
        if (length(arrow.width) > 1) {
            awidth <- arrow.width[loops.e]
        }
        xx0 <- layout[loops.v, 1] + cos(la) * vs
        yy0 <- layout[loops.v, 2] - sin(la) * vs
        mapply(loop, xx0, yy0, color = ec, angle = -la, label = loop.labels, 
            lty = lty, width = ew, arr = arr, arrow.size = asize, 
            arr.w = awidth, lab.x = loop.labx, lab.y = loop.laby)
    }
    if (length(x0) != 0) {
        if (length(edge.color) > 1) {
            edge.color <- edge.color[nonloops.e]
        }
        if (length(edge.width) > 1) {
            edge.width <- edge.width[nonloops.e]
        }
        if (length(edge.lty) > 1) {
            edge.lty <- edge.lty[nonloops.e]
        }
        if (length(arrow.mode) > 1) {
            arrow.mode <- arrow.mode[nonloops.e]
        }
        if (length(arrow.size) > 1) {
            arrow.size <- arrow.size[nonloops.e]
        }
        if (length(arrow.width) > 1) {
            arrow.width <- arrow.width[nonloops.e]
        }
        if (length(curved) > 1) {
            curved <- curved[nonloops.e]
        }
        if (length(unique(arrow.mode)) == 1) {
            lc <- igraph:::igraph.Arrows(x0, y0, x1, y1, h.col = edge.color, 
                sh.col = edge.color, sh.lwd = edge.width, h.lwd = 1, 
                open = FALSE, code = arrow.mode[1], sh.lty = edge.lty, 
                h.lty = 1, size = arrow.size, width = arrow.width, 
                curved = curved)
            lc.x <- lc$lab.x
            lc.y <- lc$lab.y
        }
        else {
            curved <- rep(curved, length = ecount(graph))[nonloops.e]
            lc.x <- lc.y <- numeric(length(curved))
            for (code in 0:3) {
                valid <- arrow.mode == code
                if (!any(valid)) {
                  next
                }
                ec <- edge.color
                if (length(ec) > 1) {
                  ec <- ec[valid]
                }
                ew <- edge.width
                if (length(ew) > 1) {
                  ew <- ew[valid]
                }
                el <- edge.lty
                if (length(el) > 1) {
                  el <- el[valid]
                }
                lc <- igraph:::igraph.Arrows(x0[valid], y0[valid], x1[valid], 
                  y1[valid], code = code, sh.col = ec, h.col = ec, 
                  sh.lwd = ew, h.lwd = 1, h.lty = 1, sh.lty = el, 
                  open = FALSE, size = arrow.size, width = arrow.width, 
                  curved = curved[valid])
                lc.x[valid] <- lc$lab.x
                lc.y[valid] <- lc$lab.y
            }
        }
        if (!is.null(elab.x)) {
            lc.x <- ifelse(is.na(elab.x), lc.x, elab.x)
        }
        if (!is.null(elab.y)) {
            lc.y <- ifelse(is.na(elab.y), lc.y, elab.y)
        }
        text(lc.x, lc.y, labels = edge.labels, col = edge.label.color, 
            family = edge.label.family, font = edge.label.font, 
            cex = edge.label.cex)
    }
    rm(x0, y0, x1, y1)
    if (length(unique(shape)) == 1) {
        igraph:::.igraph.shapes[[shape[1]]]$plot(layout, params = params)
    }
    else {
        sapply(seq(length = vcount(graph)), function(x) {
            igraph:::.igraph.shapes[[shape[x]]]$plot(layout[x, , drop = FALSE], 
                v = x, params = params)
        })
    }
    par(xpd = TRUE)
    x <- layout[, 1] + label.dist * cos(-label.degree) * (vertex.size + 
        6 * 8 * log10(2))/200
    y <- layout[, 2] + label.dist * sin(-label.degree) * (vertex.size + 
        6 * 8 * log10(2))/200
    if (length(label.family) == 1) {
        text(x, y, labels = labels, col = label.color, family = label.family, 
            font = label.font, cex = label.cex)
    }
    else {
        if1 <- function(vect, idx) if (length(vect) == 1) 
            vect
        else vect[idx]
        sapply(seq_len(vcount(graph)), function(v) {
            text(x[v], y[v], labels = if1(labels, v), col = if1(label.color, 
                v), family = if1(label.family, v), font = if1(label.font, 
                v), cex = if1(label.cex, v))
        })
    }
    rm(x, y)
    invisible(NULL)
}
