
plotPhylogram.mjl = function (tree, colors, fsize, ftype, lwd, pts, node.numbers, 
    mar, add, offset, direction, setEnv, xlim, ylim, placement, 
    tips, split.vertical, lend, asp, x_offset) 
{
    if (split.vertical && !setEnv) {
        cat("split.vertical requires setEnv=TRUE. Setting split.vertical to FALSE.\n")
        spit.vertical <- FALSE
    }
    offsetFudge <- 1.37
    cw <- reorderSimmap(tree)
    pw <- reorderSimmap(tree, "postorder")
    n <- Ntip(cw)
    m <- cw$Nnode
    Y <- matrix(NA, m + n, 1)
    if (is.null(tips)) 
        Y[cw$edge[cw$edge[, 2] <= n, 2]] <- 1:n
    else Y[cw$edge[cw$edge[, 2] <= n, 2]] <- if (is.null(names(tips))) 
        tips[sapply(1:Ntip(cw), function(x, y) which(y == x), 
            y = cw$edge[cw$edge[, 2] <= n, 2])]
    else tips[gsub(" ", "_", cw$tip.label)]
    nodes <- unique(pw$edge[, 1])
    for (i in 1:m) {
        if (placement == "intermediate") {
            desc <- pw$edge[which(pw$edge[, 1] == nodes[i]), 
                2]
            Y[nodes[i]] <- (min(Y[desc]) + max(Y[desc]))/2
        }
        else if (placement == "centered") {
            desc <- getDescendants(pw, nodes[i])
            desc <- desc[desc <= Ntip(pw)]
            Y[nodes[i]] <- (min(Y[desc]) + max(Y[desc]))/2
        }
        else if (placement == "weighted") {
            desc <- pw$edge[which(pw$edge[, 1] == nodes[i]), 
                2]
            n1 <- desc[which(Y[desc] == min(Y[desc]))]
            n2 <- desc[which(Y[desc] == max(Y[desc]))]
            v1 <- pw$edge.length[which(pw$edge[, 2] == n1)]
            v2 <- pw$edge.length[which(pw$edge[, 2] == n2)]
            Y[nodes[i]] <- ((1/v1) * Y[n1] + (1/v2) * Y[n2])/(1/v1 + 
                1/v2)
        }
        else if (placement == "inner") {
            desc <- getDescendants(pw, nodes[i])
            desc <- desc[desc <= Ntip(pw)]
            mm <- which(abs(Y[desc] - median(Y[1:Ntip(pw)])) == 
                min(abs(Y[desc] - median(Y[1:Ntip(pw)]))))
            if (length(mm > 1)) 
                mm <- mm[which(Y[desc][mm] == min(Y[desc][mm]))]
            Y[nodes[i]] <- Y[desc][mm]
        }
    }
    H <- nodeHeights(cw)

    H = H + x_offset
    par(mar = mar)
    if (is.null(offset)) 
        offset <- 0.2 * lwd/3 + 0.2/3
    if (!add) 
        plot.new()
    if (is.null(xlim)) {
        pp <- par("pin")[1]
        sw <- fsize * (max(strwidth(cw$tip.label, units = "inches"))) + 
            offsetFudge * fsize * strwidth("W", units = "inches")
        alp <- optimize(function(a, H, sw, pp) (a * 1.04 * max(H) + 
            sw - pp)^2, H = H, sw = sw, pp = pp, interval = c(0, 
            1e+06))$minimum
        xlim <- if (direction == "leftwards") 
            c(min(H) - sw/alp, max(H))
        else c(min(H), max(H) + sw/alp)
    }
    if (is.null(ylim)) 
        ylim = range(Y)
    if (direction == "leftwards") 
        H <- max(H) - H
    plot.window(xlim = xlim, ylim = ylim, asp = asp)
    if (!split.vertical) {
        for (i in 1:m) lines(H[which(cw$edge[, 1] == nodes[i]), 
            1], Y[cw$edge[which(cw$edge[, 1] == nodes[i]), 2]], 
            col = colors[names(cw$maps[[match(nodes[i], cw$edge[, 
                1])]])[1]], lwd = lwd)
    }
    for (i in 1:nrow(cw$edge)) {
        x <- H[i, 1]
        for (j in 1:length(cw$maps[[i]])) {
            if (direction == "leftwards") 
                lines(c(x, x - cw$maps[[i]][j]), c(Y[cw$edge[i, 
                  2]], Y[cw$edge[i, 2]]), col = colors[names(cw$maps[[i]])[j]], 
                  lwd = lwd, lend = lend)
            else lines(c(x, x + cw$maps[[i]][j]), c(Y[cw$edge[i, 
                2]], Y[cw$edge[i, 2]]), col = colors[names(cw$maps[[i]])[j]], 
                lwd = lwd, lend = lend)
            if (pts) 
                points(c(x, x + cw$maps[[i]][j]), c(Y[cw$edge[i, 
                  2]], Y[cw$edge[i, 2]]), pch = 20, lwd = (lwd - 
                  1))
            x <- x + if (direction == "leftwards") 
                -cw$maps[[i]][j]
            else cw$maps[[i]][j]
            j <- j + 1
        }
    }
    if (node.numbers) {
        symbols(if (direction == "leftwards") 
            max(H)
        else 0, mean(Y[cw$edge[cw$edge[, 1] == (Ntip(cw) + 1), 
            2]]), rectangles = matrix(c(1.2 * fsize * strwidth(as.character(Ntip(cw) + 
            1)), 1.4 * fsize * strheight(as.character(Ntip(cw) + 
            1))), 1, 2), inches = FALSE, bg = "white", add = TRUE)
        text(if (direction == "leftwards") 
            max(H)
        else 0, mean(Y[cw$edge[cw$edge[, 1] == (Ntip(cw) + 1), 
            2]]), Ntip(cw) + 1, cex = fsize)
        for (i in 1:nrow(cw$edge)) {
            x <- H[i, 2]
            if (cw$edge[i, 2] > Ntip(cw)) {
                symbols(x, Y[cw$edge[i, 2]], rectangles = matrix(c(1.2 * 
                  fsize * strwidth(as.character(cw$edge[i, 2])), 
                  1.4 * fsize * strheight(as.character(cw$edge[i, 
                    2]))), 1, 2), inches = FALSE, bg = "white", 
                  add = TRUE)
                text(x, Y[cw$edge[i, 2]], cw$edge[i, 2], cex = fsize)
            }
        }
    }
    if (direction == "leftwards") 
        pos <- if (par()$usr[1] > par()$usr[2]) 
            4
        else 2
    if (direction == "rightwards") 
        pos <- if (par()$usr[1] > par()$usr[2]) 
            2
        else 4
    for (i in 1:n) if (ftype) 
        text(H[which(cw$edge[, 2] == i), 2], Y[i], cw$tip.label[i], 
            pos = pos, offset = offset, cex = fsize, font = ftype)
    if (setEnv) {
        PP <- list(type = "phylogram", use.edge.length = TRUE, 
            node.pos = 1, show.tip.label = if (ftype) TRUE else FALSE, 
            show.node.label = FALSE, font = ftype, cex = fsize, 
            adj = 0, srt = 0, no.margin = FALSE, label.offset = offset, 
            x.lim = xlim, y.lim = ylim, direction = direction, 
            tip.color = "black", Ntip = Ntip(cw), Nnode = cw$Nnode, 
            edge = cw$edge, xx = sapply(1:(Ntip(cw) + cw$Nnode), 
                function(x, y, z) y[match(x, z)], y = H, z = cw$edge), 
            yy = Y[, 1])
        assign("last_plot.phylo", PP, envir = .PlotPhyloEnv)
    }
    if (split.vertical) 
        splitEdgeColor(cw, colors, lwd)
}

plotSimmap.mjl = function (tree, colors = NULL, fsize = 1, ftype = "reg", lwd = 2, 
    pts = FALSE, node.numbers = FALSE, mar = NULL, add = FALSE, 
    offset = NULL, direction = "rightwards", type = "phylogram", 
    setEnv = TRUE, part = 1, xlim = NULL, ylim = NULL, nodes = "intermediate", 
    tips = NULL, maxY = NULL, hold = TRUE, split.vertical = FALSE, 
    lend = 2, asp = NA, x_offset=0) 
{
    if (inherits(tree, "multiPhylo")) {
        par(ask = TRUE)
        for (i in 1:length(tree)) plotSimmap(tree[[i]], colors = colors, 
            fsize = fsize, ftype = ftype, lwd = lwd, pts = pts, 
            node.numbers = node.numbers, mar, add, offset, direction, 
            type, setEnv, part, xlim, ylim, nodes, tips, maxY, 
            hold, split.vertical, lend)
    }
    else {
        if (!inherits(tree, "phylo")) 
            stop("tree should be object of class \"phylo\"")
        if (is.null(tree$maps)) 
            stop("tree should contain mapped states on edges.")
        ftype <- which(c("off", "reg", "b", "i", "bi") == ftype) - 
            1
        if (!ftype) 
            fsize = 0
        if (is.null(colors)) {
            st <- sort(unique(unlist(sapply(tree$maps, names))))
            colors <- palette()[1:length(st)]
            names(colors) <- st
            if (length(st) > 1) {
                cat("no colors provided. using the following legend:\n")
                print(colors)
            }
        }
        tree$tip.label <- gsub("_", " ", tree$tip.label)
        if (is.null(mar)) 
            mar = rep(0.1, 4)
        if (hold) 
            null <- dev.hold()
        if (type == "phylogram") {
            if (direction %in% c("upwards", "downwards")) 
                updownPhylogram(tree, colors, fsize, ftype, lwd, 
                  pts, node.numbers, mar, add, offset, direction, 
                  setEnv, xlim, ylim, nodes, tips, split.vertical, 
                  lend, asp)
            else plotPhylogram.mjl(tree, colors, fsize, ftype, lwd, 
                pts, node.numbers, mar, add, offset, direction, 
                setEnv, xlim, ylim, nodes, tips, split.vertical, 
                lend, asp, x_offset)
        }
        else if (type == "fan") {
            plotFan(tree, colors, fsize, ftype, lwd, mar, add, 
                part, setEnv, xlim, ylim, tips, maxY, lend)
        }
        if (hold) 
            null <- dev.flush()
    }
}

