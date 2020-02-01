# R dependencies
library(ape)
library(cowplot)
library(data.table)
library(gginnards)
library(ggplot2)
library(ggtree)
library(grid)
library(igraph)
library(RevGadgets)



### GGPLOT UTILS

# used to manipulate plot legend
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# used to rearrange plot layers
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



### PHYLO FORMATTING

# generates information about clade assignments
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

# assigns colors to sequenced (black) and unsequenced (gray) taxon labels
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

# corrects various taxon labels
fix_vib_tip = function(phy) {
    lbl = phy$tip.label
    lbl = as.vector(sapply(lbl, function(x) { strsplit(x, "_")[[1]][2] } ))
    lbl[ lbl=="rigidum" ] = "rugosum"
    lbl[ lbl=="stellato" ] = "stellato-tomentosum"
    phy$tip.label = lbl
    return(phy)
}

# formats nodes to plot HPDs (modified from code by Sebastian Duchene)
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

# creates HPDs for ggtree plot
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

# adds epoch boxes to plots
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



### AREA + BIOME FORMATTING

# makes list of colors with state labels
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

# returns the (area,biome) to (area,biome) transition event (base-1)
make_compound_tx = function( from_bg_idx, from_biome_idx, to_bg_idx, to_biome_idx ) {
    ret = matrix(NA, nrow=0, ncol=2)
    
    is_biome_event = !identical(from_biome_idx, to_biome_idx)
    is_bg_event = !identical(from_bg_idx, to_bg_idx)
    if (is_biome_event && is_bg_event) {
        cat("from_bg_idx    =", from_bg_idx, "\n")
        cat("to_bg_idx      =", to_bg_idx, "\n")
        cat("from_biome_idx =", from_biome_idx, "\n")
        cat("to_biome_idx   =", to_biome_idx, "\n")
        error("Event identified as both a biome and bg event!")
    }
    is_dispersal = FALSE
    if (is_bg_event && (length(from_bg_idx) < length(to_bg_idx))) {
        is_dispersal = TRUE
    }
    
    for (from_bg in from_bg_idx) {
        for (to_bg in to_bg_idx) {
            for (from_biome in from_biome_idx) {
                for (to_biome in to_biome_idx) {
                    if (is_biome_event && from_bg == to_bg) {
                        # only identify biome shifts within each region
                        from_state = get_compound_state(from_bg, from_biome)
                        to_state = get_compound_state(to_bg, to_biome)
                        ret = rbind(ret, c(from_state, to_state))
                    } else if (is_dispersal && from_bg != to_bg) {
                        from_state = get_compound_state(from_bg, from_biome)
                        to_state = get_compound_state(to_bg, to_biome)
                        ret = rbind(ret, c(from_state, to_state))
                    }
                }
            }
        }
    }
    colnames(ret) = c("from","to")
    return(ret)
}

# returns the (area,biome) state pair as a vector (base-1)
get_area_biome = function(s, n_areas=6, n_biomes=4) {
    biome = as.integer( (s-1) / n_areas)
    area = s - biome*n_areas
    return( c(area, biome+1) )
}

# returns the range-vector as integers (base-1)
get_bg_state = function(s) {
    if (s==1)       return(c(1))
    else if (s==2)  return(c(2))
    else if (s==3)  return(c(3))
    else if (s==4)  return(c(4))
    else if (s==5)  return(c(5))
    else if (s==6)  return(c(6))
    else if (s==7)  return(c(1,2))
    else if (s==8)  return(c(1,3))
    else if (s==9)  return(c(2,3))
    else if (s==10) return(c(1,4))
    else if (s==11) return(c(2,4))
    else if (s==12) return(c(3,4))
    else if (s==13) return(c(1,5))
    else if (s==14) return(c(2,5))
    else if (s==15) return(c(3,5))
    else if (s==16) return(c(4,5))
    else if (s==17) return(c(1,6))
    else if (s==18) return(c(2,6))
    else if (s==19) return(c(3,6))
    else if (s==20) return(c(4,6))
    else if (s==21) return(c(5,6))
}

# returns the biome state as an integer (base-1)
get_biome_state = function(s) {
    return(s+1)
}

# returns the biome-area state as an integer (base-1)
get_compound_state = function(area, biome, n_areas=6, n_biomes=4) {
    return( (n_areas)*(biome-1) + area )
}



### SENSITIVITY ANALYSIS

# probability of freezing origin across replicates
pi_temp_prob = function(dat) {
    x = dat[,ncol(dat)]
    p = sum(x==3) / length(x)
    return(p)
}

# computes entropy between prior and estimate
compute_entropy = function(dat) {
    x = dat[,ncol(dat)]
    p = c()
    for (j in 1:4) {
        p[j] = sum(x==j)
    }
    p = p / length(x)
    p = p[ p != 0 ]
    
    return( -sum( p * log(p, 2) ) )
}

# computes entropy under (expected) uniform dist
compute_entropy_exp = function(x) {
    z = 1-x
    p = c(z/3, z/3, z/3, x)
    p = p[ p != 0 ]
    
    return( -sum( p * log(p, 2) ) )
}



### STATES THROUGH TIME

# Used in Figs 7,8
# combines biome and biogeographic states across two ancestral state .tsv files
make_joint_state = function(b1, b2) {
    
    # create joint events from independent events
    b1_tmp = cbind(b1, start_state_2=NA, end_state_2=NA, transition_char=1)
    names( b1_tmp )[ match(c("start_state", "end_state"), names(b1_tmp) ) ] = c("start_state_1","end_state_1")
    b2_tmp = cbind(b2, start_state_1=NA, end_state_1=NA, transition_char=2)
    names( b2_tmp )[ match(c("start_state", "end_state"), names(b2_tmp) ) ] = c("start_state_2","end_state_2")
    b = rbind( as.data.frame(b1_tmp, stringsAsFactors=F),
               as.data.frame(b2_tmp, stringsAsFactors=F))

    # sort joint events by transition time
    b = b[order(b$transition_time, decreasing=T),]

    # get the first valid start state for each character
    start_state_1 = na.omit(b$start_state_1)[1]
    start_state_2 = na.omit(b$start_state_2)[1]
    
    # there must be at least two events; drop first event from list if it's no_change
    if (b$transition_type[1]=="no_change") {
        b = b[-c(1),]    
    }
    
    # update missing states for the first event
    if (b$transition_char[1]==1) {
        b[1, c("start_state_2", "end_state_2")] = c( start_state_2, start_state_2 )            
    } else if (b$transition_char[1]==2) {
        b[1, c("start_state_1", "end_state_1")] = c( start_state_1, start_state_1 )
    }
    
    # if the next remaining event is the last event, return
    if (nrow(b) > 1) {
        
        # if there's more than one event, propagate joint states
        for (i in 2:nrow(b)) {
            
            # update end states for next iteration
            end_state_1 = b$end_state_1[i-1]
            end_state_2 = b$end_state_2[i-1]
            
            # does the transition occur on character 1 or 2?
            char_idx = b$transition_char[i]
            if (char_idx==1) {
                b[i,c("start_state_2","end_state_2")] = c( end_state_2, end_state_2)            
            } else if (char_idx==2) {
                b[i,c("start_state_1","end_state_1")] = c( end_state_1, end_state_1)
            }
        }
    }
    if (nrow(b)==1) {
        x1 = b$branch_start_time
        x2 = b$branch_end_time
    } else {
        x1 = c( b$branch_start_time[1], b$transition_time )
        x2 = c( b$transition_time, b$branch_end_time[1] )
        b = rbind(as.data.frame(b, stringsAsFactors=F), 
                  as.data.frame(b[nrow(b),], stringsAsFactors=F))
        b[nrow(b),c("start_state_1","start_state_2","transition_time")] = b[nrow(b),c("end_state_1","end_state_2","branch_end_time")]
    }
    b = cbind(b, x1=x1, x2=x2, s1=b$start_state_1, s2=b$start_state_2)
    return(b)
}

# Prepares Stage 1 dataframe for Stage 2 in Figs 7,8
make_stoch_bg_biome = function(stoch_bg, stoch_biome, iterations) {
    branches = 1:max(unique(stoch_bg$parent_index), na.rm=T)
    # loop over iterations
    stoch_list = list()
    for (i in 1:length(iterations)) {
        
        # get biome and biogeography stochastic mappings per iteration
        it = iterations[i]
        cat("Stage 1, processing iteration ",it," / ", max(iterations), "\n", sep="")
        sample_bg = stoch_bg[ stoch_bg$iteration==it, ]
        sample_biome = stoch_biome[ stoch_biome$iteration==it, ]
        
        # loop over branches
        tmp_branch_list = list()
        for (j in 1:length(branches)) {
            
            # get biome and biogeography stochastic mappings per branch
            nd_idx = branches[j]
            branch_bg = sample_bg[ sample_bg$node_index==nd_idx, ]
            branch_biome = sample_biome[ sample_biome$node_index==nd_idx, ]
            
            # interleave biome and biogeography stochastic mappings
            tmp_branch_list[[j]] = as.data.frame( make_joint_state(branch_bg, branch_biome), stringsAsFactors=F )
        }
        stoch_list[[i]] = rbindlist(tmp_branch_list) 
    }
    stoch_bg_biome = rbindlist(stoch_list)
    return(stoch_bg_biome)
}

# Prepares Stage 2 dataframe for Stage 3 in Fig 7
make_bg_biome_bins = function(stoch_bg_biome, bin_width, iterations, n_areas, n_biomes) {
        
    # bins
    # index ( bg x biome x time )
    bin_width = 1
    max_time = 90
    n_bins = max_time / bin_width
    state_bins = array(0, dim=c(n_areas, n_biomes, n_bins))
    
    # 0.0 to 0.5, 0.5 to 1.0, etc
    ages = seq(0.0, max_time, by=bin_width)
    
    dat_plot_colnames = c( names(stoch_bg_biome[c(),]), "age", "joint_state" )
    dat_plot = data.frame(matrix(ncol=length(dat_plot_colnames), nrow=0), stringsAsFactors=F)
    colnames(dat_plot) = dat_plot_colnames
    
    dat_tmp = data.frame(matrix(ncol=length(dat_plot_colnames), nrow=1e3), stringsAsFactors=F)
    colnames(dat_tmp) = dat_plot_colnames
    
    idx_tmp = 1
    curr_it = -1
    for (i in 1:nrow(stoch_bg_biome)) {
        
    #for (i in 1:2500) {
        if (curr_it != stoch_bg_biome$iteration[i]) {
            curr_it = stoch_bg_biome$iteration[i]
            cat("Stage 2, processing iteration ",curr_it," / ", max(stoch_bg_biome$iteration), "\n", sep="")
        }
        
        bg_idx = get_bg_state( stoch_bg_biome$start_state_1[i] )
        biome_idx = get_biome_state( stoch_bg_biome$start_state_2[i] )
        
        #age_bins = seq(floor(stoch_bg_biome$x2[i]), ceiling(stoch_bg_biome$x1[i]), by=bin_width ) / bin_width
        start_age = floor(stoch_bg_biome$x2[i])
        end_age = ceiling(stoch_bg_biome$x1[i])
        age_bins = start_age:end_age * bin_width
        time_idx = age_bins + 1
        
        #print(i)
        for (j in 1:length(age_bins)) {
            for (k in 1:length(bg_idx)) {
                joint_state = paste(bg_idx[k],biome_idx,sep="_")
                #print(joint_state)
                #dat_plot = rbind(dat_plot, c( stoch_bg[i,], age=age_bins[j], joint_state=joint_state))
                dat_tmp[idx_tmp,] = c( stoch_bg[i,], age=age_bins[j], joint_state=joint_state)
                
                if (idx_tmp == nrow(dat_tmp)) {
                    dat_plot = rbind(dat_plot, dat_tmp)
                    idx_tmp = 1
                } else if (idx_tmp < nrow(stoch_bg_biome)) {
                    idx_tmp = idx_tmp + 1
                }
            }
        }
        #print(time_idx)
        #time_idx = n_bins - seq(floor(stoch_bg_biome$x2[i]), ceiling(stoch_bg_biome$x1[i]), by=bin_width ) / bin_width + 1
        state_bins[ bg_idx, biome_idx, time_idx ] = state_bins[ bg_idx, biome_idx, time_idx ] + 1
    }
    return(state_bins)
}

# Prepares Stage 3 dataframe for plotting in Fig 7
make_lstt_dat = function(state_bins, ages, min_sample=510) {
    
    ret = list()
    
    # create a melted data frame with Count/Support for Area/Biome over time (Age)
    d1 = matrix(nrow=0, ncol=dim(state_bins)[1])
    colnames(d1) = c("age","count","Area","Biome","Area_Biome","Support")
    for (i in 1:dim(state_bins)[1]) {
        for (j in 1:dim(state_bins)[2]) {
            for (k in 1:dim(state_bins)[3]) {
                d1 = rbind(d1, c( ages[k], state_bins[i,j,k], i, j, paste(i, j, sep="_"), 0))
            }
        }
    }
    
    # prepare column values
    d2         = data.frame(d1, stringsAsFactors=FALSE)
    d2$age     = as.numeric(d2$age)
    d2$count   = as.numeric(d2$count)
    d2$Support = as.numeric(d2$Support)
    
    # compute confidence in state for each time step using
    # multinomial confidence metric (SK Ernst)
    biome_conf = t(apply( state_bins, 3, colSums))
    bg_conf    = t(apply( state_bins, 3, rowSums))
    min_sample = 510
    for (i in 1:n_bins) {
        for (j in 1:n_biomes) {
            if (biome_conf[i,j] > min_sample) { 
                biome_conf[i,j] = 1
            } else {
                biome_conf[i,j] = 0
            }
        }
        for (j in 1:n_areas) {
            if (bg_conf[i,j] > min_sample) {
                bg_conf[i,j] = 1
            } else {
                bg_conf[i,j] = 0
            }
        }
    }
    
    # only show time-bins that have contain more samples than min_sample
    d2_ages = unique(d2$age)
    d2_bg_trunc = d2
    d2_biome_trunc = d2
    
    for (i in 1:length(d2_ages)) {
        for (j in 1:n_areas) {
            c_ij = d2_bg_trunc[ d2_bg_trunc$age==i & d2_bg_trunc$Area==j, ]$count
            if (length(c_ij) == 0) { 
                # do nothing
            } else {
                d2_bg_trunc[ d2_bg_trunc$age==i & d2_bg_trunc$Area==j, ]$Support = bg_conf[i,j]
            }
        }
        
        for (j in 1:n_biomes) {
            c_ij = d2_biome_trunc[ d2_biome_trunc$age==i & d2_biome_trunc$Biome==j, ]$count
            if (length(c_ij) == 0) { 
                # do nothing
            } else {
                d2_biome_trunc[ d2_biome_trunc$age==i & d2_biome_trunc$Biome==j, ]$Support = biome_conf[i,j]
            }
        }
    }
    ret$bg = d2_bg_trunc
    ret$biome = d2_biome_trunc
    return(ret)
}

# Computes various statistics for minimum number of samples needed to accurately
# approximate multinomial frequences, per SK Ernest 1978
calculate_ske = function(s, k, alpha=0.05, D=0.05) {
    # What is the value of d for the number of samples?
    # SK Ernest 1978 -- compute smallest sample size, n, that satisfies inequality
    #   \prod_i Pr( |f_i - \pi_i| < d ) > 1 - \alpha
    # where d determines error in obs f_i from true parameter \pi_i
    #
    # We can also solve for alpha given D and N, and D given N and alpha.
    
    # We compute the posterior mean for each bin/region/biome thing, S_ijt,
    # where i is region, j is biome, and t is the time bin.
    # We then compute S_it = sum_j S_ijt and S_jt = sum_i S_ijt.
    # Some time bins have very few samples. We measure support in two ways using
    # a frequentist significance test for the multinomial distribution.
    #
    # 1) SK Ernest computes the minimum number of samples, $N'$, needed to 
    #    guarantee that all 95% of replicates simulated under the true values
    #    of the approximated multinomial distribution would produce sampled
    #    frequencies within +/- 0.05 units of the true simulating parameters.
    #    where N' = 510 regardless of the number of categories.
    #
    # 2) We use this same theory to solve for what is the alpha significance
    #    for the number of actual samples we have -- i.e. how often
    #    do we expect to generate sampled proportions +/- d of the true values??
    #
    
    # We report values of alpha for each time bin. Alpha values correspond to
    # the probability that we would expect proportions with values, f_i, that
    # similar to multinomial probabilities, p_i, such that |f_i-p_i| <= 0.05
    # provided the sample size is equals what we sampled for each LSTT time bin
    # by MCMC.
    
    # We report our support metric as the proportion of samples N_it = 
    # and the proportion of th
    # 
    # 
    # The support statistic is defined for the $S_{i,t} = S_{i,j,t}$ samples
    # drawn for time bin $t$ for biome/region $i$. The statistic
    # reports the probability, $\alpha$, that all category frequencies are
    # within the true underlying posterior proportions $\pi_{i,j,t}$
    # 

    # s     : number of samples
    # k     : number of bins
    # alpha : significance level
    # D     : estimate difference must satisfy | f_i - pi_i | > D for all i in k
    # N     : min num samples needed for
    n = c()
    d = c()
    a = c()

    if (s<=1) {
        return( list(n=Inf, d=1, a=1, alpha=alpha, s=s, D=D) )
    }
    for (m in 1:k) {
        # significance value for alpha/m
        z = qnorm(p=(1-alpha/(2*m)),mean=0,sd=1)
        # how many samples, n, needed for z-values tolerating error D
        n[m] = ceiling(z^2*(1/m) * (1 - 1/m) / D^2)
        # how accurately can you estimate pi_i for sample size s w/ z-values
        d[m] = sqrt(z^2*(1/m) * (1 - 1/m) / s )
        # what alpha level can you obtain with s samples for tolerance D
        aval = seq(0,1,by=0.001)
        aval = aval[2:(length(aval)-1)]
        # get z-values for possible alpha values
        zz = qnorm(p=(1-aval/(2*m)),mean=0,sd=1)
        # compute n values for possible alpha values under D
        nn = (zz^2*(1/m) * (1 - 1/m) / D^2)
        # for values of n greater than the actual sample, which is the smallest?
        a_idx = min( sum(s < nn) + 1, length(nn) )
        # choose that alpha value
        a[m] = aval[a_idx]
    }
    return( list(n=max(n), d=max(d), a=max(a, na.rm=T), alpha=alpha, s=s, D=D) )
}


### TRANSITION GRAPH

# check binned states to ensure no 'double-events' occurred, in which
# both the biome and area change simultaneously
check_bin_diag_zero = function(x) {
    n_states = dim(x)[1]
    n_bins = dim(x)[3]
    for (i in 1:n_states) {
        for (j in 1:n_bins) {
            if (x[i,i,j] != 0) {
                error("bin",j,"+","state",i,"=",x[i,i,j],"> 0\n")
            }
        }
    }   
}

# gets radian for xy coords w/r/t (0,0)
get_radian = function( xy, origin=c(0,0) ) {
    dx = xy[1] - origin[1]
    dy = xy[2] - origin[2]
    rad = atan2(dy, dx)
    #xy[1]-origin[0]  
    return(rad)
}

# rotates coordinate system by angle
coords_rotate = function(x, angle=0, origin=c(0,0)) {
    for (i in 1:nrow(x)) {
        xy = x[i,]
        xy_new = xy
        xy_new[1] = origin[1] + cos(angle) * (xy[1] - origin[1]) - sin(angle) * (xy[2] - origin[2])
        xy_new[2] = origin[2] + sin(angle) * (xy[1] - origin[1]) + cos(angle) * (xy[2] - origin[2])
        x[i,] = xy_new
    }
    return(x)
}

# label elements in x by name if they have values >= p
make_count_class = function(x, p) {
    n = length(p)
    y = rep(0, length(x))
    for (i in 1:length(x)) {
        z = sum(x[i] >= p)
        if (z == 0) {
            y[i] = NA #paste("<", p[1], sep="")
        } else {
            y[i] = paste(">", p[z], sep="")
        }
    }
   return(y)
}

# create custom circle shape
mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}

# register new 'mycircle' shape with igraph
add.vertex.shape("fcircle",
                 clip=shapes("circle")$clip,
                 plot=mycircle,
                 parameters=list(vertex.frame.color=1, vertex.frame.width=1))


# Prepares Stage 2 dataframe for plotting in Fig 8
make_bg_biome_tx = function(stoch_bg_biome, bin_width, iterations, state_lbl) {
    
    n_states = length(state_lbl)
    n_it = length(iterations)
    bin_width  = 1
    state_bins = array(0, dim=c(n_states, n_states, n_iter))
    hit_bins   = array(0, dim=c(n_states, n_states, n_iter))
    #ages       = seq(0.0, max_time, by=bin_width)
    
    # populate the LSTT matrix
    curr_it = -1
    for (i in 1:nrow(stoch_bg_biome)) {
        
        if (curr_it != stoch_bg_biome$iteration[i]) {
            curr_it = stoch_bg_biome$iteration[i]
            cat("Stage 2, processing iteration ",curr_it," / ", max(stoch_bg_biome$iteration), "\n", sep="")
        }
        it_idx = match(curr_it, iterations)
        
        if (stoch_bg_biome$transition_type[i] == "anagenetic") {
    
            from_bg_idx      = get_bg_state( stoch_bg_biome$start_state_1[i] )
            from_biome_idx   = get_biome_state( stoch_bg_biome$start_state_2[i] )
            to_bg_idx        = get_bg_state( stoch_bg_biome$end_state_1[i] )
            to_biome_idx     = get_biome_state( stoch_bg_biome$end_state_2[i] )
    
            # get the matrix of state transitions
            state_tx = make_compound_tx( from_bg_idx, from_biome_idx, to_bg_idx, to_biome_idx )
            
            # store each event's info as count or hit
            if (nrow(state_tx) > 0) {
                for (i in 1:nrow(state_tx)) {
                    state_bins[ state_tx[i,1], state_tx[i,2], it_idx ] = state_bins[ state_tx[i,1], state_tx[i,2], it_idx ] + 1
                    hit_bins[ state_tx[i,1], state_tx[i,2], it_idx ] = 1
                }
            }
        }
    }
    
    # validate bins don't include self transitions
    check_bin_diag_zero(state_bins)
    
    # convert from posterior samples
    mean_event = rowSums(state_bins,dims=2) * (1/n_iter)
    prob_event = rowSums(hit_bins,dims=2) * (1/n_iter)
    
    # format & filter
    rownames(mean_event)=state_lbl; colnames(mean_event)=state_lbl
    rownames(prob_event)=state_lbl; colnames(prob_event)=state_lbl
    
    return(mean_event)
}



# STOCHASTIC MAPPING

# Plots phylogeny for stochastic mapping for Figs 3,4,5
plotPhylogram.vib_div = function (tree, colors, fsize, ftype, lwd, pts, node.numbers, 
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

# Plots phylogeny with stochastic mapping for Figs 3,4,5
plotSimmap.vib_div = function (tree, colors = NULL, fsize = 1, ftype = "reg", lwd = 2, 
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
            else plotPhylogram.vib_div(tree, colors, fsize, ftype, lwd, 
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

