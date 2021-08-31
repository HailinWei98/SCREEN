#' @export
#' @import Gviz
#' @import ensembldb
ciceroPlot<- function(score_dir, pval_dir, species = "Hs", version = "v75", annotations = NULL, 
                      p_val_cut = 0.05, score_cut = 0, selected, upstream = 2000000, downstream = 2000000, 
                      track_size = c(1,.3,.2,.3), include_axis_track = TRUE, connection_color = "#7F7CAF",
                      connection_color_legend = TRUE, connection_width = 2, connection_ymax = NULL, 
                      gene_model_color = "#81D2C7", alpha_by_coaccess = FALSE, 
                      gene_model_shape = c("smallArrow", "box")){
    
    color_names = NULL
    
    # get genome annotations
    
    if(is.null(annotations)){
        if(species == "Hs"){
            if(version == "v75"){
                a <- genes(EnsDb.Hsapiens.v75)
            }else if(version == "v79"){
                a <- genes(EnsDb.Hsapiens.v79)
            }else if(version == "v86"){
                a <- genes(EnsDb.Hsapiens.v86)
            }  
        }else if(species == "Mm"){
            if(version == "v75"){
                a <- genes(EnsDb.Mmusculus.v75)
            }else if(version == "v79"){
                a <- genes(EnsDb.Mmusculus.v79)
            }
        }
        gene_anno<- data.frame(a)
        gene_anno$chromosome <- paste0("chr", gene_anno$seqnames)
        gene_anno$transcript <- gene_anno$symbol
    }else{
        gene_anno <- annotations
            }
        
    #get score and p-value
    
    if(is.character(score_dir)){
        score <- read.table(score_dir, header = T, row.names = 1)
    }else{
        score <- score_dir
    }
    if(is.character(pval_dir)){
        pval <- read.table(pval_dir, header = T, row.names = 1)
    }else{
        pval <- pval_dir
    }
    
    #get information from selected enhancer
    
    chr<- unlist(strsplit(selected, "[.|:|-]"))[1]
    start<- as.numeric(unlist(strsplit(selected, "[.|:|-]"))[2])
    end<- as.numeric(unlist(strsplit(selected, "[.|:|-]"))[3])
        
    #extend region
    
    minbp<- start - upstream
    maxbp<- end + downstream
    
    #get gene model

    gene_model <- gene_anno
    gene_model <- gene_model[!is.na(gene_model$chromosome) & 
                             !is.na(gene_model$start) &
                             !is.na(gene_model$end) &
                             !is.na(gene_model$strand) &
                             !is.na(gene_model$transcript), ]
    gene_model <- gene_model[gene_model$chromosome == chr &
                   ((gene_model$start > minbp & gene_model$start < maxbp) |
                    (gene_model$end > minbp & gene_model$end < maxbp) |
                    (gene_model$start < minbp & gene_model$end > maxbp)), ]
    gene_model <- gene_model[gene_model$transcript %in% rownames(score), ]
    
    #get score and p-value of selected enhancer
    
    conn_input <- data.frame(score = score[,selected], pval = pval[, selected])
    rownames(conn_input) <- rownames(score)
    conn_input <- conn_input[rownames(conn_input) %in% gene_model$transcript, ]
    conn_input <- subset(conn_input, ((score > score_cut) | (score < -score_cut)) & pval < p_val_cut)
    gene_model <- subset(gene_model, transcript %in% rownames(conn_input))
    grtrack <- Gviz::GeneRegionTrack(gene_model, chromosome = chr, geneSymbols = TRUE,
                                     name = "", fill = "#81D2C7",
                                     col = "#81D2C7",  fontcolor = "black",
                                     fontcolor.group = "black", fontsize.group = 6,
                                     fontsize = 6, shape = gene_model_shape, 
                                     collapseTranscripts = "longest",cex.group = 1)
        
    #generate peak connection data frame
        
    peak <- gsub("\\.", "-", selected)
    connection_df <- data.frame()
    for(i in 1:nrow(gene_model)){
        gene <- gene_model[i, ]
        grange <- paste(chr, gene$start, gene$end, sep = "-")
        coaccess <- conn_input[gene$transcript, "score"]
        connection_df <- rbind(connection_df, c(grange, peak, as.numeric(coaccess)))
    }
    colnames(connection_df) <- c("Peak1", "Peak2", "coaccess")
    connection_df$peak_color <- "black"
    sub <- generate_plotting_subset(connection_df, chr, minbp, maxbp)
    sub$coaccess <- as.numeric(sub$coaccess)
    
    #get grang of peaks and generate dataTrack
    gr <- make_peak_track(sub)
    bk <- c(seq(-1, -0.1, by = 0.01),seq(0, 1, by = 0.01))
    dk <- c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
           colorRampPalette(colors = c("white","red"))(length(bk)/2))
    dtrack <- DataTrack(gr, type = c("heatmap"), chromosome = chr, gradient = dk, 
                   ylim = c(-1 , 1), yTicksAt = c(-1,-0.5,0,0.5), name = "scMAGeCK score")
    
    #rename sub
    
    if (!nrow(sub) == 0) {
        if (connection_color %in% names(sub)) {
          color_levs <- levels(as.factor(sub[,connection_color]))
          color_names <- rep("temp", length(color_levs))
          names(color_names) <- color_levs
          new_connection_color <- get_colors(sub[,connection_color])
          for(n in color_levs) {
              color_names[n] <-
                  new_connection_color[which(sub[,connection_color] == n)[1]]
          }
          connection_color <- new_connection_color
        }
        sub$color <- connection_color

        sub$width <- connection_width

        sub <- sub[,c("chr", "bp1", "bp2", "chr_2", "bp1_2", "bp2_2", "coaccess",
                      "width", "color")]
        names(sub) <- c("chrom1","start1","stop1","chrom2","start2","stop2",
                        "height", "width", "color")
      }else {
        warning("No connections above score_cutoff")
      }
    sub$height <- abs(sub$height)
    ctrack <- CustomTrack(plottingFunction = function(GdObject, prepare) {
        Gviz::displayPars(GdObject) <- list(ylim = c(0, max(abs(as.numeric(connection_df$coaccess)))))
        if(!prepare) {
          plotBedpe(sub, chrom = chr, chromstart = minbp, chromend = maxbp,
                    max(abs(as.numeric(connection_df$coaccess))), score_cut,
                    connection_width, alpha_by_coaccess,
                    color_names)
         }
        return(invisible(GdObject))}, name = "regulatory potential", fontsize.group = 6,fontsize = 12, , cex.title = 1.3)
        
    #in order to show all the gene names
        
    minbp <- minbp - 500000
    
    #get atrack
    
    atrack <- Gviz::GenomeAxisTrack(fontsize = 6, name = "")
    
    #plot
        
    if(include_axis_track == TRUE){
        gg <- as.ggplot(function()plotTracks(list(ctrack,  dtrack, atrack, grtrack), 
                                            title.width = 1.3, showTitle = TRUE, from = minbp, to = maxbp, 
                                            chromosome = chr, sizes = track_size, 
                                            transcriptAnnotation = "symbol", background.title = "transparent",
                                            col.border.title="transparent", lwd.border.title = "transparent",
                                            col.axis = "black", fontsize.group = 6, col.title="black",
                                            fontcolor.legend = "black"))
        gg <- gg + labs(title = selected) + 
                       theme(plot.title=element_text(hjust=0.5), text=element_text(size=12,face = "bold"))
    }else{
        gg <- as.ggplot(function()plotTracks(list(ctrack,  dtrack, grtrack), 
                                            title.width = 1.3, showTitle = TRUE, from = minbp, to = maxbp, 
                                            chromosome = chr, sizes = track_size, 
                                            transcriptAnnotation = "symbol", background.title = "transparent",
                                            col.border.title="transparent", lwd.border.title = "transparent",
                                            col.axis = "black", fontsize.group = 6, col.title="black",
                                            fontcolor.legend = "black"))
        gg <- gg + labs(title = selected) + 
                       theme(plot.title=element_text(hjust=0.5), text=element_text(size=12,face = "bold"))
    }
    return(gg)
}
    
generate_plotting_subset <- function(connections, chr, minbp, maxbp) {
  connections$Peak1 <- as.character(connections$Peak1)
  connections$Peak2 <- as.character(connections$Peak2)

  pcolor_map <- data.frame(Peak1 = connections$Peak1,
                           peak_color = connections$peak_color)
  pcolor_map <- pcolor_map[!duplicated(pcolor_map),]
  connections$peak_color <- NULL

  if(sum(!c("chr_1", "chr_2", "bp1_1", "bp2_1", "bp2_1", "bp2_2") %in%
         names(connections)) != 0 ) {
    suppressWarnings(connections$chr <- NULL)
    suppressWarnings(connections$bp1 <- NULL)
    suppressWarnings(connections$bp2 <- NULL)
    suppressWarnings(connections$chr_2 <- NULL)
    suppressWarnings(connections$bp1_2 <- NULL)
    suppressWarnings(connections$bp2_2 <- NULL)

    connections <- cbind(connections, df_for_coords(connections$Peak1)[,c(1, 2, 3)])
    cons2 <- df_for_coords(connections$Peak2) #slow
    cons2$Peak <- NULL
    names(cons2) <- c("chr_2", "bp1_2", "bp2_2")
    connections <- cbind(connections, cons2) #slow
  } else {
    if(!grepl("chr", connections$chr_1[1])) {
      connections$chr_1 <- paste0("chr", connections$chr_1)
    }
    if(!grepl("chr", connections$chr_2[1])) {
      connections$chr_2 <- paste0("chr", connections$chr_2)
    }
    names(connections)[names(connections) == "chr_1"] <- "chr"
    names(connections)[names(connections) == "bp1_1"] <- "bp1"
    names(connections)[names(connections) == "bp2_1"] <- "bp2"
  }

  sub <- connections[connections$chr_2 == chr & connections$bp1 <= maxbp &
                       connections$bp2 <= maxbp & connections$bp1 >= minbp &
                       connections$bp2 >= minbp & connections$bp1_2 <= maxbp &
                       connections$bp2_2 <= maxbp & connections$bp1_2 >= minbp &
                       connections$bp2_2 >= minbp,]


  sub <- sub[!duplicated(sub),]

  sub <- merge(sub, pcolor_map, all.x = TRUE)
  sub$peak_color <- as.character(sub$peak_color)
  sub$peak_color[is.na(sub$peak_color)] <- "black"

  return(sub)
}
                            
df_for_coords <- function(coord_strings) {
  coord_strings <- gsub(",", "", coord_strings)
  coord_cols <- as.data.frame(split_peak_names(coord_strings),
                              stringsAsFactors = FALSE )
  names(coord_cols) <- c("chr", "bp1", "bp2")
  coord_cols$Peak <- coord_strings
  coord_cols$bp1 <- as.numeric(coord_cols$bp1)
  coord_cols$bp2 <- as.numeric(coord_cols$bp2)
  coord_cols
}
                        
split_peak_names <- function(inp) {
  out <- stringr::str_split_fixed(stringi::stri_reverse(inp), 
                                  ":|-|_", 3)
  out[,1] <- stringi::stri_reverse(out[,1])
  out[,2] <- stringi::stri_reverse(out[,2])
  out[,3] <- stringi::stri_reverse(out[,3])
  out[,c(3,2,1), drop=FALSE]
}

make_peak_track <- function(df) {
  df <- df[!duplicated(df[,c("chr", "bp1", "bp2", "peak_color")]),]

  if (sum(duplicated(df[,c("chr", "bp1", "bp2")])) > 0)
    stop(paste("Multiple peak colors correspond to a single peak. Be sure that",
               "your peak_color column name assigns colors for Peak1 only",
               collapse = " "))
    
  df2 <- df[!duplicated(df[,"Peak1"]), ]
  gr <-  GenomicRanges::GRanges(as.character(df2$chr),
                 IRanges::IRanges(as.numeric(as.character(df2$bp1)),
                         as.numeric(as.character(df2$bp2))), score = df2$coaccess)
  return(gr)
}
    
plotBedpe <- function(bedpedata,
                      chrom,
                      chromstart,
                      chromend,
                      ymax,
                      score_cut,
                      width,
                      alpha_by_coaccess,
                      color_names = NULL)
{ ###### All borrowed and modified from Sushi package.

  if (nrow(bedpedata) == 0) {
    warning("Nothing to plot")
    return()
  }

  bedpedata  <- bedpedata[,c("chrom1","start1","stop1","chrom2","start2",
                             "stop2", "height", "width", "color")]

  # normalize height
  maxheight <- ymax

  bedpedata$alpha <- .6
  if(alpha_by_coaccess) {
    bedpedata$alpha <- (bedpedata$height-score_cut)/maxheight
  }
  bedpedata$height <- bedpedata$height/maxheight
  # remove any rows with 0 height
  bedpedata <- bedpedata[abs(bedpedata$height) > 0,]

  # reclass data
  if (any(class(bedpedata) == "data.table")) {
    for(i in c("start1", "stop1", "start2", "stop2")) {
      bedpedata[[i]] <- as.numeric(as.character((bedpedata[[i]])))
    }
  } else {
    for(i in c("start1", "stop1", "start2", "stop2")) {
      bedpedata[,i] <- as.numeric(as.character((bedpedata[,i])))
    }
  }

  # add position columns
  bedpedata$pos1 = apply(bedpedata[,c("start1","stop1")],1,mean)
  bedpedata$pos2 = apply(bedpedata[,c("start2","stop2")],1,mean)

  totalrange <- as.numeric(as.character(chromend)) -
    as.numeric(as.character(chromstart))
  if (nrow(bedpedata) == 0) warning("Nothing to plot")

  #legFactors <- sort(names(which(apply(legInfo, 2, any))))
  #boxSize <-  if(length(setdiff(legFactors, c("col", "cex")))==0) 0.1 else 0.3
  #pcols <- Gviz:::.getPlottingFeatures(GdObject)

  if (!is.null(color_names)) {
    boxSize <- .3
    spacing <- 0.2
    vspace <- .05
    for (i in seq_len(length(color_names))) {
      grid::grid.lines(unit(c(spacing,spacing + boxSize), "inches"),
                       c(1 - vspace*i, 1 - vspace*i),
                       gp=grid::gpar(col=color_names[i], lwd=width))
      grid::grid.text(x=unit(.1 + (boxSize + spacing), "inches"),
                      y=1 - vspace*i, just=c(0, 0.5),
                      label=names(color_names)[i])
    }
  }

  # plot the data
  grid::grid.function(function(x) list(x=x, y=(score_cut/(ymax))), 
                      gp=grid::gpar(col="black", lty="dashed", lwd=width)) #
  for (row in (seq_len(nrow(bedpedata)))) {
    x1     = bedpedata$pos1[row]
    x2     = bedpedata$pos2[row]
    height = bedpedata$height[row]
    width  = bedpedata$width[row]
    color  = bedpedata$color[row]
    alpha  = bedpedata$alpha[row]
    plotpair(x1,x2,height,totalrange,width,color, chromstart, alpha)
  }
}

plotpair <- function(start, end, height, totalrange,
                     width, color, chromstart, alpha) {
  #scale values for plotting
  x1 = (min(start,end) - as.numeric(as.character(chromstart)))/totalrange
  x2 = (max(start,end) - as.numeric(as.character(chromstart)))/totalrange
  hx1 <- (x1 + x2)/2
  hy1 <- height/.725

  grid::grid.bezier(x = c(x1, hx1, hx1, x2), y = c(0, hy1, hy1, 0),
                    default.units = "npc",
                    gp=grid::gpar(col=color, lwd=width,
                                  alpha = (alpha*.9 + .1),fontsizecex = 10))
}

#' @export    
    
ATACciceroPlot<- function(object, score_dir, pval_dir, species = "Hs", version = "v75", gene_annotations = NULL, 
                          pro_annotations = NULL, pro_up = 2000, pro_down = 200, overlap_cut = 0, p_val_cut = 0.05,
                          score_cut = 0, p_adj_cut = 0.05, logFC_cut = 1, selected, NTC = "NTC", min.pct = 0.2,
                          upstream = 2000000, downstream = 2000000, test.use = "wilcox", track_size = c(1,.3,.2,.3), 
                          include_axis_track = TRUE, connection_color = "#7F7CAF", connection_color_legend = TRUE,
                          connection_width = 2,connection_ymax = NULL, gene_model_color = "#81D2C7", 
                          alpha_by_coaccess = FALSE, gene_model_shape = c("smallArrow", "box")){
    
    color_names = NULL
    
    #get peak matrix
    
    if(is.character(object)){
        peak <- readRDS(object)
    }else{
        peak <- object
    }
    
    #rename ident
    
    if("perturbations" %in% colnames(peak@meta.data)){
        if(selected %in% peak$perturbations){
            peak@active.ident<- peak$perturbations
        }else{
            stop(paste("Cannot find ", selected, " in perturbations", sep = "'"))
        }
    }else{
        stop("Cannot find 'perturbations' in object, please renamed it or using 'Add_meta_data' function to add it")
    }
    
    #find DA peaks
    
    peak<- RunTFIDF(peak)
    da_peaks <- FindMarkers(
        object = peak,
        ident.1 = selected,
        ident.2 = NTC,
        min.pct = min.pct,
        test.use = test.use
        )
    da_peak <- subset(da_peaks, p_val_adj <= p_adj_cut & (avg_log2FC >= logFC_cut | avg_log2FC <= -logFC_cut))
    if(nrow(da_peak) == 0){
        stop("Cannot find DA peaks pass threshold")
    }
    
    da_peak$chromosome<- t(data.frame(strsplit(rownames(da_peak), "[.|:|-]")))[ ,1]
    da_peak$start<- t(data.frame(strsplit(rownames(da_peak), "[.|:|-]")))[ ,2]
    da_peak$end<- t(data.frame(strsplit(rownames(da_peak), "[.|:|-]")))[ ,3]
    
    # get genome annotations
    
    if(is.null(gene_annotations)){
        if(species == "Hs"){
            if(version == "v75"){
                a <- genes(EnsDb.Hsapiens.v75)
            }else if(version == "v79"){
                a <- genes(EnsDb.Hsapiens.v79)
            }else if(version == "v86"){
                a <- genes(EnsDb.Hsapiens.v86)
            }  
        }else if(species == "Mm"){
            if(version == "v75"){
                a <- genes(EnsDb.Mmusculus.v75)
            }else if(version == "v79"){
                a <- genes(EnsDb.Mmusculus.v79)
            }
        }
        gene_anno<- data.frame(a)
        gene_anno$chromosome <- paste0("chr", gene_anno$seqnames)
        gene_anno$transcript <- gene_anno$symbol
    }else{
        gene_anno <- gene_annotations
            }
    
    #get promoters annotations
        
    if(is.null(pro_annotations)){
        if(species == "Hs"){
            if(version == "v75"){
                b <- promoters(EnsDb.Hsapiens.v75)
            }else if(version == "v79"){
                b <- promoters(EnsDb.Hsapiens.v79)
            }else if(version == "v86"){
                b <- promoters(EnsDb.Hsapiens.v86)
            }  
        }else if(species == "Mm"){
            if(version == "v75"){
                b <- promoters(EnsDb.Mmusculus.v75)
            }else if(version == "v79"){
                b <- promoters(EnsDb.Mmusculus.v79)
            }
        }
        pro_anno<- data.frame(b)
        pro_anno$chromosome <- paste0("chr", pro_anno$seqnames)
    }else{
        pro_anno <- pro_annotations
            }
    
    #get enhancer list
        
    enhancer_list<- data.frame()
    for(i in 1:nrow(da_peak)){
        chr <- da_peak[i, "chromosome"]
        start <- as.numeric(da_peak[i, "start"])
        end <- as.numeric(da_peak[i, "end"])
        minbp <- start - 2*(pro_up + pro_down)
        maxbp <- end + 2*(pro_up + pro_down)
        peak_model <- pro_anno
        peak_model <- peak_model[!is.na(peak_model$chromosome) & 
                             !is.na(peak_model$start) &
                             !is.na(peak_model$end) &
                             !is.na(peak_model$strand),]
        peak_model <- peak_model[peak_model$chromosome == chr &
                   ((peak_model$start > minbp & peak_model$start < maxbp) |
                    (peak_model$end > minbp & peak_model$end < maxbp) |
                    (peak_model$start < minbp & peak_model$end > maxbp)),]
        if(nrow(peak_model) == 0){
            enhancer_list <- rbind(enhancer_list, c(chr, start, end))
        }else{
            pro_seq <- apply(peak_model[,c("start", "end")], 1, function(x){seq(x[1], x[2])})
            overlap <- length(intersect(seq(start, end), pro_seq))
            if(overlap <= overlap_cut){
                enhancer_list <- rbind(enhancer_list, c(chr, start, end))
            }
        }
    }
    if(nrow(enhancer_list != 0)){
        colnames(enhancer_list) <- c("chromosome", "start", "end")
        rownames(enhancer_list) <- paste(enhancer_list$chromosome, enhancer_list$start, enhancer_list$end, sep = "-")
    }else{
        enhancer_list <- da_peak[, c("chromosome", "start", "end")]
        warning("All DA peaks has overlap with promoters, using all DA peaks passed threshold as input list")
    }
    
    #get score and p-value
    
    if(is.character(score_dir)){
        score <- read.table(score_dir, header = T, row.names = 1)
    }else{
        score <- score_dir
    }
    if(is.character(pval_dir)){
        pval <- read.table(pval_dir, header = T, row.names = 1)
    }else{
        pval <- pval_dir
    }    
    
    #get score and p-value of selected enhancer
    
    conn_input <- data.frame(score = score[, selected], pval = pval[, selected])
    rownames(conn_input) <- rownames(score)
    #conn_input <- conn_input[rownames(conn_input) %in% gene_model$transcript, ]
    conn_input <- subset(conn_input, ((score > score_cut) | (score < -score_cut)) & pval < p_val_cut)    
    if(nrow(conn_input) == 0){
        stop("Cannot find genes passed threshold")
    }
    
    #create store list
    
    results <- list()
        
    #get results for all candidate enhancer
    
    for(i in 1:nrow(enhancer_list)){
        chr <- enhancer_list[i, "chromosome"]
        start <- as.numeric(enhancer_list[i, "start"])
        end <- as.numeric(enhancer_list[i, "end"])
        
        #extend region
    
        minbp<- start - upstream
        maxbp<- end + downstream
    
        #get gene model

        gene_model <- gene_anno
        gene_model <- gene_model[!is.na(gene_model$chromosome) & 
                                 !is.na(gene_model$start) &
                                 !is.na(gene_model$end) &
                                 !is.na(gene_model$strand) &
                                 !is.na(gene_model$transcript), ]
        gene_model <- gene_model[gene_model$chromosome == chr &
                       ((gene_model$start > minbp & gene_model$start < maxbp) |
                        (gene_model$end > minbp & gene_model$end < maxbp) |
                        (gene_model$start < minbp & gene_model$end > maxbp)), ]
        gene_model <- gene_model[gene_model$transcript %in% rownames(score), ]
    

        gene_model <- subset(gene_model, transcript %in% rownames(conn_input))
        grtrack <- Gviz::GeneRegionTrack(gene_model, chromosome = chr, geneSymbols = TRUE,
                                         name = "", fill = "#81D2C7",
                                         col = "#81D2C7",  fontcolor = "black",
                                         fontcolor.group = "black", fontsize.group = 6,
                                         fontsize = 6, shape = gene_model_shape, 
                                         collapseTranscripts = "longest",cex.group = 1)
        
        #generate peak connection data frame
        
        peak <- paste(chr, start, end, sep = "-")
        connection_df <- data.frame()
        if(nrow(gene_model) == 0){
            results[[i]] <- "No results"
            names(results)[i] <- peak
            next
        }
        for(j in 1:nrow(gene_model)){
            gene <- gene_model[j, ]
            grange <- paste(chr, gene$start, gene$end, sep = "-")
            coaccess <- conn_input[gene$transcript, "score"]
            connection_df <- rbind(connection_df, c(grange, peak, as.numeric(coaccess)))
        }
        colnames(connection_df) <- c("Peak1", "Peak2", "coaccess")
        connection_df$peak_color <- "black"
        sub <- generate_plotting_subset(connection_df, chr, minbp, maxbp)
        sub$coaccess <- as.numeric(sub$coaccess)
    
        #get grang of peaks and generate dataTrack
        gr <- make_peak_track(sub)
        bk <- c(seq(-1, -0.1, by = 0.01),seq(0, 1, by = 0.01))
        dk <- c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
               colorRampPalette(colors = c("white","red"))(length(bk)/2))
        dtrack <- DataTrack(gr, type = c("heatmap"), chromosome = chr, gradient = dk, 
                            ylim = c(-1 , 1), yTicksAt = c(-1,-0.5,0,0.5), name = "scMAGeCK score")
    
        #rename sub
    
        if (!nrow(sub) == 0) {
            if (connection_color %in% names(sub)) {
              color_levs <- levels(as.factor(sub[,connection_color]))
              color_names <- rep("temp", length(color_levs))
              names(color_names) <- color_levs
              new_connection_color <- get_colors(sub[,connection_color])
              for(n in color_levs) {
                  color_names[n] <-
                      new_connection_color[which(sub[,connection_color] == n)[1]]
              }
              connection_color <- new_connection_color
            }
            sub$color <- connection_color

            sub$width <- connection_width

            sub <- sub[,c("chr", "bp1", "bp2", "chr_2", "bp1_2", "bp2_2", "coaccess",
                      "width", "color")]
            names(sub) <- c("chrom1","start1","stop1","chrom2","start2","stop2",
                        "height", "width", "color")
          }else {
            warning("No connections above score_cutoff")
          }
        sub$height <- abs(sub$height)
        ctrack <- CustomTrack(plottingFunction = function(GdObject, prepare) {
            Gviz::displayPars(GdObject) <- list(ylim = c(0, max(abs(as.numeric(connection_df$coaccess)))))
            if(!prepare) {
              plotBedpe(sub, chrom = chr, chromstart = minbp, chromend = maxbp,
                        max(abs(as.numeric(connection_df$coaccess))), score_cut,
                        connection_width, alpha_by_coaccess,
                        color_names)
             }
            return(invisible(GdObject))}, name = "regulatory potential", fontsize.group = 6,fontsize = 12, , cex.title = 1.3)
        
        #in order to show all the gene names
        
        minbp <- minbp - 500000
    
        #get atrack
    
        atrack <- Gviz::GenomeAxisTrack(fontsize = 6, name = "")
        
        #plot
        
        if(include_axis_track == TRUE){
        gg <- as.ggplot(function()plotTracks(list(ctrack,  dtrack, atrack, grtrack), 
                                            title.width = 1.3, showTitle = TRUE, from = minbp, to = maxbp, 
                                            chromosome = chr, sizes = track_size, 
                                            transcriptAnnotation = "symbol", background.title = "transparent",
                                            col.border.title="transparent", lwd.border.title = "transparent",
                                            col.axis = "black", fontsize.group = 6, col.title="black",
                                            fontcolor.legend = "black"))
        gg <- gg + labs(title = paste(selected, peak, sep = ":")) + 
                       theme(plot.title=element_text(hjust=0.5), text=element_text(size=12,face = "bold"))
        }else{
        gg <- as.ggplot(function()plotTracks(list(ctrack,  dtrack, grtrack), 
                                            title.width = 1.3, showTitle = TRUE, from = minbp, to = maxbp, 
                                            chromosome = chr, sizes = track_size, 
                                            transcriptAnnotation = "symbol", background.title = "transparent",
                                            col.border.title="transparent", lwd.border.title = "transparent",
                                            col.axis = "black", fontsize.group = 6, col.title="black",
                                            fontcolor.legend = "black"))
        gg <- gg + labs(title = paste(selected, peak, sep = ":")) + 
                       theme(plot.title=element_text(hjust=0.5), text=element_text(size=12,face = "bold"))
        }
        results[[i]] <- gg
        names(results)[i] <- peak
        }
        return(results)
}