#' @import Seurat
#' @export

EnhancerGeneExpression <- function(sg_dir, mtx_dir, selected = NULL) {
    
    #read files
    
    if (is.character(mtx_dir)) {
        message(paste("Reading RDS file:", mtx_dir))
        mtx <- readRDS(mtx_dir)
    } else {
        mtx <- mtx_dir
    }
    if (is.character(sg_dir)) {
        message(paste("Reading sgRNA lib file:", sg_dir))
        sg_lib <- read.table(sg_dir, header = T)
    } else {
        sg_lib <- sg_dir
    }

    if(is.null(selected)){
        selected <- sg_lib$gene
        selected <- selected[grep("^chr", selected)]
    }
    
    dir <- file.path(prefix, "enhancer_gene_expression")
    dir.create(path = dir)
    
    if(is.character(selected)){
        if(length(selected) == 1){
                                  
            #get score and p-value of selected enhancer
                
            conn_input <- data.frame(score = score[, selected], pval = pval[, selected])
            rownames(conn_input) <- rownames(score)
            conn_input <- subset(conn_input, ((score > score_cut) | (score < -score_cut)) & pval < p_val_cut)
            if(nrow(conn_input) == 0){
                stop(paste("Cannot find genes passed threshold in ", selected, " scMAGeCK results", sep = "'"))
            }
            #get information from selected enhancer
    
            all <- unlist(strsplit(selected, "[.|:|-]"))
            chr <- all[1]
            start <- as.numeric(all[2])
            end <- as.numeric(all[3])
        
            #extend region
    
            minbp <- start - upstream
            maxbp <- end + downstream
            
            #get results
            
            gg <- get_results(chr, start, end, minbp, maxbp, gene_anno, track_size, gene_model_shape,
                              conn_input, connection_color, include_axis_track, score, 
                              score_cut, connection_width, alpha_by_coaccess, color_names)
            if(is.null(gg)){
                stop("Cannot find gene model close to selected region.")
            }
            gg <- gg + labs(title = selected) + 
                            theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 12,face = "bold"))
            
            #save plot
            
            pdf(file = file.path(dir, paste(selected, ".pdf", sep = "")))
            print(gg)
            dev.off()
        
            png(file.path(file.path(prefix, "img/cicero"), paste(selected, ".png", sep = "")), 
                width = 600 * 3, height = 600 * 3, res = 72 * 3)
            print(gg)
            dev.off()
            
            return(gg)
        }else{
                
            #get results for all perturbations
                
            j <- 0
            results <- list()
            for(perturb in selected){
                
                #get score and p-value of selected enhancer
                
                conn_input <- data.frame(score = score[, perturb], pval = pval[, perturb])
                rownames(conn_input) <- rownames(score)
                conn_input <- subset(conn_input, ((score > score_cut) | (score < -score_cut)) & pval < p_val_cut)
                if(nrow(conn_input) == 0){
                    warning(paste("Cannot find genes passed threshold in ", perturb, " scMAGeCK results", sep = "'"))
                        next
                }
                
                #get information from selected enhancer
                
                all <- unlist(strsplit(perturb, "[.|:|-]"))
                chr <- all[1]
                start <- as.numeric(all[2])
                end <- as.numeric(all[3])
        
                #extend region
    
                minbp <- start - upstream
                maxbp <- end + downstream
            
                #get results
            
                gg <- get_results(chr, start, end, minbp, maxbp, gene_anno, track_size, gene_model_shape,
                                  conn_input, connection_color, include_axis_track, score, 
                                  score_cut, connection_width, alpha_by_coaccess, color_names)
                if(is.null(gg)){
                    next
                }
                j <- j + 1
                gg <- gg + labs(title = perturb) + 
                theme(plot.title = element_text(hjust = 0.5), 
                      text = element_text(size = 12, face = "bold"))
                results[[j]] <- gg
                names(results)[j] <- perturb
                
                #save plot
            
                pdf(file = file.path(dir, paste(perturb, ".pdf", sep = "")))
                print(gg)
                dev.off()
        
                png(file.path(file.path(prefix, "img/cicero"), paste(perturb, ".png", sep = "")), 
                    width = 600 * 3, height = 600 * 3, res = 72 * 3)
                print(gg)
                dev.off()
            }
            return(results)
        }
    }else{
        stop("Please input correct format of selected perturbations")
    }
    p1 <- VlnPlot(
        object = mtx,
        features = c("MAGI3", "RBM15", "DRAM2"),
        pt.size = 0.1,
        idents = c("chr14:55569715-55569738", "NTC")
    )
}