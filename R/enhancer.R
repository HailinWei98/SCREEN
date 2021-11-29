#' @import Seurat
#' @export

EnhancerGeneExpression <- function(sg_dir, mtx_dir, selected = NULL, prefix = "./", upstream = 2000000, 
                                   downstream = 2000000, gene_annotations = NULL, species = "Hs", version = "v75",
                                   NTC = "NTC") {
    
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
    if (!dir.exists(dir)) {
        dir.create(path = dir)
    }
    
    img_dir <- file.path(prefix, "img")
    if (!dir.exists(img_dir)) {
        dir.create(path = img_dir)
    }
    img_dir <- file.path(img_dir, "enhancer_gene_expression")
    if (!dir.exists(img_dir)) {
        dir.create(path = img_dir)
    }
    
    if(is.character(selected)){
        if(length(selected) == 1){
                                  
            #get information from selected enhancer
    
            all <- unlist(strsplit(selected, "[.|:|-]"))
            chr <- all[1]
            start <- as.numeric(all[2])
            end <- as.numeric(all[3])
        
            #extend region
    
            minbp <- start - upstream
            maxbp <- end + downstream
            
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
            gene_model <- gene_model[gene_model$transcript %in% rownames(mtx), ]
            
            if(nrow(gene_model) == 0){
                stop("Cannot find genes close to selected region.")
            }
            gene_list <- gene_model$transcript
            
            results <- list()
            select_sg <- subset(sg_lib, gene == selected)
            
            #rename enhancer
            
            new_select <- gsub("[:|-]", ".", selected)
            
            mtx_new <- subset(mtx, cells = c(select_sg$cell, colnames(subset(mtx, perturbations == NTC))))
            mtx_new$perturbations <- gsub(selected,
                                          "single enhancer",
                                          mtx_new$perturbations)
            mtx_new$perturbations <- gsub("multiple", "multiple enhancer", mtx_new$perturbations)
            mtx_new$perturbations <- factor(mtx_new$perturbations, 
                                            levels = c("single enhancer", 
                                                       "multiple enhancer", NTC))
            mtx_new@active.ident<- mtx_new$perturbations           
            
            img_dir <- file.path(img_dir, new_select)
            if (!dir.exists(img_dir)) {
                dir.create(path = img_dir)
            }
            
            for(i in 1:nrow(gene_model)) {
                gene <- gene_list[i]
                p <- VlnPlot(object = mtx_new,
                             features = gene,
                             pt.size = 0.1) +
                theme(plot.title = element_text(hjust = 0.5, size = 25),
                      axis.text.x = element_text(size = 25),
                      axis.title.x = element_text(size = 22), 
                      axis.title.y = element_text(size = 22),
                      axis.text.y = element_text(hjust = 0.5, size = 25)) + 
                NoLegend()
                results[[i]] <- p
                names(results[i]) <- gene
            }
            
            #save plot

            pdf(file = file.path(prefix, paste(new_select, ".pdf", sep = "")), height = 20, width = 18)
            for(i in 1:(ceiling(nrow(gene_model)/12))){
                if(i == ceiling(nrow(gene_model)/12)) {
                    up <- nrow(gene_model)
                } else {
                    up <- i * 12
                }
                p <- plot_grid(plotlist = results[(i * 12 - 11) : up], 
                               ncol = 4, align = "hv") + 
                theme(plot.margin = margin(20, 20, 20, 40))
                p <- ggpubr::annotate_figure(p, left = ggpubr::text_grob(new_select,
                                                                         size = 40, 
                                                                         face = "bold", 
                                                                         rot = 90,
                                                                         vjust = 1)) + 
                theme(plot.margin = margin(20, 20, 20, 30))
                print(p)
            }
            dev.off()
        
            for(i in 1:(ceiling(nrow(gene_model)/12))){
                if(i == ceiling(nrow(gene_model)/12)) {
                    up <- nrow(gene_model)
                } else {
                    up <- i * 12
                }
                p <- plot_grid(plotlist = results[(i * 12 - 11) : up], 
                               ncol = 4, align = "hv") + 
                theme(plot.margin = margin(20, 20, 20, 40))
                p <- ggpubr::annotate_figure(p, left = ggpubr::text_grob(new_select, 
                                                                         size = 40,
                                                                         face = "bold", 
                                                                         rot = 90, 
                                                                         vjust = 1)) + 
                theme(plot.margin = margin(20, 20, 20, 30))
                png(file.path(img_dir, paste(i, ".png", sep = "")), 
                    height = 1400, width = 1300, res = 72)
                print(p)
                dev.off()
            }
            
        }else{
                
            #get results for all perturbations
                
            j <- 0
            results <- list()
            for(perturb in selected){
                
                #get information from selected enhancer
                
                all <- unlist(strsplit(perturb, "[.|:|-]"))
                chr <- all[1]
                start <- as.numeric(all[2])
                end <- as.numeric(all[3])
        
                #extend region
    
                minbp <- start - upstream
                maxbp <- end + downstream
                
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
                gene_model <- gene_model[gene_model$transcript %in% rownames(mtx), ]
            
                if(nrow(gene_model) == 0){
                    next
                }
                gene_list <- gene_model$transcript
            
                #rename enhancer
                
                new_perturb <- gsub("[:|-]", ".", perturb)
                
                result <- list()
                select_sg <- subset(sg_lib, gene == perturb)
                
                mtx_new <- subset(mtx, cells = c(select_sg$cell, colnames(subset(mtx, perturbations == NTC))))
                mtx_new$perturbations <- gsub(perturb,
                                              "single enhancer",
                                              mtx_new$perturbations)
                mtx_new$perturbations <- gsub("multiple", "multiple enhancer", mtx_new$perturbations)
                mtx_new$perturbations <- factor(mtx_new$perturbations, 
                                                levels = c("single enhancer", 
                                                           "multiple enhancer", NTC))
                mtx_new@active.ident<- mtx_new$perturbations
                
                #create save path
                
                dir <- file.path(img_dir, new_perturb)
                if (!dir.exists(dir)) {
                    dir.create(path = dir)
                }
                
                for(i in 1:nrow(gene_model)) {
                    gene <- gene_list[i]
                    p <- VlnPlot(object = mtx_new,
                                 features = gene,
                                 pt.size = 0.1) +
                    theme(plot.title = element_text(hjust = 0.5, size = 25),
                          axis.text.x = element_text(size = 25),
                          axis.title.x = element_text(size = 22), 
                          axis.title.y = element_text(size = 22),
                          axis.text.y = element_text(hjust = 0.5, size = 25)) + 
                    NoLegend()
                    
                    result[[i]] <- p
                    names(result[i]) <- gene
                
                }
                
                #save plot

                pdf(file = file.path(prefix, paste(new_perturb, ".pdf", sep = "")), height = 20, width = 18)
                for(i in 1:(ceiling(nrow(gene_model)/12))){
                    if(i == ceiling(nrow(gene_model)/12)) {
                        up <- nrow(gene_model)
                    } else {
                        up <- i * 12
                    }
                    p <- plot_grid(plotlist = result[(i * 12 - 11) : up], 
                                   ncol = 4, align = "hv") + 
                    theme(plot.margin = margin(20, 20, 20, 40))
                    p <- ggpubr::annotate_figure(p, left = ggpubr::text_grob(new_perturb,
                                                                             size = 40, 
                                                                             face = "bold", 
                                                                             rot = 90,
                                                                             vjust = 1)) + 
                theme(plot.margin = margin(20, 20, 20, 30))
                print(p)
            }
            dev.off()
        
            for(i in 1:(ceiling(nrow(gene_model)/12))){
                if(i == ceiling(nrow(gene_model)/12)) {
                    up <- nrow(gene_model)
                } else {
                    up <- i * 12
                }
                p <- plot_grid(plotlist = result[(i * 12 - 11) : up], 
                               ncol = 4, align = "hv") + 
                theme(plot.margin = margin(20, 20, 20, 40))
                p <- ggpubr::annotate_figure(p, left = ggpubr::text_grob(new_perturb, 
                                                                         size = 40,
                                                                         face = "bold", 
                                                                         rot = 90, 
                                                                         vjust = 1)) + 
                theme(plot.margin = margin(20, 20, 20, 30))
                png(file.path(dir, paste(i, ".png", sep = "")), 
                    height = 1400, width = 1300, res = 72)
                print(p)
                dev.off()
            }
                
                j <- j + 1
                results[[j]] <- result
                names(results)[j] <- new_perturb
            }
        }
        return(results)
    }else{
        stop("Please input correct format of selected perturbations")
    }
}