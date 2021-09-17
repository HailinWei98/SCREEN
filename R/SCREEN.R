#' SCREEN
#' @export

SCREEN <- function(sg_dir, mtx_dir, fragments = NULL, species = "Hs", version = "v75", data_type = "RNA", 
                   Mixscape = TRUE, prefix = "./", label = "", gene_type = "Symbol", protein_coding = TRUE, 
                   frac = 0.01, cal.mt = TRUE, nFeature = c(200, 5000), nCount = 1000, FRiP = 0.1,
                   mt = 10, blank_NTC = FALSE, lambda = 0.01, permutation = NULL,
                   p_val_cut = 0.05, score_cut = 0.5, cicero_p_val_cut = 0.05, cicero_score_cut = 0,
                   ylimit = c(-600, 600, 200), project = "perturb", NTC = "NTC", replicate = 1, 
                   select_gene = NULL, selected =  NULL, gene_annotations = NULL, pro_annotations = NULL, 
                   pro_up = 3000, pro_down = 0, overlap_cut = 0, p_adj_cut = 0.05, logFC_cut = 1,
                   min.pct = 0.2, upstream = 2000000, downstream = 2000000, test.use = "wilcox", 
                   track_size = c(1,.3,.2,.3), include_axis_track = TRUE, connection_color = "#7F7CAF",
                   connection_color_legend = TRUE, connection_width = 2, connection_ymax = NULL,
                   gene_model_color = "#81D2C7", alpha_by_coaccess = FALSE,
                   gene_model_shape = c("smallArrow", "box")) {
    #get matrix
    
    if (is.character(mtx_dir)) {
        message(paste("Reading RDS file:", mtx_dir))
        mtx <- readRDS(mtx_dir)
    } else {
        mtx <- mtx_dir
    }
    
    #get sgRNA library
    
    if (is.character(sg_dir)) {
        message(paste("Reading sgRNA lib file:", sg_dir))
        sg_lib <- read.table(sg_dir, header =T )
    } else {
        sg_lib <- sg_dir
    }
    
    #for RNA input
    
    if (data_type == "RNA") {
        mtx <- Add_meta_data(sg_lib, mtx, cal.mt, species, replicate)
        mtx_QC <- scQC(mtx, prefix, label, species, frac, nFeature, nCount, mt, blank_NTC)
        p <- sgRNA_quality_plot(sg_lib, mtx,
                               LABEL = paste(label, "before_QC", sep = ""), prefix)
        q <- sgRNA_quality_plot(sg_lib, mtx_QC,
                                LABEL = paste(label, "after_QC", sep = ""), prefix)
        if(Mixscape == TRUE){
            message("Running Mixscape")
            mixscape <- IntegratedMixscape(sg_lib, mtx_QC, NTC, prefix, label)
        }
        mtx_QC <- normalize_scale(mtx_QC)
        scaled <- GetAssayData(object = mtx_QC, slot = "scale.data")
        
        message("Save RDS file")
        saveRDS(mtx, file = file.path(prefix, "perturb.rds"))
        saveRDS(mtx_QC, file = file.path(prefix, "perturb_QC.rds"))
        
        message("Save sgRNA library to save path")
        write.table(sg_lib, file = file.path(prefix, "sg_lib.txt"),
                    row.names = FALSE, quote = FALSE)
        
        message("Running scMAGeCK")
        result <- improved_scmageck_lr(sg_lib, scaled, NTC, select_gene,
                                        LABEL = paste(label, "improved",sep = ""),
                                        permutation, prefix, lambda)
        score <- result[[1]]
        p_val <- result[[2]]
        
        message("Running DE_gene_plot")
        y <- DE_gene_plot(score, p_val, project,
                          prefix, label, p_val_cut, score_cut,
                          ylimit)
        if(length(grep("^chr", sg_lib$gene)) != 0){
            message("Detect enhancer perturbation in sgRNA library, running ciceroPlot")
            cicero_result <- ciceroPlot(score, p_val, selected, species, versions, gene_annotations, 
                                        cicero_p_val_cut, cicero_score_cut, upstream, downstream, 
                                        track_size, include_axis_track, connection_color,
                                        connection_color_legend, connection_width, connection_ymax, 
                                        gene_model_color, alpha_by_coaccess, gene_model_shape)
            message("Save cicero_results")
            dir.create(file.path(prefix, "cicero_results"))
            path <- file.path(prefix, "cicero_results")
            if(is.list(cicero_result)){
                for(i in 1:length(cicero_result)){
                    a <- cicero_result[i]
                    for(j in 1:length(a)){
                        pdf(paste(file.path(path, names(a)[j]), ".pdf", sep = ""))
                        print(a[j])
                        dev.off()
                    }
                }
            }else{
                pdf(paste(file.path(path, selected), ".pdf", sep = ""))
                print(cicero_result)
                dev.off()
            }
        } 
        results <- list()
        results[1] <- result
        results[2] <- y
        results[3] <- cicero_result
        names(results) <- c("scMAGeCK_lr", "DE_gene_plot", "cicero_results")
    } else if(data_type == "ATAC"){
        mtx <- ATAC_Add_meta_data(sg_lib, mtx, fragments, replicate)
        mtx_QC <- ATAC_scQC(mtx, prefix, label, frac, nFeature, nCount, FRiP, blank_NTC)
        p <- sgRNA_quality_plot(sg_lib, mtx,
                               LABEL = paste(label, "before_QC", sep = ""), prefix)
        q <- sgRNA_quality_plot(sg_lib, mtx_QC,
                                LABEL = paste(label, "after_QC", sep = ""), prefix)
        
        message("Calculate gene activity")
        RNA <- CalculateGeneActivity(mtx, fragments, species, version, gene_type, protein_coding, pro_up, pro_down)
        
        if(Mixscape == TRUE){
            message("Running Mixscape")
            mixscape <- IntegratedMixscape(sg_lib, RNA, NTC, prefix, label)
        }
        RNA <- normalize_scale(RNA)
        scaled <- GetAssayData(object = RNA, slot = "scale.data")
        
        message("Save RDS file")
        saveRDS(mtx, file = file.path(prefix, "peak.rds"))
        saveRDS(RNA, file = file.path(prefix, "RNA.rds"))
        
        message("Save sgRNA library to save path")
        write.table(sg_lib, file = file.path(prefix, "sg_lib.txt"),
                    row.names = FALSE, quote = FALSE)
        
        message("Running scMAGeCK")
        result <- improved_scmageck_lr(sg_lib, scaled, NTC, select_gene,
                                        LABEL = paste(label, "improved",sep = ""),
                                        permutation, prefix, lambda)
        score <- result[1]
        p_val <- result[2]
        
        message("Running DE_gene_plot of gene activity matrix")
        y <- DE_gene_plot(score, p_val, project,
                          prefix, label, p_val_cut, score_cut,
                          ylimit)
        message("Detect enhancer perturbation in sgRNA library, running ATACciceroPlot")
        cicero_result <- ATACciceroPlot(mtx_QC, score, p_val, selected, species, versions, gene_annotations, 
                                        pro_annotations, pro_up, pro_down, overlap_cut, cicero_p_val_cut, 
                                        cicero_score_cut, p_adj_cut, logFC_cut, NTC, min.pctupstream, downstream, 
                                        test.use, track_size, include_axis_track, connection_color,
                                        connection_color_legend, connection_width, connection_ymax, 
                                        gene_model_color, alpha_by_coaccess, gene_model_shape)
        
        message("Save cicero_results")
        dir.create(file.path(prefix, "cicero_results"))
        path <- file.path(prefix, "cicero_results")
        if(is.list(cicero_result)){
            for(i in 1:length(cicero_result)){
                a <- cicero_result[i]
                for(j in 1:length(a)){
                    pdf(paste(file.path(path, names(a)[j]), ".pdf", sep = ""))
                    print(a[j])
                    dev.off()
                }
            }
        }
        results <- list()
        results[1] <- result
        results[2] <- y
        results[3] <- cicero_result
        names(results) <- c("scMAGeCK_lr", "DE_gene_plot", "cicero_results")
    }
    return(results)
}
