#' @export

config_generation <- function(mtx, mtx_QC, sg_lib, score, p_val, project, prefix = "./", label = "", 
                              species = "Hs", version = "v75", type = "ATAC", NTC = "NTC", article = "", data = "",
                              article_name = "", data_name = "",gene_type = "Symbol", p_val_cut = 0.05, 
                              score_cut = 0.5, DA = NULL, cicero = NULL, enhancer = NULL) {
    
    if (is.character(mtx)) {
        message(paste("Reading RDS file:", mtx))
        mtx <- readRDS(mtx)
    }
    
    if (is.character(mtx_QC)) {
        message(paste("Reading RDS file:", mtx_QC))
        mtx_QC <- readRDS(mtx_QC)
    }
    
    if (is.character(sg_lib)) {
        message(paste("Reading sgRNA lib file:", sg_lib))
        sg_lib <- read.table(sg_lib, header = T)
    }
    
    PerturbGene <- length(unique(sg_lib$gene)) - length(NTC)
    sg_num <- length(unique(sg_lib$barcode))
    NTC_num <- length(unique(subset(sg_lib, gene %in% NTC)))
    
    #generate config of html

    if (type == "ATAC") {
        if(species == "Hs"){
            if(version == "v75"){
                gene_anno <- "EnsDb.Hsapiens.v75(hg19)"
            }else if(version == "v79"){
                gene_anno <- "EnsDb.Hsapiens.v79(hg38)"
            }else if(version == "v86"){
                gene_anno <- "EnsDb.Hsapiens.v86(hg38)"
            }
        }else if(species == "Mm"){
            if(version == "v75"){
                gene_anno <- "EnsDb.Mmusculus.v75(mm9)"
            }else if(version == "v79"){
                gene_anno <- "EnsDb.Mmusculus.v79(mm10)"
            }
        }
        
        if(gene_type == "Symbol") {
            gene_label <- "Gene Symbol"
        } else if(gene_type == "Ensembl") {
            gene_label <- "Ensembl ID"
        }
            
        if (article != "") {
            if (article_name == "") {
                article_name <- article
            }
            article <- paste('<a href=\"', article, '\" target="view_window">', article_name, '</a>', sep = "")
        }
            
        if (data != "") {
            if (data_name == "") {
                data_name <- data
            }
            data <- paste('<a href=\"', data, '\" target="view_window">', data_name, '</a>', sep = "")
        }
        config <- c("", "", PerturbGene, sg_num, NTC_num, prefix, "single-cell CRISPR ATAC-seq screens", 
                    project, gene_anno, gene_label,
                    p_val_cut, score_cut)
        names(config) <- c("ArticlePath", "GEO", "PerturbGene", "sgNumbers", "NTC", 
                           "ResultsPath", "DataType", "Project", "ReferenceVersion",
                           "GeneVersion", "Score_cut", "pval_cut")
        
    } else if (type == "enhancer") {
        config <- c("", "", PerturbGene, sg_num, NTC_num, prefix, "single-cell CRISPR RNA-seq screens", 
                    project, gene_anno, gene_label,
                    p_val_cut, score_cut)
        names(config) <- c("ArticlePath", "GEO", "PerturbGene", "sgNumbers", "NTC", 
                           "ResultsPath", "DataType", "Project", "ReferenceVersion",
                           "GeneVersion", "Score_cut", "pval_cut")
    } else if (type == "RNA") {
        config <- c("", "", PerturbGene, sg_num, NTC_num, prefix, "single-cell CRISPR RNA-seq screens", 
                    project, gene_anno, gene_label,
                    p_val_cut, score_cut)
        names(config) <- c("ArticlePath", "GEO", "PerturbGene", "sgNumbers", "NTC", 
                           "ResultsPath", "DataType", "Project", "ReferenceVersion",
                           "GeneVersion", "Score_cut", "pval_cut") 
    }

    config <- paste(names(config), config, collapse = "\" , \"", sep = "\" : \"")
    config <- paste("\"Config\" : {\"", config, "\"}", sep = "")
        
        
    #generate single-cell quality of html
    
    if (type == "ATAC") {
        vln <- file.path("img", "ATAC_quality", paste(label, "raw_matrix_quality_vlnplot.png", sep = ""))
        frag <- file.path("img", "ATAC_quality", paste(label, "FragmentsSize.png", sep = ""))
        umap1 <- file.path("img", "ATAC_quality", "UMAP", paste(label, "umap_seurat_clusters.png", sep = ""))
        umap2 <- file.path("img", "ATAC_quality", "UMAP", paste(label, "umap_perturbations.png", sep = ""))
        sg_figure <- c(vln, frag, umap1, umap2)
        names(sg_figure) <- c("vlnplot", "fragments", "UMAP1", "UMAP2")
    } else {
        vln <- file.path("img", "quality", paste(label, "raw_matrix_quality_vlnplot.png", sep = ""))
        umap1 <- file.path("img", "quality", "UMAP", paste(label, "umap_seurat_clusters.png", sep = ""))
        umap2 <- file.path("img", "quality", "UMAP", paste(label, "umap_perturbations.png", sep = ""))
        sg_figure <- c(vln, umap1, umap2)
        names(sg_figure) <- c("vlnplot", "UMAP1", "UMAP2")
    }


    #get information of quality plot
    
    sg_figure <- paste(names(sg_figure), sg_figure, collapse = "\" , \"", sep = "\" : \"")
    sg_figure <- paste("\"figure\" : {\"", sg_figure, "\"}", sep = "")

    #get information of cells and genes
    
    info <- c(project, ncol(mtx), nrow(mtx), ncol(mtx_QC), nrow(mtx_QC))
    names(info) <- c("SampleName", "RawCell", "RawGene", "QCCell", "QCGene")
    info <- paste(names(info), info, collapse = "\" , \"", sep = "\" : \"")
    info <- paste("\"quality\" : {\"", paste(info, sg_figure, sep = "\" , "), "}", sep = "")

    #generate sgRNA information of html
    
    img_prefix <- file.path("img", "sgRNA", "sgRNA_quality_of_gene")
    sg_right <- c(file.path("img", "sgRNA", paste(label, "cell_numbers.png", sep = "")))
    sg_lib_filtered <- subset(sg_lib, cell %in% intersect(sg_lib$cell, colnames(mtx)))
    for(genes in sort(unique(sg_lib_filtered$gene))) {
        sg_right <- c(sg_right, file.path(img_prefix, paste(genes, ".png", sep = "")))
    }
    names(sg_right) <- c("Top10", sort(unique(sg_lib_filtered$gene)))
    sg_right <- paste(names(sg_right), sg_right, collapse = "\" , \"", sep = "\" : \"")
    sg_right <- paste("\"right\" : {\"", sg_right, "\"}", sep = "")
    sg_left <- paste("left\" : \"", file.path("img", "sgRNA", paste(label, "sgRNA_numbers.png", sep = "")))
    sg <- paste("\"sgRNA\" : {\"", paste(sg_left, sg_right, sep = "\" , "), "}", sep = "")

    #generate Mixscape information of html
    
    if (file.exists(file.path(prefix, "img", "perturbation_efficiency", "mixscape", 
                              paste(label, "mixscape_after.png", sep = "")))){
        cluster <- c(file.path("img", "perturbation_efficiency", "mixscape", 
                               paste(label, "umap_seurat_clusters.png", sep = "")), 
                     file.path("img", "perturbation_efficiency", "mixscape", 
                                paste(label, "umap_perturbations.png", sep = "")))
        names(cluster) <- c("UMAP1", "UMAP2")
        cluster <- paste(names(cluster), cluster, collapse = "\" , \"", sep = "\" : \"")
        cluster <- paste("{\"", cluster, "\"}", sep = "")
        umap <- c(file.path("img", "perturbation_efficiency", "mixscape", 
                            paste(label, "mixscape_before.png", sep = "")), 
                  file.path("img", "perturbation_efficiency", "mixscape", 
                            paste(label, "mixscape_after.png", sep = "")),
                  cluster)
        mixscape_heatmap <- c(file.path("img", "perturbation_efficiency", 
                                        paste(label, "before_perturb_ratio.png", sep = "")),
                              file.path("img", "perturbation_efficiency", 
                                        paste(label, "after_perturb_ratio.png", sep = "")))

    } else {
        umap <- c(file.path("img", "perturbation_efficiency", "mixscape", 
                            paste(label, "mixscape_before.png", sep = "")), "", "")
        mixscape_heatmap <- c("", "")
    }

    names(umap) <- c("before", "after", "clustering")
    umap <- paste(names(umap), umap, collapse = "\" , \"", sep = "\" : \"")
    umap <- paste("\"UMAP\" : {\"", umap, "\"}", sep = "")
    names(mixscape_heatmap) <- c("before", "after")
    mixscape_heatmap <- paste(names(mixscape_heatmap), mixscape_heatmap, collapse = "\" , \"", sep = "\" : \"")
    mixscape_heatmap <- paste("\"heatmap\" : {\"", mixscape_heatmap, "\"}", sep = "")

    if (dir.exists(file.path(prefix, "img", "perturbation_efficiency", "mixscape", "KO_percent", ""))) {
        l <- length(list.files(file.path(prefix, "img", "perturbation_efficiency", "mixscape", "KO_percent", "")))
        ko <- c()
        for (i in 1 : l) {
            ko <- c(ko, file.path("img", "perturbation_efficiency", "mixscape", "KO_percent", 
                                  paste(i, ".png", sep = "")))
        }
        names(ko) <- seq(1 : l)
        ko <- paste(names(ko), ko, collapse = "\" , \"", sep = "\" : \"")
        ko <- paste("\"KO\" : {\"", ko, "\"}", sep = "")
    } else {
        ko <- "\"KO\" : \"\""
    }


    mixscape <- paste("\"Mixscape\" : {", paste(umap, mixscape_heatmap, ko, sep = " , "), "}", sep = "")

    #generate scMAGeCK information of html
    
    prefix_volcano <- file.path("img", "perturbation_efficiency", "volcano")
    volcano <- c()
    for (genes in sort(colnames(score))) {
        volcano <- c(volcano, file.path(prefix_volcano, paste(genes, ".png", sep = "")))
    }
    names(volcano) <- sort(colnames(score))
    volcano <- paste(names(volcano), volcano, collapse = "\" , \"", sep = "\" : \"")
    volcano <- paste("\"volcano\" : {\"", volcano, "\"}", sep = "")

    if (file.exists(file.path(prefix, "img", "perturbation_efficiency", 
                              paste(label, "correlation_heatmap.png", sep = "")))) {
        sc <- c(file.path("img", "perturbation_efficiency", 
                          paste(label, "DE_gene_cutoff", score_cut, "_p", p_val_cut, ".png", sep = "")),
                file.path("img", "perturbation_efficiency", paste(label, "correlation_heatmap.png", sep = "")))
    } else {
        sc <- c(file.path("img", "perturbation_efficiency", 
                          paste(label, "DE_gene_cutoff", score_cut, "_p", p_val_cut, ".png", sep = "")), 
                "")
    }
    names(sc) <- c("barplot", "heatmap")
    sc <- paste(names(sc), sc, collapse = "\" , \"", sep = "\" : \"")
    scmageck <- paste("\"scMAGeCK\" : {\"", paste(sc, volcano, sep = "\" , "), "}", sep = "")
        
    #generate DE genes information of html
        
    results <- data.frame()
    GO <- c()
    prefix_GO <- file.path(prefix, "img", "potential_target_gene", "GO")
    go <- list.files(file.path(prefix_GO))
    j <- 0
    k <- 0
    DE <- c()
    for (i in sort(colnames(score))) {
        a_scmageck <- data.frame(target = rownames(score), score = score[, i], p_val = p_val[, i])
        a_scmageck <- subset(a_scmageck, abs(score) > score_cut & p_val < p_val_cut)
        if (paste(i, ".png", sep = "") %in% go) {
            k <- k + 1
            GO <- c(GO, file.path("img", "potential_target_gene", "GO", paste(i, ".png", sep = "")))
            names(GO)[k] <- i
        }
        if (nrow(a_scmageck) == 0) {
            next
        } else {
            j <- j + 1
            de <- data.frame(Perturbations = rep(i, nrow(a_scmageck)), 
                             target = a_scmageck$target, 
                             score = a_scmageck$score,
                             p_val = a_scmageck$p_val)
            results <- rbind(results, de)
            de <- paste("<tr><td>", de$Perturbations, 
                        "</td><td>", de$target, 
                        "</td><td>", de$score, 
                        "</td><td>", de$p_val, 
                        "</td></tr>", sep = "")
            de <- paste(de, collapse = "", sep = "")
            DE <- c(DE, de)
            names(DE)[j] <- i
        }
    }
    if (nrow(results) != 0) {
        all_de <- paste("<tr><td>", results$Perturbations, 
                        "</td><td>", results$target, 
                        "</td><td>", results$score, 
                        "</td><td>", results$p_val, 
                        "</td></tr>", sep = "")
        all_de <- paste(all_de, collapse = "", sep = "")
        DE <- c(all_de, DE)
        names(DE)[1] <- "All_Results"
        DE <- paste(names(DE), DE, collapse = "\" , \"", sep = "\" : \"")
        DE <- paste("\"table\" : {\"", DE, "\"}", sep = "")
    } else {
        DE <- paste("\"All_Results\" : ", "\"\"", sep = "")
        DE <- paste("\"table\" : {", DE, "}", sep = "")
    }
        

    if (length(GO) != 0) {
        GO <- paste(names(GO), GO, collapse = "\" , \"", sep = "\" : \"")
        GO <- paste("\"GO\" : {\"", GO, "\"}", sep = "")
    } else {
        GO <- "\"GO\" : \"\""
    }

    DEgenes <- paste("\"DEgenes\" : {", paste(DE, GO, sep = " , "), "}", sep = "")
        
    #generate enhancer results from output files
        
    if (is.null(enhancer)) {
        if (type == "enhancer") {
            en <- list.files(file.path(prefix, "img", "enhancer_function", "enhancer_gene_expression", ""))
            all_enhancer <- c()
            for (perturb in en) {
                enhancer <- c()
                
                #generate enhancer directory of html
                
                enhancer_prefix <- file.path(file.path(prefix, "img", "enhancer_function", 
                                                       "enhancer_gene_expression", ""), perturb)

                l <- length(list.files(file.path(enhancer_prefix)))
                for(i in 1 : l){
                    enhancer <- c(enhancer, paste(file.path("img", "enhancer_function", 
                                                            "enhancer_gene_expression", perturb, i), ".png", sep = ""))
                    names(enhancer)[i] <- i
                }
                
                enhancer <- paste(names(enhancer), enhancer, collapse = "\" , \"", sep = "\" : \"")
                enhancer <- paste("\"", perturb, "\" : {\"", enhancer, "\"}", sep = "")
                all_enhancer <- c(all_enhancer, enhancer)
                
            }
            all_enhancer <- paste("", all_enhancer, collapse = ", ", sep = "")
            all_enhancer <- paste("\"enhancer\" : {", all_enhancer, "}", sep = "")
            enhancer <- all_enhancer
        }
    }
        
    #generate cicero results from output files
    
    if (is.null(cicero)) {
        if (type == "enhancer") {
            ci <- list.files(file.path(prefix, "img", "enhancer_function", "cicero", ""))
            ci <- gsub(".png", "", ci)
            k <- 0
            cicero <- c()
            for(j in ci) {    
                
                #get html config
                
                k <- k + 1
                cicero <- c(cicero, file.path("img", "enhancer_function", "cicero", paste(j, ".png", sep = "")))
                names(cicero)[k] <- j
            }
            cicero <- paste(names(cicero), cicero, collapse = "\" , \"", sep = "\" : \"")
            cicero <- paste("\"cicero\" : {\"", cicero, "\"}", sep = "")
        }
    }
        
    #generate results information of html
        
    if (type == "enhancer") {
        results <- paste(info, sg, mixscape, scmageck, DEgenes, enhancer, cicero, sep = " , ")
    } else if (type == "ATAC") {
        results <- paste(info, sg, mixscape, scmageck, DEgenes, DA, cicero, sep = " , ")
    } else if (type == "RNA") {
        results <- paste(info, sg, mixscape, scmageck, DEgenes, sep = " , ")
    }
        
    #generate final config value
       
    results <- paste("\"Results\" : {", results, "}", sep = "")
    final <- paste(config, results, sep = " , ")
    final <- paste("\nvar SCREEN_result = {", final, "};\n", sep = "")
        
    return(final)
        
}


#' @import xml2

html_output <- function(html_dir, config, replace = 3, prefix = ".", label = "") {
    html <- xml2::read_html(html_dir)
    xml2::xml_set_text(xml2::xml_child(html, replace), config)
    xml2::write_html(html, file.path(prefix, paste(label, "summary.html", sep = "")))
}