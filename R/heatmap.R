#' @import pheatmap
#' @import psych
#' @export

heatmap <- function(score_dir, pval_dir, score_cut = 0.5, p_val_cut = 0.05, 
                    num = 0, prefix= "./", label = "", project = "perturb") {
    
    #get score and p_val
    
    if(is.character(score_dir)){
        score <- read.table(score_dir, header = T, row.names = 1)
    }else{
        score <- score_dir
    }
    if(is.character(pval_dir)){
        p_val <- read.table(pval_dir, header = T, row.names = 1)
    }else{
        p_val <- pval_dir
    }
    
    #filter data frame
    
    remove_neg <- score[, -which(colnames(score) == "NegCtrl")]
    p_val_filter <- apply(p_val[, colnames(remove_neg)], 1, function(x) sum(x < p_val_cut))
    score_filter <- apply(remove_neg, 1, function(x) sum(abs(x) > score_cut))
    gene_count <- p_val_filter[p_val_filter > num & score_filter > num]
    filtered <- remove_neg[names(gene_count), ]
    
    #calculate correlation
                          
    sg_cor <- corr.test(filtered)
    
    bk <- c(seq(-1, -0.1, by = 0.01), seq(0, 1, by = 0.01))
                          
    dir <- file.path(prefix, "pdf")
    if (!(dir.exists(dir))) {
        dir.create(dir)
    }
    
    dir <- file.path(dir, "perturbation_efficiency")
    if (!(dir.exists(dir))) {
        dir.create(dir)
    }
    
    img_dir <- file.path(prefix, "img")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }
    
    img_dir <- file.path(img_dir, "perturbation_efficiency")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }
                          
    file <- file.path(dir, paste(label, "correlation_heatmap.pdf", sep = ""))
    pheatmap(sg_cor$r, treeheight_col = F, treeheight_row = F, 
             display_numbers = matrix(ifelse(sg_cor$p <= 0.01, "**", ifelse(sg_cor$p <= 0.05 , "*", " ")),
                                      nrow(sg_cor$p)), 
             main = "Perturbations Correlation",
             c(colorRampPalette(colors = c("#00BFC4", "white"))(length(bk)/2),
               colorRampPalette(colors = c("white", "#F8766D"))(length(bk)/2)),
             legend_breaks = seq(-1, 1, 0.5), annotation_row = names(sg_cor$r),
             breaks = bk, show_rownames = T, show_colnames = T, fontsize = 15, 
             cellwidth = 15, cellheight = 15, filename = file)
    
    img_file <- file.path(img_dir, paste(label, "correlation_heatmap.png", sep = ""))
    pheatmap(sg_cor$r, treeheight_col = F, treeheight_row = F, 
             display_numbers = matrix(ifelse(sg_cor$p <= 0.01, "**", ifelse(sg_cor$p <= 0.05 , "*", " ")), 
                                      nrow(sg_cor$p)), 
             main = "Perturbations Correlation",
             c(colorRampPalette(colors = c("#00BFC4", "white"))(length(bk)/2),
               colorRampPalette(colors = c("white", "#F8766D"))(length(bk)/2)),
             legend_breaks = seq(-1, 1, 0.5), annotation_row = names(sg_cor$r),
             breaks = bk, show_rownames = T, show_colnames = T, fontsize = 15, 
             cellwidth = 15, cellheight = 15, filename = img_file)

    p <- pheatmap(sg_cor$r, treeheight_col = F, treeheight_row = F, 
                  display_numbers = matrix(ifelse(sg_cor$p <= 0.01, "**", ifelse(sg_cor$p <= 0.05 , "*", " ")), 
                                           nrow(sg_cor$p)), 
                  main = "Perturbations Correlation",
                  c(colorRampPalette(colors = c("#00BFC4", "white"))(length(bk)/2),
                    colorRampPalette(colors = c("white", "#F8766D"))(length(bk)/2)),
                  legend_breaks = seq(-1, 1, 0.5), annotation_row = names(sg_cor$r),
                  breaks = bk, show_rownames = T, show_colnames = T, fontsize = 15, 
                  cellwidth = 15, cellheight = 15)
    return(p)
}