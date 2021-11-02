#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import ggplot2
#' @export

GOenrichment <- function(score_dir, pval_dir, p_val_cut = 0.05, score_cut = 0.5, DE_gene_to_use = "all",
                         database = "org.Hs.eg.db", gene_type = "Symbol", showCategory = 10, 
                         prefix = "./", label = ""){
    
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
    
    #get GO results for each perturbation
    
    results <- list()
    j = 0
    dir <- file.path(prefix, "GO_results")
    if (!(dir.exists(dir))) {
        dir.create(path = dir)
    }
    
    img_dir <- file.path(prefix, "img")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }
    
    img_dir <- file.path(img_dir, "GO")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }

    for(gene in colnames(score)){
        
        #get score and p_val for each perturbation
        
        diff_table <- data.frame(score = score[, gene], p_val = p_val[, gene])
        rownames(diff_table) <- rownames(score)
        
        #define up-regulated, down-regulated and all differential expressed gene list
        
        up_genes <- rownames(subset(diff_table, score > score_cut & p_val < p_val_cut))
        down_genes <- rownames(subset(diff_table, score < -score_cut & p_val < p_val_cut))
        DE_genes <- c(up_genes, down_genes)
        
        #prepare gene list for GO enrichment
        
        if(DE_gene_to_use == "all") {
            de <- DE_genes
        } else if(DE_gene_to_use == "up"){
            de <- up_genes
        } else if(DE_gene_to_use == "down"){
            de <- down_genes
        } else{
            stop("DE_gene_to_use must be one of c('all', 'up', 'down')")
        }
        
        #convert gene name
        
        if (gene_type == "Symbol") {
            de_genes <- bitr(de, "SYMBOL", "ENTREZID", database)
        } else if (gene_type == "Ensembl") {
            de_genes <- bitr(de, "ENSEMBL", "ENTREZID", database)
        } else {
            stop("gene_type must be one of c('Symbol', 'Ensembl')")
        }
            
        #GO enrichment
            
        MF_ego <- enrichGO(gene = de_genes[, 2],
                           OrgDb = database,
                           ont = "MF",
                           keyType = "ENTREZID",
                           pAdjustMethod = "BH",
                           minGSSize = 1,
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           readable = TRUE)
        BP_ego <- enrichGO(gene = de_genes[, 2],
                           OrgDb = database,
                           ont = "BP",
                           keyType = "ENTREZID",
                           pAdjustMethod = "BH",
                           minGSSize = 1,
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           readable = TRUE)
        CC_ego <- enrichGO(gene = de_genes[, 2],
                           OrgDb = database,
                           ont = "CC",
                           keyType = "ENTREZID",
                           pAdjustMethod = "BH",
                           minGSSize = 1,
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05,
                           readable = TRUE)
        
        if (!is.numeric(showCategory)) {
            stop("showCategory must be numeric")
        }
        MF_results <- na.omit(as.data.frame(MF_ego))[1:showCategory, ]
        BP_results <- na.omit(as.data.frame(BP_ego))[1:showCategory, ]
        CC_results <- na.omit(as.data.frame(CC_ego))[1:showCategory, ]

        all_results <- as.data.frame(rbind(MF_results, CC_results, BP_results))
        colnames(all_results) <- colnames(BP_results)

        all_results$ONTOLOGY <- factor(rep(c("MF", "CC", "BP"), each = 10), levels = c("BP", "CC", "MF"), ordered = T)

        #plot
        
        g1 <- ggplot(all_results) + geom_bar(aes(x = Description, y = Count, fill = ONTOLOGY), stat = 'identity') +
        labs(x = "GO Terms",y = "Gene Numbers",title = "GO Enrichment Results", fill = "Ontology") + 
        coord_flip() + scale_x_discrete(limits = all_results$Description) + theme_bw() +
        theme(panel.grid = element_blank()) +  scale_fill_manual(values = c("#8DA1CB", "#FD8D62", "#66C3A5")) +
        theme(plot.title = element_text(hjust = 0.5, size = 20), 
              axis.text.y = element_text(size = 12),
              axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 16),
              axis.title.x = element_text(size = 16))
        
        #save results and return
        
        j = j + 1
        results[[j]] <- g1
        names(results)[j] <- gene
        
        pdf(file = file.path(dir, paste(gene, ".pdf", sep = "")), width = 18, height = 6)
        print(g1)
        dev.off()
        
        png(file.path(img_dir, paste(gene, ".png", sep = "")), 
            width = 1800 * 3, height = 600 * 3, res = 72 * 3)
        print(g1)
        dev.off()
    }
    return(results)
}