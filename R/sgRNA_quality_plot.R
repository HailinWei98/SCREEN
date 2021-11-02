#' function definitions ##### Plot for sgRNA information.
#' @export

sgRNA_quality_plot<- function(sg_dir, mtx_dir, LABEL = "", prefix = "./"){
    
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

    #remove cells in sgRNA library that are not included in matrix
    
    sg_lib_filtered <- subset(sg_lib, cell %in% intersect(sg_lib$cell, colnames(mtx)))

    #count sgRNA in each cell;count cells of each sgRNA
    
    sg_count <- plyr::count(subset(sg_lib_filtered, cell %in% colnames(mtx))$barcode)
    colnames(sg_count) <- c("sgRNA","freq")
    sg_count <- sg_count[order(-sg_count$freq), ]
    sg_count$order <- seq(1, nrow(sg_count))

    #prepare for plot each gene
    
    sg_gene <- unique(sg_lib[ ,c("barcode", "gene")])
    rownames(sg_gene) <- sg_gene$barcode
    sg_count$gene <- sg_gene[sg_count$sgRNA, "gene"]

    g1 <- ggplot(data = sg_count, mapping = aes(x = order, y = freq)) +
    geom_line(size = 1) + theme_classic() +
    geom_point(data = head(sg_count, 10),
               mapping = aes(x = order, y = freq, color = sgRNA), size = 5) +
    labs(x = "sgRNA",y = "cell numbers",title = "cell numbers of sgRNA") +
    theme(plot.title = element_text(hjust = 0.5, size = 20), 
          legend.text = element_text(size = 12),
          text = element_text(hjust = 0.5, face = "bold"),
          axis.ticks.x = element_blank(), 
          axis.text.x = element_blank(), 
          legend.title = element_text(size = 18),
          axis.title.x = element_text(size = 16), 
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(angle = 90, hjust = 0.5, size = 12)) +
    scale_color_discrete(name = "top 10 sgRNA", breaks = head(sg_count, 10)$sgRNA) +
    guides(colour = guide_legend(title.hjust = 0.5))

    sg_num_count <- plyr::count(mtx[["sgRNA_num"]])
    if (max(mtx$sgRNA_num) > 10) {
        over_10 <- subset(sg_num_count, sgRNA_num > 10)
        sg_num_count <- subset(sg_num_count, sgRNA_num <= 10)
        over_10 <- sum(over_10$freq)
        sg_num_count <- rbind(sg_num_count, c(">10", over_10))
    }
    colnames(sg_num_count) <- c("sgRNA_num", "freq")
    sg_num_count$freq <- as.numeric(sg_num_count$freq)
    sg_num_count$sgRNA_num <- factor(sg_num_count$sgRNA_num, levels = sg_num_count$sgRNA_num, ordered = T)

    g2 <- ggplot(sg_num_count,mapping = aes(x = sgRNA_num, y = freq)) +
    geom_bar(stat = "identity",color = "#e9ecef", fill = "#69b3a2") + theme_classic() + 
    labs(x = "sgRNA numbers", y = "cell numbers", title = "sgRNA numbers in each cell") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          text = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
          axis.title.x = element_text(size = 16), 
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(angle = 90, hjust = 0.5, size = 12)) +
    geom_text(aes(label = freq), stat = "identity", vjust = -0.5)

    #save plot
    
    pdf(file = file.path(prefix, paste(LABEL, "sgRNA_quality.pdf", sep = "")))
    print(g1)
    print(g2)
    dev.off()
    
    img_dir <- file.path(prefix, "img")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }
    
    img_dir <- file.path(img_dir, "sgRNA")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }

    
    png(file.path(img_dir, paste(LABEL, "cell_numbers.png", sep = "")), 
        width = 500, height = 500)
    print(g1)
    dev.off()
    
    png(file.path(img_dir, paste(LABEL, "sgRNA_numbers.png", sep = "")), 
        width = 500, height = 500)
    print(g2)
    dev.off()

    dir2 <- file.path(prefix, "sgRNA_quality_of_gene")
    dir <- file.path(img_dir, "sgRNA_quality_of_gene")
    if (!(dir.exists(dir))) {
        dir.create(path = dir)
    }
    if (!(dir.exists(dir2))) {
        dir.create(path = dir2)
    }
    
    for(select_gene in unique(sg_count$gene)){
        pdf(file = file.path(dir2, paste(select_gene, ".pdf", sep = "")))
        p1 <- ggplot(data = sg_count, mapping = aes(x = order, y = freq)) +
        geom_line(size = 1) + theme_classic() +
        geom_point(data = subset(sg_count, gene == select_gene),
                   mapping = aes(x = order, y = freq, color = sgRNA), size = 5)+
        labs(x = "sgRNA",y = "cell numbers",title = "cell numbers of sgRNA") +
        theme(plot.title = element_text(hjust = 0.5, size = 20), 
              legend.text = element_text(size = 12),
              text = element_text(hjust = 0.5, face = "bold"),
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank(), 
              legend.title = element_text(size = 18),
              axis.title.x = element_text(size = 16), 
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(angle = 90, hjust = 0.5, size = 12)) +
        scale_color_discrete(name = paste("sgRNA of", select_gene, sep = " "),
                             breaks = head(subset(sg_count, gene == select_gene), 10)$sgRNA) +
        guides(colour = guide_legend(title.hjust = 0.5))
        print(p1)
        dev.off()
        
        png(file.path(dir, paste(select_gene, ".png", sep = "")), 
            width = 500, height = 500)
        print(p1)
        dev.off()
    }
    return(list(g1, g2))
}
