#' @import ggplot2
#' @import plyr
#' @export

fragmentsSize <- function(mtx_dir, fragments, CBCindex = 4, startIndex = 2, 
                          endIndex = 3, maxSize = 1000, prefix = "./", label = ""){
    
    #get fragments information
    
    if (is.character(fragments)) {
        message(paste("Reading fragments file:", fragments))
        fragments <- read.table(gzfile(fragments), header = F)
    } else {
        fragments <- fragments
    }
    
    #get peak matrix
    
    if (is.character(mtx_dir)) {
        message(paste("Reading RDS file:", mtx_dir))
        peak <- readRDS(mtx_dir)
    } else {
        peak <- mtx_dir
    }

    #filter cells not in peak matrix
    
    fragments <- fragments[(fragments[, CBCindex] %in% colnames(peak)), ]
    
    #calculate fragments size
    
    fragments$size <- fragments[, 3] - fragments[, 2]
    
    #plot
    
    fragments_count <- plyr::count(fragments$size)
    fragments_count <- subset(fragments_count, x <= maxSize)
    g1 <- ggplot(data = fragments_count, mapping = aes(x = x, y = freq)) +
    geom_line(size = 1, color = "red") + theme_classic() +
    labs(x = "Fragments size", y = "Frequency",title = "Fragments Size Distribution") +
    theme(plot.title = element_text(hjust = 0.5, size = 20), 
          axis.text.y = element_text(angle = 90, hjust = 0.5, size = 12),
          axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 16),
          axis.title.x = element_text(size = 16))
    
    #save plot
    
    dir <- file.path(prefix, "img")
    if (!(dir.exists(dir))) {
        dir.create(dir)
    }
    
    pdf(file = file.path(prefix, paste(label, "FragmentsSize.pdf", sep = "")))
    print(g1)
    dev.off()
    
    png(file.path(dir, paste(label, "FragmentsSize.png", sep = "")), 
        width = 600 * 3, height = 600 * 3, res = 72 * 3)
    print(g1)
    dev.off()
    return(g1)
}