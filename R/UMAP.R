#' @import Signac
#' @import Seurat
#' @export

umap <- function(mtx_dir, nfeature = 2000, dims = 1:20, assays = "RNA", algorithm = 1,
                 prefix = "./", label = "", label.cut = 20, resolution = 0.8){
    
    #read file
    
    if (is.character(mtx_dir)) {
        message(paste("Reading RDS file:", mtx_dir))
        perturb_QC <- readRDS(mtx_dir)
    } else {
        perturb_QC <- mtx_dir
    }

    #perturbation information
    
    if(!("perturbations" %in% colnames(perturb_QC@meta.data))){
        stop("Please add perturbation information to the SeuratObject first.")
    }
    
    #visualize sgRNA label of each cell in UMAP
    
    perturb_QC <- FindVariableFeatures(perturb_QC, assay = assays, selection.method = "vst", nfeatures = nfeature)
    perturb_QC <- RunPCA(object = perturb_QC, assay = assays, features = perturb_QC@assays[[assays]]@var.features)
    perturb_QC <- FindNeighbors(perturb_QC, assay = assays, dims = dims)
    perturb_QC <- FindClusters(perturb_QC, resolution = resolution, algorithm = algorithm)
    perturb_QC <- RunUMAP(perturb_QC, dims = dims, assay = assays)
    
    custom_theme <- theme(
        plot.title = element_text(size = 25, hjust = 0.5),
        legend.key.size = unit(0.7, "cm"),
        legend.text = element_text(size = 14), 
        axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(hjust = 0.5, size = 16), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20))
    
    dir <- file.path(prefix, "pdf")
    if (!(dir.exists(dir))) {
        dir.create(path = dir)
    }
    
    dir <- file.path(dir, "quality")
    if (!(dir.exists(dir))) {
        dir.create(path = dir)
    }
    
    dir <- file.path(dir, "UMAP")
    if (!(dir.exists(dir))) {
        dir.create(path = dir)
    }
    
    img_dir <- file.path(prefix, "img")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }
    
    img_dir <- file.path(img_dir, "quality")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }
    
    img_dir <- file.path(img_dir, "UMAP")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }
    
    p1 <- DimPlot(perturb_QC, reduction = "umap", group.by = "perturbations") + 
    ggtitle("Perturbations") +
    ylab("UMAP 2") +
    xlab("UMAP 1") +
    custom_theme
    if(length(unique(perturb_QC$perturbations)) > label.cut) {
        p1 <- p1 + NoLegend()
    }
    
    p2 <- DimPlot(perturb_QC, reduction = "umap") +
    ggtitle("Seurat clusters") +
    ylab("UMAP 2") +
    xlab("UMAP 1") +
    custom_theme
    
    pdf(file = file.path(dir, paste(label, "umap_perturbations.pdf", sep = "")))
    print(p1)
    dev.off()
    
    
    pdf(file = file.path(prefix, paste(label, "umap_seurat_clusters.pdf", sep = "")))
    print(p2)
    dev.off()
    
    png(file.path(img_dir, paste(label, "umap_perturbations.png", sep = "")), 
        width = 600, height = 600)
    print(p1)
    dev.off()
    png(file.path(img_dir, paste(label, "umap_seurat_clusters.png", sep = "")), 
        width = 600, height = 600)
    print(p2)
    dev.off()
    
    print(p1/p2)
    return(perturb_QC)
}

#' @export

ATACumap <- function(mtx_dir, min.cutoff = "q5", n = 20, assays = NULL, 
                     reduction.key = 'LSI_', reduction.name = 'lsi', algorithm = 1,
                     prefix = "./", label = "", label.cut = 20, resolution = 0.8) {
    
    #read file
    
    if (is.character(mtx_dir)) {
        message(paste("Reading RDS file:", mtx_dir))
        peak <- readRDS(mtx_dir)
    } else {
        peak <- mtx_dir
    }
    
    peak<- RunTFIDF(peak)
    peak<- FindTopFeatures(peak, assay = assays, min.cutoff = "q0")
    peak <- RunSVD(
        object = peak,
        assay = assays,
        reduction.key = reduction.key,
        reduction.name = reduction.name
    )
    
    a <- DepthCor(peak, assay = assays, reduction = reduction.name, n = n)
    
    peak <- FindNeighbors(object = peak, reduction = reduction.name, assay = assays, 
                          dims = subset(a$data, abs(counts) < 0.5)$Component)
    peak <- FindClusters(peak, resolution = resolution, algorithm = algorithm)
    
    peak <- RunUMAP(peak, dims = subset(a$data, abs(counts) < 0.5)$Component, assay = assays)

        custom_theme <- theme(
        plot.title = element_text(size = 25, hjust = 0.5),
        legend.key.size = unit(0.7, "cm"),
        legend.text = element_text(size = 14), 
        axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(hjust = 0.5, size = 16), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20))
    
    dir <- file.path(prefix, "pdf")
    if (!(dir.exists(dir))) {
        dir.create(path = dir)
    }
    
    dir <- file.path(dir, "ATAC_quality")
    if (!(dir.exists(dir))) {
        dir.create(path = dir)
    }
    
    dir <- file.path(dir, "UMAP")
    if (!(dir.exists(dir))) {
        dir.create(path = dir)
    }
    
    img_dir <- file.path(prefix, "img")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }
    
    img_dir <- file.path(img_dir, "ATAC_quality")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }
    
    img_dir <- file.path(img_dir, "UMAP")
    if (!(dir.exists(img_dir))) {
        dir.create(img_dir)
    }
    
    p1 <- DimPlot(perturb_QC, reduction = "umap", group.by = "perturbations") + 
    ggtitle("Perturbations") +
    ylab("UMAP 2") +
    xlab("UMAP 1") +
    custom_theme
    if(length(unique(perturb_QC$perturbations)) > label.cut) {
        p1 <- p1 + NoLegend()
    }
    
    p2 <- DimPlot(perturb_QC, reduction = "umap") +
    ggtitle("Seurat clusters") +
    ylab("UMAP 2") +
    xlab("UMAP 1") +
    custom_theme
    
    pdf(file = file.path(dir, paste(label, "umap_perturbations.pdf", sep = "")))
    print(p1)
    dev.off()
    
    pdf(file = file.path(prefix, paste(label, "umap_seurat_clusters.pdf", sep = "")))
    print(p2)
    dev.off()
    
    png(file.path(img_dir, paste(label, "umap_perturbations.png", sep = "")), 
        width = 600, height = 600)
    print(p1)
    dev.off()
    png(file.path(img_dir, paste(label, "umap_seurat_clusters.png", sep = "")), 
        width = 600, height = 600)
    print(p2)
    dev.off()
    
    print(p1/p2)
    
    return(peak)
}
