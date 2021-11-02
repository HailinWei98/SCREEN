#' @import Signac
#' @import Seurat
#' @export

umap <- function(mtx_dir, nfeature = 2000, dims = 1:20, data_type = "RNA",  
                 assays = "RNA", prefix = "./", label = ""){
    
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
    
    perturb_QC <- FindVariableFeatures(perturb_QC, selection.method = "vst", nfeatures = nfeature)
    perturb_QC <- RunPCA(object = perturb_QC, features = perturb_QC@assays[[assays]]@var.features)
    perturb_QC <- FindNeighbors(perturb_QC, dims = dims)
    perturb_QC <- FindClusters(perturb_QC)
    perturb_QC <- RunUMAP(perturb_QC, dims = dims)
    
    pdf(file = file.path(prefix, paste(label, "umap_perturbations.pdf", sep = "")))
    p1 <- DimPlot(perturb_QC, reduction = "umap", group.by = "perturbations")
    print(p1)
    dev.off()
    pdf(file = file.path(prefix, paste(label, "umap_seurat_clusters.pdf", sep = "")))
    p2 <- DimPlot(perturb_QC, reduction = "umap")
    print(p2)
    dev.off()
    
    dir <- file.path(prefix, "img")
    if (!(dir.exists(dir))) {
        dir.create(dir)
    }
    dir <- file.path(dir, "UMAP")
    if (!(dir.exists(dir))) {
        dir.create(dir)
    }
    
    png(file.path(dir, paste(label, "umap_perturbations.png", sep = "")), 
        width = 500, height = 500)
    print(p1)
    dev.off()
    png(file.path(dir, paste(label, "umap_seurat_clusters.png", sep = "")), 
        width = 500, height = 500)
    print(p2)
    dev.off()
    return(p1/p2)
}