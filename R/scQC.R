#' function definitions ##### Common single-cell quality control.
#' This step will only save the VlnPlot before and after QC but not show.
#' Users can do this step by themselves using Seurat.
#' @export

scQC<- function(mtx_dir, prefix = "./", label = "", species = "Hs", gene_frac = 0.01,
                nFeature = c(200, 5000), nCount = 1000, mt = 10, blank_NTC = FALSE){
  #read file
  if (is.character(mtx_dir)) {
    message(paste("Reading RDS file:", mtx_dir))
    perturb <- readRDS(mtx_dir)
  } else {
    perturb <- mtx_dir
  }

    perturb$nFeature_RNA <- perturb[[paste("nFeature_", perturb@active.assay, sep = "")]][, 1]
    perturb$nCount_RNA <- perturb[[paste("nCount_", perturb@active.assay, sep = "")]][, 1]
    
  #QC plot of the single cell matrix
  pdf(file = file.path(prefix, paste(label, "raw_matrix_quality_vlnplot.pdf", sep = "")))
  p1 <- VlnPlot(perturb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
  print(p1)
  dev.off()


  #filter cells with low quality
  if(blank_NTC == TRUE){
    perturb_QC <- subset(perturb,
                        nFeature_RNA <= nFeature[2] &
                          nFeature_RNA >= nFeature[1] &
                          nCount_RNA >= nCount &
                          percent.mt <= mt)
  }else{
    perturb_QC <- subset(perturb,
                       nFeature_RNA <= nFeature[2] &
                         nFeature_RNA >= nFeature[1] &
                         nCount_RNA >= nCount &
                         percent.mt <= mt &
                         perturbations != 'blank')
  }
  perturb_QC <- CreateSeuratObject(counts = GetAssayData(object = perturb_QC, slot = "counts"),
                                  min.cells = gene_frac * ncol(perturb_QC), project = perturb@project.name)
    if("replicate" %in% colnames(perturb@meta.data)){
        perturb_QC$replicate <- replicate
    }
    
    if("perturbations" %in% colnames(perturb@meta.data)){
        perturb_QC$perturbations <- perturb$perturbations
    }
  #calculate percent.mt
    
  if(species == "Hs"){
    perturb_QC[["percent.mt"]] <- PercentageFeatureSet(perturb_QC, pattern = "^MT-")
  }else{
    perturb_QC[["percent.mt"]] <- PercentageFeatureSet(perturb_QC, pattern = "^mt-")
  }
  pdf(file = file.path(prefix, paste(label, "QC_matrix_quality_vlnplot.pdf", sep = "")))
  p2 <- VlnPlot(perturb_QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
  print(p2)
  dev.off()

  #saveRDS(perturb_QC,file=file.path(prefix, "perturb_QC.rds"))
  return(perturb_QC)
}
