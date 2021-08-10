#' function definitions ##### Normalize and scale data.
#' @export

normalize_scale<- function(mtx_dir, label = "", prefix = "./"){
  #read file
  if (is.character(mtx_dir)) {
    message(paste("Reading RDS file:", mtx_dir))
    perturb_QC <- readRDS(mtx_dir)
  } else {
    perturb_QC <- mtx_dir
  }

  #normalize and scale on the data
  perturb_QC <- NormalizeData(
    object = perturb_QC,
    normalization.method = "LogNormalize",
    scale.factor = 10000)
  perturb_QC <- FindVariableFeatures(perturb_QC, selection.method = "vst", nfeatures = 2000)
  perturb_QC <- ScaleData(perturb_QC, features = rownames(perturb_QC),vars.to.regress = c("nCount_RNA","percent.mt"))

  saveRDS(perturb_QC, file=file.path(prefix, paste(label, "perturb_QC.rds", sep = "")))
  saveRDS(GetAssayData(object = perturb_QC, slot = "scale.data"), file=file.path(prefix, paste(label, "scale_data.rds", sep = "")))
  return(perturb_QC)
}
