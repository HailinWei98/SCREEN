#' function definitions ##### Run all function with one command
#' @export

SCREEN<- function(sg_dir, mtx_dir, species = "Hs",
                  prefix = "./", label = "", gene_frac = 0.01,
                  nFeature = c(200, 5000), nCount = 1000,
                  mt = 10, blank_NTC = FALSE, lambda = 0.01, permutation = NULL,
                  p_val_cut = 0.05, score_cut = 0.5,
                  ylimit = c(-600, 600, 200), project = "perturb",
                  NTC = "NTC", replicate = 1, select_gene = NULL){
  if (is.character(mtx_dir)) {
    message(paste("Reading RDS file:", mtx_dir))
    mtx <- readRDS(mtx_dir)
  } else {
    mtx <- mtx_dir
  }
  if (is.character(sg_dir)) {
    message(paste("Reading sgRNA lib file:", sg_dir))
    sg_lib <- read.table(sg_dir,header=T)
  } else {
    sg_lib <- sg_dir
  }
  mtx<- Add_meta_data(sg_lib, mtx, species, prefix, label)
  mtx_QC<- scQC(mtx, prefix, species, gene_frac, nFeature, nCount, mt, blank_NTC)
  p<- sgRNA_quality_plot(sg_lib, mtx,
                         label = paste(label, "before_QC", sep = ""), prefix)
  q<- sgRNA_quality_plot(sg_lib, mtx_QC,
                         label = paste(label, "after_QC", sep = ""), prefix)
  mixscape<- IntegratedMixscape(sg_lib, mtx_QC, NTC, replicate, prefix)
  mtx_QC<- normalize_scale(mtx_QC, label, prefix)
  mtx_QC<- GetAssayData(object = mtx_QC, slot = "scale.data")
  results<- improved_scmageck_lr(sg_lib, mtx_QC, NTC, select_gene,
                                 LABEL = paste(lael, "improved",sep = ""),
                                 permutation, prefix, lambda)
  y<- DE_gene_plot(results[1], results[2], project,
                   prefix, p_val_cut, score_cut, label,
                   ylimit)
  return(results)
}
