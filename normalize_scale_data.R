library(Seurat)

#get parameter
args<- commandArgs(T)
rds_dir<- args[1]

#read file
perturb_QC<- readRDS(rds_dir)

#normalize and scale on the data
perturb_QC <- NormalizeData(
  object = perturb_QC,
  normalization.method = "LogNormalize",
  scale.factor = 10000)
perturb_QC <- FindVariableFeatures(perturb_QC, selection.method = "vst", nfeatures = 2000)
perturb_QC <- ScaleData(perturb_QC, features = rownames(perturb_QC),vars.to.regress = c("nCount_RNA","percent.mt"))

saveRDS(perturb_QC,file="perturb_QC.rds")