library(Seurat)

#get parameter
args<- commandArgs(T)
rds_dir<- args[1]

#read file
perturb<- readRDS(rds_dir)

#QC plot of the single cell matrix
pdf("raw_matrix_quality_vlnplot.pdf")
VlnPlot(perturb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()

#filter cells with low quality
perturb_QC<- subset(perturb, 
                     nFeature_RNA <= 5000 & 
                     nFeature_RNA >= 200 & 
                     nCount_RNA >= 1000 &
                     percent.mt <= 10)
pdf("QC_matrix_quality_vlnplot.pdf")
VlnPlot(perturb_QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()

saveRDS(perturb_QC,file="perturb_QC.rds")
saveRDS(perturb,file="perturb.rds")