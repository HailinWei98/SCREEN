library(Seurat)
library(ggplot2)
library(plyr)

#get parameter
args<- commandArgs(T)
rds_dir<- args[1]
sg_dir<- args[2]

#read files
perturb_QC<- readRDS(rds_dir)
sg_lib<- read.table(sg_dir,header=T)

#do pca on the data
perturb_QC <- RunPCA(
  object = perturb_QC,
  features = perturb_QC@assays[["RNA"]]@var.features)

sg_lib_filtered<- subset(sg_lib,cell %in% intersect(sg_lib$cell,colnames(perturb_QC)))

#label each cells
label<- rep("blank",times=ncol(perturb_QC))
sg_num<- rep(0,times=ncol(perturb_QC))
names(label)<- colnames(perturb_QC)
names(sg_num)<- colnames(perturb_QC)
sg_in_cell<- data.frame(count(sg_lib_filtered$cell))

#find unique label and multiple label
for(i in 1:nrow(sg_in_cell)){
    x<- sg_in_cell[i,]
    if(x[,2]==1){
        label[x[,1]]<- sg_lib_filtered[which(sg_lib_filtered$cell==x[,1]),3]}else{
        label[x[,1]]<- "multiple"}
}

for(i in 1:nrow(sg_in_cell)){
    if(sg_in_cell[i,1] %in% names(sg_num)){
        sg_num[sg_in_cell[i,1]]<- sg_in_cell[i,2]
    }
}
#add the label information to Seurat object
label<- as.factor(label)
mtx<- AddMetaData(perturb_QC,label,col.name = "perturbations")
mtx<- AddMetaData(perturb_QC,sg_num,col.name = "sgRNA_num")

#visualize sgRNA label of each cell in UMAP
perturb_QC <- FindNeighbors(perturb_QC, dims = 1:20)
perturb_QC <- FindClusters(perturb_QC)
perturb_QC<- RunUMAP(perturb_QC,dims = 1:20)

pdf("umap_perturbations.pdf")
DimPlot(perturb_QC, reduction="umap", group.by = "perturbations")
dev.off()
pdf("umap_seurat_clusters.pdf")
DimPlot(perturb_QC, reduction="umap")
dev.off()
pdf("umap_NTC_vs_blank.pdf")
DimPlot(perturb_QC, reduction="umap", group.by = "perturbations", cells = colnames(subset(perturb_QC, perturbations %in% c("NTC","blank"))))
dev.off()

saveRDS(perturb_QC,file="perturb_QC.rds")