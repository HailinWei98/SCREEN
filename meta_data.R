#library packages
library(Seurat)
library(stringi)
library(stringr)
library(plyr)

#get parameter
args<- commandArgs(T)
mtx_dir<- args[1]
sg_dir<- args[2]
species<- args[3]

#read files
sg_lib<- read.table(sg_dir,header=T)
mtx<- readRDS(mtx_dir)

#remove cells in sgRNA library that are not included in matrix
sg_lib_filtered<- subset(sg_lib,cell %in% intersect(sg_lib$cell,colnames(mtx)))

#label each cells
label<- rep("blank",times=ncol(mtx))
sg_num<- rep(0,times=ncol(mtx))
names(label)<- colnames(mtx)
names(sg_num)<- colnames(mtx)
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
mtx<- AddMetaData(mtx,label,col.name = "perturbations")
mtx<- AddMetaData(mtx,sg_num,col.name = "sgRNA_num")

#calculate percent.mt
if(species=="Hs"){
    mtx[["percent.mt"]] <- PercentageFeatureSet(mtx, pattern = "^MT-")
}else{
    mtx[["percent.mt"]] <- PercentageFeatureSet(mtx, pattern = "^mt-")
}

saveRDS(mtx,file="perturb.rds")