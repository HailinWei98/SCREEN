#library packages

library(plyr)
library(dplyr)
library(Seurat)
library(stringr)

#get parameter
args<- commandArgs(T)
mtx_dir<- args[1]
project<- args[2]
species<- args[3]

#read file
crop.data<- read.table(gzfile(mtx_dir),row.names=1,header=T)

#process dataset and get sgRNA information
sg_info<- as.vector(t(crop.data[2,]))
sg_names<- as.vector(t(crop.data[3,]))
mtx<- data.frame(crop.data[6:nrow(crop.data),])
gene<- as.vector(t(crop.data[4,]))
colnames(mtx)<- sg_info
gene<- gsub("CTRL","NTC",gene)

#create sgRNA lib
rename<- data.frame(unique(cbind(sg_names,gene)))
rownames(rename)<- rename$sg_names
for(i in unique(gene)){
    a<- subset(rename,gene==i)
    a$sg_names<- paste(i,which(a$gene==i),sep="_sgRNA")
    for(j in rownames(a)){
        sg_names<- str_replace(sg_names,j,a[j,1])
    }
}
sg_lib<- data.frame(t(rbind(sg_info,sg_names,gene)))
names(sg_lib)<- c("cell","barcode","gene")

#convert matrix to Seurat object
perturb<- CreateSeuratObject(counts = mtx, project = project)

#label each cells
sg_num<- rep(1,times=ncol(mtx))
names(sg_num)<- colnames(mtx)

#add the label information to Seurat object
sg_names<- as.factor(sg_names)
perturb<- AddMetaData(perturb,sg_names,col.name = "perturbations")
perturb<- AddMetaData(perturb,sg_num,col.name = "sgRNA_num")
                     
#QC plot of the single cell matrix
if(species=="Hs"){
    perturb[["percent.mt"]] <- PercentageFeatureSet(perturb, pattern = "^MT-")
}else{
    perturb[["percent.mt"]] <- PercentageFeatureSet(perturb, pattern = "^mt-")
}
pdf("raw_matrix_quality_vlnplot.pdf")
VlnPlot(perturb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
saveRDS(perturb,file="perturb.rds")
write.table(sg_lib,file="sg_lib_all.txt",row.names=FALSE)
write.table(sg_lib,file = "sg_lib.txt",row.names = FALSE)
