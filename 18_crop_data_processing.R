#library packages
library(Seurat)
library(ggplot2)
library(stringr)
library(plyr)

#get parameter
args<- commandArgs(T)
mtx_dir<- args[1]
sg_dir<- args[2]
project<- args[3]
species<- args[4]

#read files
mtx<- Read10X(data.dir = mtx_dir)
sg_dict<- read.csv(gzfile(sg_dir))

#create sg_lib
gene_symbol<- t(data.frame(lapply(strsplit(sg_dict[,2],"[.]"),function(x) x[3])))
sg_lib<- data.frame(cbind(sg_dict[,1:2],gene_symbol))
names(sg_lib)<- c("cell","barcode","gene")
sg_lib$cell<- paste0(sg_lib$cell,"-1")
sg_lib$gene<- gsub("NonTarget","NTC",sg_lib$gene)
                                  
#rename sgRNAs
rename<- data.frame(unique(sg_lib[,2:3]))
rownames(rename)<- rename$barcode
for(i in unique(sg_lib$gene)){
    a<- subset(rename,gene==i)
    a$barcode<- paste(i,which(a$gene==i),sep="_sgRNA")
    for(j in rownames(a)){
        sg_lib$barcode<- str_replace(sg_lib$barcode,j,a[j,1])
    }
}

#remove cells in sgRNA library that are not included in matrix
sg_lib_filtered<- subset(sg_lib,cell %in% intersect(sg_lib$cell,colnames(mtx)))
                                  
#convert matrix to Seurat object
perturb<- CreateSeuratObject(counts = mtx, project = project)
                    
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
        label[x[,1]]<- sg_lib_filtered[which(sg_lib_filtered$cell==x[,1]),3]}
    else{
        label[x[,1]]<- "multiple"}
}

for(i in 1:nrow(sg_in_cell)){
    if(sg_in_cell[i,1] %in% names(sg_num)){
        sg_num[sg_in_cell[i,1]]<- sg_in_cell[i,2]
    }
}

#add the label information to Seurat object
label<- as.factor(label)
perturb<- AddMetaData(perturb,label,col.name = "perturbations")
perturb<- AddMetaData(perturb,sg_num,col.name = "sgRNA_num")

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
sg_lib<- subset(sg_lib,cell %in% subset(sg_in_cell,freq==1)[,1])
write.table(sg_lib,file = "sg_lib.txt",row.names = FALSE)
                                  
