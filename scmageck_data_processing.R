#library packages

library(Seurat)
library(stringi)
library(stringr)
library(plyr)

#get parameter
args<- commandArgs(T)
mtx_dir<- args[1]
project<- args[2]
species<- args[3]

#read files
mtx<- read.delim(gzfile(mtx_dir),header=T,row.names=1)

#split matrix into two parts,sgRNA information and count matrix
sgRNA<- mtx[(nrow(mtx)-54):nrow(mtx),4:ncol(mtx)]
mtx<- mtx[2:(nrow(mtx)-55),4:ncol(mtx)]

#create sg_lib
sg_num<- apply(data.frame(sgRNA),2,function(x) sum(x!=0))
sg_info<- data.frame()
for(i in 1:54){
    sg<- sgRNA[i,]
    l<- colnames(sg)[which(sg!=0)]
    for(j in l){
        sg_info<- rbind(sg_info,c(j,rownames(sg)))
    }
}
colnames(sg_info)<- c("cell","barcode")
gene<- data.frame(lapply(strsplit(sg_info$barcode,"-"),function(x) head(x,1)))
gene<- as.vector(t(gene[1,]))
gene<- str_replace(gene,"Negative","NTC")
sg_names<- sg_info$barcode
                         
#rename sgRNA
rename<- data.frame(unique(cbind(sg_names,gene)))
rownames(rename)<- rename$sg_names
for(i in unique(gene)){
    a<- subset(rename,gene==i)
    a$sg_names<- paste(i,which(a$gene==i),sep="_sgRNA")
    for(j in rownames(a)){
        sg_names<- str_replace(sg_names,j,a[j,1])
    }
}
sg_lib<- data.frame(t(rbind(sg_info$cell,sg_names,gene)))
names(sg_lib)<- c("cell","barcode","gene")

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
sg_lib<- subset(sg_lib,cell %in% subset(sg_in_cell,freq==1)[,1])
write.table(sg_lib,file = "sg_lib.txt",row.names = FALSE)
