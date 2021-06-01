#library packages
library(Matrix)
library(Seurat)
library(scMAGeCK)
library(stringi)
library(ggplot2)
library(plyr)
library(stringr)

#get parameter
args<- commandArgs(T)
mtx_dir<- args[1]
barcode_dir<- args[2]
gene_dir<- args[3]
sg_dir<- args[4]
project<- args[5]
species<- args[6]

mtx<- readMM(gzfile(mtx_dir))
#if input barcode is a csv file
barcodes<- read.csv(gzfile(barcode_dir),header=T,row.names=1)
#if input genes is a csv file
genes<- read.csv(gzfile(gene_dir),header=T,row.names=1)
#if input sgRNA dictionary is a csv file
sg_dict<- read.csv(gzfile(sg_dir),header=F,row.names=1)

#get the gene name 
gene_symbol<- strsplit(genes[,1],"_")
gene_symbol<- data.frame(lapply(gene_symbol,function(x) x[2]))

#create matrix with cellnames and features
rownames(mtx)<- gene_symbol[1,]
colnames(mtx)<- barcodes[,1]

#create sg_lib
cell<- strsplit(sg_dict[,1],", ")
names(cell)<- rownames(sg_dict)
l_list<- lapply(cell, function(x) length(x))
sg_names<- unlist(lapply(names(l_list), function(x) rep(x,time=l_list[x])))
sg_info<- unlist(cell)
gene<- strsplit(sg_names,"_")
gene<- unlist(lapply(gene,function(x) x[2]))
gene<- gsub("INTERGENIC\\w*","NTC",gene)
gene<- gsub("sg","",gene)

#rename sgRNAs
rename<- data.frame(unique(cbind(sg_names,gene)))
rownames(rename)<- rename$sg_names
for(i in unique(gene)){
    a<- subset(rename,gene==i)
    a$sg_names<- paste(i,which(a$gene==i),sep="_sgRNA")
    for(j in rownames(a)){
        sg_names<- str_replace(sg_names,j,a[j,1])
    }
}
sg_names<- gsub("MouseNTC","NTC",sg_names)
sg_lib<- data.frame(t(rbind(sg_info,sg_names,gene)))
names(sg_lib)<- c("cell","barcode","gene")
sg_lib$gene<- str_replace(sg_lib$gene,"MouseNTC","NTC")

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