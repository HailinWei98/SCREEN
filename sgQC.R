library(Seurat)
library(ggplot2)
library(plyr)

#get parameter
args<- commandArgs(T)
rds_dir<- args[1]
sg_dir<- args[2]

#read files
sg_lib<- read.table(sg_dir,header=T)
mtx<- readRDS(rds_dir)

#remove cells in sgRNA library that are not included in matrix
sg_lib_filtered<- subset(sg_lib,cell %in% intersect(sg_lib$cell,colnames(mtx)))

QC_mtx<- subset(mtx, 
                     nFeature_RNA <= 5000 & 
                     nFeature_RNA >= 200 & 
                     nCount_RNA >= 1000 &
                     percent.mt <= 10)
QC_sg_lib_filtered<- subset(sg_lib,cell %in% intersect(sg_lib$cell,colnames(QC_mtx)))

#count sgRNA in each cell;count cells of each sgRNA

#after QC filtered
sg_count_QC<- count(subset(QC_sg_lib_filtered,cell %in% colnames(QC_mtx))$barcode)
colnames(sg_count_QC)<- c("sgRNA","freq")

g1<- ggplot(sg_count_QC,mapping = aes(x=reorder(sgRNA,-freq),y=freq)) + 
  geom_bar(stat="identity",color="#e9ecef",fill="#69b3a2") + theme_classic() + 
  labs(x="sgRNA",y="cell numbers",title = "cell numbers of sgRNA") + 
  theme(plot.title=element_text(hjust=0.5),axis.text.x = element_text(angle = 90,hjust = 0.5,vjust = 0.5))

g2<- ggplot(QC_mtx[["sgRNA_num"]],aes(x = sgRNA_num)) + 
  geom_histogram(binwidth = 1,color="#e9ecef",fill="#69b3a2") + 
  theme_classic() + 
  labs(x="perturbed numbers",y="cell numbers",title = "sgRNA number in each cell") +
  theme(plot.title=element_text(hjust=0.5)) + 
  geom_text(aes(label=as.character(..count..)),stat="bin",binwidth=1,vjust=-0.5) + 
  scale_x_continuous(breaks=seq(min(QC_mtx[["sgRNA_num"]]),max(QC_mtx[["sgRNA_num"]]),1))

#before QC filtered
sg_count<- count(subset(sg_lib_filtered,cell %in% colnames(mtx))$barcode)
colnames(sg_count)<- c("sgRNA","freq")

g3<- ggplot(sg_count,mapping = aes(x=reorder(sgRNA,-freq),y=freq)) + 
  geom_bar(stat="identity",color="#e9ecef",fill="#69b3a2") + theme_classic() + 
  labs(x="sgRNA",y="cell numbers",title = "cell numbers of sgRNA") + 
  theme(plot.title=element_text(hjust=0.5),axis.text.x = element_text(angle = 90,hjust = 0.5,vjust = 0.5))

g4<- ggplot(mtx[["sgRNA_num"]],aes(x = sgRNA_num)) + 
  geom_histogram(binwidth = 1,color="#e9ecef",fill="#69b3a2") + 
  theme_classic() + 
  labs(x="perturbed numbers",y="cell numbers",title = "sgRNA number in each cell") +
  theme(plot.title=element_text(hjust=0.5)) + 
  geom_text(aes(label=as.character(..count..)),stat="bin",binwidth=1,vjust=-0.5) + 
  scale_x_continuous(breaks=seq(min(mtx[["sgRNA_num"]]),max(mtx[["sgRNA_num"]]),1))

#save plot
pdf("sgRNA_quality_before_QC.pdf")
g3;g4
dev.off()

pdf("sgRNA_quality_after_QC.pdf")
g1;g2
dev.off()