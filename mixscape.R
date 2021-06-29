#library packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(reshape2)
library(plyr)

#get parameters
args<- commandArgs(T)
mtx_dir<- args[1]
sg_dir<- args[2]
NTC<- args[3]

custom_theme <- theme(
  plot.title = element_text(size=16, hjust = 0.5), 
  legend.key.size = unit(0.7, "cm"), 
  legend.text = element_text(size = 14))

eccite<- readRDS(mtx_dir)
eccite <- subset(eccite, 
                     nFeature_RNA <= 5000 & 
                     nFeature_RNA >= 200 & 
                     nCount_RNA >= 1000 &
                     percent.mt <= 10)
#add barcode information
sg_lib<- read.table(sg_dir,header=T)
sg_in_cell<- data.frame(plyr::count(sg_lib$cell))
sg_lib<- subset(sg_lib,cell %in% subset(sg_in_cell,freq==1)[,1])
rownames(sg_lib)<- sg_lib$cell
eccite<- eccite[,colnames(eccite) %in% sg_lib$cell]
NT<- data.frame(cell = colnames(eccite),nt = "a")
NT$nt<- sg_lib[NT$cell,]$barcode
eccite <- AddMetaData(eccite,as.factor(NT$nt),col.name = "NT")

# Prepare RNA assay for dimensionality reduction: 
# Normalize data, find variable features and scale data.
DefaultAssay(object = eccite) <- 'RNA'
eccite <- NormalizeData(object = eccite) %>% FindVariableFeatures() %>% ScaleData()

# Run Principle Component Analysis (PCA) to reduce the dimensionality of the data.
eccite <- RunPCA(object = eccite)

# Run Uniform Manifold Approximation and Projection (UMAP) to visualize clustering in 2-D.
eccite <- RunUMAP(object = eccite, dims = 1:40)

#calculate cell cycle gene
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
eccite<- CellCycleScoring(eccite, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#transform label
eccite$gene <- eccite$perturbations
eccite$perturbations<- as.character(eccite$perturbations)
eccite$perturbations[which(eccite$perturbations!=NTC)] <- "perturbed"

# Generate plots to check if clustering is driven by biological replicate ID, 
# cell cycle phase or target gene class.
p1 <- DimPlot(
  object = eccite, 
  group.by = 'replicate', 
  label = F, 
  pt.size = 0.2, 
  reduction = "umap", cols = "Dark2", repel = T) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("Biological Replicate") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  custom_theme

p2 <- DimPlot(
  object = eccite, 
  group.by = 'Phase', 
  label = F
, pt.size = 0.2, 
  reduction = "umap", repel = T) + 
  ggtitle("Cell Cycle Phase") +
  ylab("UMAP 2") +
  xlab("UMAP 1") +
  custom_theme

p3 <- DimPlot(
  object = eccite, 
  group.by = 'perturbations', 
  pt.size = 0.2, 
  reduction = "umap", 
  split.by = "perturbations", 
  ncol = 1, 
  cols = c("grey39","goldenrod3")) + 
  ggtitle("Perturbation Status") +
  ylab("UMAP 2") +
  xlab("UMAP 1") +
  custom_theme

#save plots.
pdf("mixscape_before.pdf")
((p1 / p2 + plot_layout(guides = 'auto')) | p3 )
dev.off()

#calculate perturb signature
eccite<- CalcPerturbSig(
  object = eccite, 
  assay = "RNA", 
  gd.class ="gene", 
  nt.cell.class = NTC, 
  reduction = "pca", 
  ndims = 40, 
  num.neighbors = 20, 
  split.by = "replicate", 
  new.assay.name = "PRTB")

#### Prepare PRTB assay for dimensionality reduction: 
# Normalize data, find variable features and center data.
DefaultAssay(object = eccite) <- 'PRTB'

# Use variable features from RNA assay.
VariableFeatures(object = eccite) <- VariableFeatures(object = eccite[["RNA"]])
eccite <- ScaleData(object = eccite, do.scale = F, do.center = T)

# Run PCA to reduce the dimensionality of the data.
eccite <- RunPCA(object = eccite, reduction.key = 'prtbpca', reduction.name = 'prtbpca')

# Run UMAP to visualize clustering in 2-D.
eccite <- RunUMAP(
  object = eccite, 
  dims = 1:40, 
  reduction = 'prtbpca', 
  reduction.key = 'prtbumap', 
  reduction.name = 'prtbumap')

# Generate plots to check if clustering is driven by biological replicate ID, 
# cell cycle phase or target gene class.
q1 <- DimPlot(
  object = eccite, 
  group.by = 'replicate', 
  reduction = 'prtbumap', 
  pt.size = 0.2, cols = "Dark2", label = F, repel = T) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("Biological Replicate") +
  ylab("UMAP 2") +
  xlab("UMAP 1") +
  custom_theme

q2 <- DimPlot(
  object = eccite, 
  group.by = 'Phase', 
  reduction = 'prtbumap', 
  pt.size = 0.2, label = F, repel = T) +
  ggtitle("Cell Cycle Phase") +
  ylab("UMAP 2") +
  xlab("UMAP 1") + 
  custom_theme

q3 <- DimPlot(
  object = eccite,
  group.by = 'perturbations',
  reduction = 'prtbumap', 
  split.by = "perturbations", 
  ncol = 1, 
  pt.size = 0.2, 
  cols = c("grey39","goldenrod3")) +
  ggtitle("Perturbation Status") +
  ylab("UMAP 2") +
  xlab("UMAP 1") +
  custom_theme

#save plots.
pdf("mixscape_after.pdf")
(q1 / q2 + plot_layout(guides = 'auto') | q3)
dev.off()

# Run mixscape.
eccite <- RunMixscape(
  object = eccite, 
  assay = "PRTB", 
  slot = "scale.data", 
  labels = "gene", 
  nt.class.name = NTC, 
  min.de.genes = 5, 
  iter.num = 10, 
  de.assay = "RNA", 
  verbose = F,
  prtb.type = "KO")

saveRDS(eccite,"mixscape.rds")

if(length(unique(eccite$mixscape_class.global))==3){
    # Calculate percentage of KO cells for all target gene classes.
    df <- prop.table(table(eccite$mixscape_class.global, eccite$NT),2)
    df2 <- reshape2::melt(df)
    df2$Var2 <- as.character(df2$Var2)
    test <- df2[which(df2$Var1 == "KO"),]
    test <- test[order(test$value, decreasing = T),]
    new.levels <- test$Var2
    df2$Var2 <- factor(df2$Var2, levels = new.levels )
    df2$Var1 <- factor(df2$Var1, levels = c(NTC, "NP", "KO"))
    df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "_sgRNA")[[1]][1])
    df2$guide_number <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "_sgRNA")[[1]][2])
    df3 <- df2[-c(which(df2$gene == NTC)),]
    df4 <- data.frame(unique(df3[,c(1,3,4)]))
    df5 <- subset(df3,gene %in% subset(df4,Var1=="KO" & value !=0)$gene)
    #only remain genes with non-zero KO
    l <- length(unique(df5$gene))
    pdf("mixscape_KO_percent.pdf")
    if(l < 12){
    p <- ggplot(df5, aes(x = guide_number, y = value*100, fill= Var1)) + 
                               geom_bar(stat= "identity") + 
                               theme_classic()+ 
                               scale_fill_manual(values = c("grey49", "grey79","coral1")) + 
                               ylab("% of cells") +
                               xlab("sgRNA")
    p + theme(axis.text.x = element_text(size = 18, hjust = 1), 
           axis.text.y = element_text(size = 18), 
           axis.title = element_text(size = 16), 
           strip.text = element_text(size=8, face = "bold")) + 
  facet_wrap(vars(gene),ncol = 3, scales = "free") +
  labs(fill = "mixscape class") +theme(legend.title = element_text(size = 14),
                                       legend.text = element_text(size = 12))
        }else{
        for(i in 1:(l%%12 + 1)){
            df6<- subset(df5,gene %in% unique(df5$gene)[(i*12-11):i*12])
            p <- ggplot(df6, aes(x = guide_number, y = value*100, fill= Var1)) + 
                               geom_bar(stat= "identity") + 
                               theme_classic()+ 
                               scale_fill_manual(values = c("grey49", "grey79","coral1")) + 
                               ylab("% of cells") +
                               xlab("sgRNA")
            p + theme(axis.text.x = element_text(size = 18, hjust = 1), 
           axis.text.y = element_text(size = 18), 
           axis.title = element_text(size = 16), 
           strip.text = element_text(size=8, face = "bold")) + 
            facet_wrap(vars(gene),ncol = 3, scales = "free") +
            labs(fill = "mixscape class") +theme(legend.title = element_text(size = 14),
                                       legend.text = element_text(size = 12))
        }
    }
    dev.off()
                           }