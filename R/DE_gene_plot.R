#' function definitions ##### Plot DE gene numbers of each perturbations
#' @export

DE_gene_plot<- function(score_dir, pval_dir, project = "perturb",
                        prefix = "./", p_val_cut = 0.05, score_cut = 0.5, label = "",
                        ylimit = c(-600, 600, 200)){
  score<- read.table(score_dir,header=T,row.names=1)
  p_val<- read.table(pval_dir,header=T,row.names=1)
  de_genes<- data.frame(non=rep(0, ncol(score)),
                        up=rep(0,ncol(score)), down=rep(0, ncol(score)))
  rownames(de_genes)<- colnames(score)
  for(i in colnames(score)){
    scmageck<- data.frame(t(rbind(score[,i],p_val[,i])))
    colnames(scmageck)<- c("score","p_val")
    scmageck$diff<- "non"
    scmageck$diff[scmageck$score > score_cut & scmageck$p_val < p_val_cut]<- "up"
    scmageck$diff[scmageck$score < score_cut & scmageck$p_val < p_val_cut]<- "down"
    a<- count(scmageck$diff)
    rownames(a)<- a$x
    for(j in a$x){
      de_genes[i,j]<- a[j,2]
    }
  }
  new<- data.frame(reshape2::melt(de_genes[,2:3]))
  new$factor<- rep(rownames(de_genes),2)
  colnames(new)<- c("diff","gene_numbers","factor")
  new1<- new
  new1$diff<- as.character(new1$diff)
  for(i in unique(new$factor)){
    a<- subset(new, factor==i)
    all<- sum(a$gene_numbers)
    new1<- rbind(new1, c("all",all,i))
  }
  new2<- subset(new1, diff=="all")
  new3<- subset(new,factor %in% head(new2[order(as.numeric(new2$gene_numbers), decreasing = T),],20)$factor)
  new3[which(new3$diff=="down"),]$gene_numbers<- -new3[which(new3$diff=="down"),]$gene_numbers
  p1<- ggplot(new3,aes(x = factor,y = gene_numbers))+
    geom_bar(stat = 'identity',aes(fill = diff)) +
    theme_classic() +
    labs(x="sgRNA",y="DE genes",title = project)+
    theme(text=element_text(size=12,face = "bold"),
          axis.text.x = element_text(size=10)) +
    theme(plot.title=element_text(hjust=0.5),
          axis.text.x = element_text(angle = 90,hjust = 0.5,vjust = 0.5))+
    scale_y_continuous(limits=c(ylimit[1],ylimit[2]),
                       breaks=seq(ylimit[1],ylimit[2],ylimit[3]),
                       labels = c(seq(ylimit[2],0,-ylimit[3]),
                                  seq(ylimit[3],ylimit[2],ylimit[3])))
  pdf(file.path(prefix,
                paste(label, "DE_gene_cutoff",score_cut,"_p",p_val_cut,".pdf",sep = "")))
  p1
  dev.off()
  return(p1)
}
