#' function definitions ##### Add perturbations, sgRNA_num and
#' percent.mt of each cells to the Seurat object.
#' For other downstream function of SCREEN, this step is always needed
#' @export

Add_meta_data<- function(sg_dir, mtx_dir, species = "Hs", prefix = "./", label = ""){
  #read files
  if (is.character(mtx_dir)) {
    message(paste("Reading RDS file:", mtx_dir))
    mtx <- readRDS(mtx_dir)
  } else {
    mtx <- mtx_dir
  }
  if (is.character(sg_dir)) {
    message(paste("Reading sgRNA lib file:", sg_dir))
    sg_lib <- read.table(sg_dir,header=T)
  } else {
    sg_lib <- sg_dir
  }

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
      label[x[,1]]<- sg_lib_filtered[which(sg_lib_filtered$cell==x[,1]),]$gene}else{
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

  saveRDS(mtx,file = file.path(prefix, paste(label, "perturb.rds", sep = "")))
  write.table(sg_lib,file = file.path(prefix, paste(label, "sg_lib_all.txt", sep = "")),row.names=FALSE,quote=FALSE)
  sg_in_cell<- data.frame(plyr::count(sg_lib$cell))
  sg_lib<- subset(sg_lib,cell %in% subset(sg_in_cell,freq==1)[,1])
  write.table(sg_lib,file = file.path(prefix, paste(label, "sg_lib.txt", sep = "")), row.names = FALSE,quote=FALSE)
  return(mtx)
}
