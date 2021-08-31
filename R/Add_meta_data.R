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
  lab<- rep("blank",times=ncol(mtx))
  sg_num<- rep(0,times=ncol(mtx))
  names(lab)<- colnames(mtx)
  names(sg_num)<- colnames(mtx)
  sg_in_cell<- data.frame(plyr::count(sg_lib_filtered$cell))

  #find unique label and multiple label
  for(i in 1:nrow(sg_in_cell)){
    x<- sg_in_cell[i,]
    if(x[,2]==1){
      lab[x[,1]]<- sg_lib_filtered[which(sg_lib_filtered$cell==x[,1]),]$gene}else{
        lab[x[,1]]<- "multiple"}
  }

  for(i in 1:nrow(sg_in_cell)){
    if(sg_in_cell[i,1] %in% names(sg_num)){
      sg_num[sg_in_cell[i,1]]<- sg_in_cell[i,2]
    }
  }
  #add the label information to Seurat object
  lab<- as.factor(lab)
  mtx<- AddMetaData(mtx, lab, col.name = "perturbations")
  mtx<- AddMetaData(mtx, sg_num, col.name = "sgRNA_num")

  #calculate percent.mt
  if(species=="Hs"){
    mtx[["percent.mt"]] <- PercentageFeatureSet(mtx, pattern = "^MT-")
  }else{
    mtx[["percent.mt"]] <- PercentageFeatureSet(mtx, pattern = "^mt-")
  }
  #saveRDS(mtx,file = file.path(prefix, paste(label, "perturb.rds", sep = "")))
  return(mtx)
}
