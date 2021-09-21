#' @export
#' @import Signac

ATAC_Add_meta_data <- function(sg_dir, mtx_dir, fragments, replicate = 1){
    
    #Add "sgRNA_num" and "perturbations"
    
    peak <- Add_meta_data(sg_dir, mtx_dir, cal.mt = FALSE, replicate = replicate)

    #calculate FRiP
    
    if(!("FRiP" %in% colnames(peak@meta.data))){
        if(is.character(fragments)){
            if(fragments %in% colnames(peak@meta.data)){
                peak <- FRiP(peak, "peaks", total.fragments = fragments)
            }else{
                frag <- CountFragments(fragments = fragments, cells = colnames(peak))
                peak$fragments <- frag$reads_count
                peak <- FRiP(peak, "peaks", total.fragments = "fragments")
            }
        }else{
            stop("Please provide the path of fragments file or meta data names of total fragments counts")
        }
    }
    return(peak)
}

#' @export

ATAC_scQC <- function(mtx_dir, prefix = "./", label = "", peak_frac = 0.01, nFeature = c(200, 500000), nCount = 1000, FRiP = 0.1, blank_NTC = FALSE){
      #read file
  if (is.character(mtx_dir)) {
    message(paste("Reading RDS file:", mtx_dir))
    perturb <- readRDS(mtx_dir)
  } else {
    perturb <- mtx_dir
  }
    
    perturb$nFeature_peak <- perturb[[paste("nFeature_", perturb@active.assay, sep = "")]][, 1]
    perturb$nCount_peak <- perturb[[paste("nCount_", perturb@active.assay, sep = "")]][, 1]

  #QC plot of the single cell matrix
  pdf(file = file.path(prefix, paste(label, "raw_matrix_quality_vlnplot.pdf", sep = "")))
  p1 <- VlnPlot(perturb, features = c("nFeature_peak", "nCount_peak", "FRiP"), ncol = 3, pt.size = 0.1)
  print(p1)
  dev.off()

    
  #filter cells with low quality
  if(blank_NTC == TRUE){
    perturb_QC <- subset(perturb,
                        nFeature_peak <= nFeature[2] &
                          nFeature_peak >= nFeature[1] &
                          nCount_peak >= nCount &
                          FRiP >= FRiP)
  }else{
    perturb_QC <- subset(perturb,
                         nFeature_peak <= nFeature[2] &
                         nFeature_peak >= nFeature[1] &
                         nCount_peak >= nCount &
                         FRiP >= FRiP &
                         perturbations != 'blank')
  }
  perturb_QC <- CreateSeuratObject(counts = GetAssayData(object = perturb_QC, slot = "counts"), assay = "peak",
                                  min.cells = peak_frac * ncol(perturb_QC), project = perturb@project.name)
    if("FRiP" %in% colnames(perturb@meta.data)){
        perturb_QC$FRiP <- perturb$FRiP
    }
    
    if("replicate" %in% colnames(perturb@meta.data)){
        perturb_QC$replicate <- perturb$replicate
    }
    
    if("perturbations" %in% colnames(perturb@meta.data)){
        perturb_QC$perturbations <- perturb$perturbations
    }
  pdf(file = file.path(prefix, paste(label, "QC_matrix_quality_vlnplot.pdf", sep = "")))
  p2 <- VlnPlot(perturb_QC, features = c("nFeature_peak", "nCount_peak", "FRiP"), ncol = 3, pt.size = 0.1)
  print(p2)
  dev.off()

  #saveRDS(perturb_QC,file=file.path(prefix, "perturb_QC.rds"))
  return(perturb_QC)
}

#' @export

CalculateGeneActivity <- function(mtx_dir, fragments, species = "Hs", version = "v75", gene_type = "Symbol", protein_coding = TRUE, pro_up = 3000, pro_down = 0){
    
    #get promoter region
    
    pro <- GetPromoter(species, version, gene_type, protein_coding, pro_up, pro_down)
    genebodyandpromoter.coords <- pro[[1]]
    gene.key <- pro[[2]]
    
    #get count matrix
    
    if (is.character(mtx_dir)) {
        message(paste("Reading RDS file:", mtx_dir))
        perturb <- readRDS(mtx_dir)
    } else {
        perturb <- mtx_dir
    }

    #get fragments
    
    if(is.character(fragments)){
        fragments<- CreateFragmentObject(fragments, cells = colnames(perturb))
    }
    
    #calculate gene activity
    
    gene.activity<- FeatureMatrix(fragments = fragments, features = genebodyandpromoter.coords, cells = colnames(perturb))
    
    #generate gene activity matrix
    
    rownames(gene.activity) <- gene.key[rownames(gene.activity)]
    perturb_RNA<- CreateSeuratObject(counts = gene.activity, project = perturb@project.name)
    if("replicate" %in% colnames(perturb@meta.data)){
        replicate <- perturb$replicate
        perturb_RNA$replicate <- replicate
    }
    
    if("perturbations" %in% colnames(perturb@meta.data)){
        perturb_RNA$perturbations <- perturb$perturbations
    }
    return(perturb_RNA)
}

#' @export

GetPromoter <- function(species = "Hs", version = "v75", gene_type = "Symbol", protein_coding = TRUE, pro_up = 3000, pro_down = 0){

    #get gene ranges from selected reference
    if(species == "Hs"){
        if(version == "v75"){
            gene.ranges <- genes(EnsDb.Hsapiens.v75)
        }else if(version == "v79"){
            gene.ranges <- genes(EnsDb.Hsapiens.v79)
        }else if(version == "v86"){
            gene.ranges <- genes(EnsDb.Hsapiens.v86)
        }  
    }else if(species == "Mm"){
        if(version == "v75"){
            gene.ranges <- genes(EnsDb.Mmusculus.v75)
        }else if(version == "v79"){
            gene.ranges <- genes(EnsDb.Mmusculus.v79)
        }
    }

    seqlevelsStyle(gene.ranges) <- 'UCSC'
    
    if(protein_coding == TRUE){
        gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
    }

    gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

    genebodyandpromoter.coords <- suppressWarnings(trim(Extend(x = gene.ranges, upstream = pro_up, downstream = pro_down)))

    if(gene_type == "Symbol"){
        gene.key <- genebodyandpromoter.coords$gene_name
    }else if(gene_type == "Ensembl"){
        gene.key <- genebodyandpromoter.coords$gene_id
    }else{
        warning("This function only support 'gene Symbol' and 'Ensembl id' as gene names, using gene Symbol instead")
        gene.key <- genebodyandpromoter.coords$gene_name
    }
    

    names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
    
    return(list(genebodyandpromoter.coords, gene.key))
}

