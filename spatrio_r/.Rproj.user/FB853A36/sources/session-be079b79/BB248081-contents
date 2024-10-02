#' Calcualte significant gene-motif combos for celltypes
#'
#' @param object Input data in SeuratObject format
#' @param motif_gene_df Dataframe of motifs and genes
#' @param aver_exp_all Dataframe of fold changes of gene expression
#' @param aver_chromvar_all Dataframe of fold changes of motif activity
#'
#' @export
#'
#' @import dplyr
gene_motif_analysis <- function(object,motif_gene_df,aver_exp_all=NULL,aver_chromvar_all=NULL){
  motif.ls = motif_gene_df$motif
  df.ls <- lapply(aver_chromvar_all$gene, function(motif) {
    #print(motif)
    gene <- motif_gene_df[motif_gene_df$motif == motif,]$gene_name
    #print(gene)

    aver_chromvar <- aver_chromvar_all[aver_chromvar_all$gene==motif,]
    aver_exp <- aver_exp_all[aver_exp_all$gene==gene,]

    # join the matrices and return null if gene expression is not detected for all cell types
    mat <- tryCatch(full_join(aver_chromvar, aver_exp, by = "cluster"),
                    error=function(e) NULL)
    if(!is.null(mat)) {
      df <- data.frame(chromvar=mat$avg_log2FC.x, rna=mat$avg_log2FC.y, celltype=mat$cluster,
                       motif=mat$gene.x, gene=mat$gene.y, chrom_pval=mat$p_val_adj.x, gene_pval=mat$p_val_adj.y)
    } else {
      return(NULL)
    }
  })
  # aggregate the motif-gene df
  df <- bind_rows(df.ls)
  df <- dplyr::arrange(df, gene, motif)

  # filter the motif-gene df for pairs with significant changes in gene expression and chromvar activity
  df.f <- dplyr::filter(df, chrom_pval < 0.05, gene_pval < 0.05)
  df.f <- na.omit(df.f)

  # create a column for each unique gene-motif combo
  df.f <- dplyr::mutate(df.f, gene_motif = paste0(gene,"_",motif))

  # compute a pearson r2 and pval for each motif-gene combo across all celltypes
  combos <- unique(df.f$gene_motif)
  library(ggpmisc)
  df.pearson <- lapply(combos, function(combo) {
    require(ggrepel)

    df <- dplyr::filter(df.f, gene_motif == combo)
    max_exp <- max(abs(df$rna), na.rm=TRUE) + 1
    max_chrom <- max(abs(df$chromvar), na.rm=TRUE) + 1
    pearson <- tryCatch(cor.test(df$chromvar, df$rna, method="pearson", conf.level=0.95), error=function(e) NULL)
    df$corr <- tryCatch(signif(pearson$estimate,2), error=function(e) NULL)
    df$pval <- tryCatch(signif(pearson$p.value,2), error=function(e) NULL)
    df$max_exp <- max_exp
    df$max_chrom <-max_chrom
    df$num_celltypes <- length(unique(df$celltype))
    return(df)
  })

  # aggregate a final df with all the available pearson R2 coefficients and pvals
  df.final <- bind_rows(df.pearson)
  df.final <- dplyr::distinct(df.final, celltype, gene_motif, .keep_all = TRUE)
  df.sig <- dplyr::filter(df.final, pval < 0.05) %>%
    arrange(-num_celltypes, -corr)

  df.final$celltype <- as.factor(df.final$celltype)
  levels(df.final$celltype) <- levels(Idents(object))


  df.sig$celltype <- as.factor(df.sig$celltype)
  levels(df.sig$celltype) <- levels(Idents(object))

  return(df.sig)
}


#' Calcualte the score of modules
#'
#' @param object Input data in SeuratObject format
#' @param input_list Modules to be plotted
#' @param nbin Number of bins of aggregate expression levels for all analyzed feature
#' @param ctrl Number of control features
#' @param seed Random seed
#' @param clean Whether to clear the previous score
#' @param name Prefix of module serial number
#'
#' @export
#'
#' @import magrittr
#' @import pheatmap
#' @import dplyr
spatrio_score<- function(object,
                         input_list,
                         nbin=10,
                         ctrl=50,
                         seed=42,
                         clean=T,
                         name="Module"){
  if(clean){
    object@meta.data[,grep(name,colnames(object@meta.data))]<-NULL
  }
  object <- AddModuleScore(object, features=input_list$group, name=name, nbin=nbin, ctrl=ctrl, seed=seed)
  return(object)
}

#' Infer the gene-gene celluar communication by SpaTrio
#'
#' @param object Input data in SeuratObject format
#' @param assay Name of assay to use
#' @param xcoord x coordinates of all cells
#' @param ycoord y coordinates of all cells
#' @param species A character indicating the species of the spatial transcriptomics data. Either ‘Human’ or ‘Mouse’
#' @param celltype A character specifying the cell type of ST data. To bypass the deconvolution step and directly infer cell-cell communication, please provide the cell type. Default is NULL
#'
#' @export
#'
#' @import SpaTalk
gene_cci_analysis<-function(object,
                            assay="RNA",
                            xcoord=NULL,
                            ycoord=NULL,
                            species="Human",
                            celltype=NULL,
                            lrpairs=NULL,
                            pathways=NULL,...
){
  if (!is(object, "Seurat")) {
    stop("Invalid class for object: must be 'Seurat'!")
  }
  if (is.null(xcoord)) {
    stop("Please enter coordinate information !")
  }
  if (is.null(ycoord)) {
    stop("Please enter coordinate information !")
  }
  if (is.null(celltype)) {
    stop("Please enter cell type information !")
  }
  meta=data.frame(cell=paste("C",1:length(object$orig.ident),sep = ""),x=xcoord,y=ycoord)
  count<-as.data.frame(object@assays$RNA@counts)
  colnames(count)<-paste("C",1:ncol(object),sep = "")
  count<-as.matrix(count)
  obj <- SpaTalk::createSpaTalk(st_data = count,
                       st_meta = meta,
                       species = species,
                       if_st_is_sc = T,
                       spot_max_cell = 1,
                       celltype = as.character(celltype))
  obj <- SpaTalk::find_lr_path(object = obj, lrpairs = lrpairs, pathways = pathways)
  obj <- SpaTalk::dec_cci_all(object = obj,...)
  return(obj)
}

#' Infer the protein-protein celluar communication
#'
#' @param object Input data in SeuratObject format
#' @param input_list Spatrio module analysis results
#' @param id Dataframe for id transfer
#' @param lr Dataframe of cell-cell communication information
#' @param celltype_sender
#' @param celltype_receiver
#' @param celltype_pair
#' @param lr_pair
#' @param module_sender
#' @param module_receiver
#' @param mode
#' @param transfer
#'
#' @export
#'
#' @import dplyr
#' @import ggplot2
adt_cci_analysis<-function(object,
                           input_list=NULL,
                           id=NULL,
                           lr=NULL,
                           celltype_sender=NA,
                           celltype_receiver=NA,
                           celltype_pair=NULL,
                           lr_pair=NULL,
                           module_sender=NA,
                           module_receiver=NA,
                           mode="1",
                           transfer=F){
  module_list<-input_list
  module_df=data.frame(module=NULL,orig=NULL)
  for (i in names(module_list$group)) {
    module_sub<-data.frame(module=i,orig=module_list$group[[i]])
    module_df <-rbind(module_df,module_sub)
  }

  if (transfer){
    module_df<-left_join(module_df,id,by="orig")
    module_df<-module_df[,c("module","orig","transfer")]
    module_df$transfer[is.na(module_df$transfer)]<-module_df$orig[is.na(module_df$transfer)]

    count<-object@assays$ADT@counts
    count<-as.data.frame(count)
    adt_name<-rownames(count)
    for (i in 1:nrow(count)) {
      count$adt[i]<-id$transfer[id$orig==rownames(count)[i]][1]
    }
    count<-aggregate(x=count[,1:c(ncol(count)-1)], by=list(count$adt),sum)
    index<-count$Group.1
    count<-count[,2:ncol(count)]
    rownames(count)<-index
    count<-as.matrix(count)
    object@assays$ADT<-CreateAssayObject(count)
    object@assays[["ADT"]]@key<-"adt_"
    DefaultAssay(object)<-"ADT"
    VariableFeatures(object)<-rownames(object)
    object<-NormalizeData(object, normalization.method = "CLR", margin = 2)%>%ScaleData()
  }else{
    module_df$transfer<-module_df$orig
  }

  if(mode=="1"){
    tmp<-lr[lr$ligand_gene_symbol%in%module_df$transfer[module_df$module==module_sender],]
    lrpair_df<-tmp[tmp$receptor_gene_symbol%in%module_df$transfer[module_df$module==module_receiver],]
    lrpair_df<-lrpair_df[,1:3]
    av<-AverageExpression(object,features = unique(c(lrpair_df$ligand_gene_symbol,lrpair_df$receptor_gene_symbol)))$ADT
    lrpair_df$ligand_exp<-av[lrpair_df$ligand_gene_symbol,celltype_sender]
    lrpair_df$receptor_exp<-av[lrpair_df$receptor_gene_symbol,celltype_receiver]
    lrpair_df$score<-sqrt(lrpair_df$ligand_exp*lrpair_df$receptor_exp)
    lrpair_df$score<-lrpair_df$score/max(lrpair_df$score)

    lrpair_df<-lrpair_df[order(lrpair_df$score,decreasing = T),]
    lrpair_df$ligand_gene_symbol=factor(lrpair_df$ligand_gene_symbol,levels = rev(unique(lrpair_df$ligand_gene_symbol)))
    lrpair_df$receptor_gene_symbol=factor(lrpair_df$receptor_gene_symbol,levels = unique(lrpair_df$receptor_gene_symbol))

  }else if(mode=="2"){

    a<-lr[lr$ligand_gene_symbol%in%module_df$transfer,]
    b<-a[a$receptor_gene_symbol%in%module_df$transfer,]
    b<-b[,1:3]
    for (i in 1:nrow(b)) {
      b$sender_module[i]<-module_df$module[which(module_df$transfer==b$ligand_gene_symbol[i])]
    }
    for (i in 1:nrow(b)) {
      b$reciever_module[i]<-module_df$module[which(module_df$transfer==b$receptor_gene_symbol[i])]
    }
    av<-as.data.frame(AverageExpression(object,features = unique(c(b$ligand_gene_symbol,b$receptor_gene_symbol)),assays = "ADT")$ADT)
    for (i in 1:nrow(b)) {
      b$sender[i]<-names(which.max(av[b$ligand_gene_symbol[i],]))
      b$reciever[i]<-names(which.max(av[b$receptor_gene_symbol[i],]))
      b$ligand_exp[i]<-as.numeric(max(av[b$ligand_gene_symbol[i],]))
      b$receptor_exp[i]<-as.numeric(max(av[b$receptor_gene_symbol[i],]))
    }
    b$score<-sqrt(b$ligand_exp*b$receptor_exp)
    lrpair_df<-as.data.frame(b)
    lrpair_df$pair<-paste(lrpair_df$sender,lrpair_df$reciever,sep = "_")
    lrpair_df$lr<-paste(lrpair_df$ligand_gene_symbol,lrpair_df$receptor_gene_symbol,sep = "_")
    if (!is.null(celltype_pair)) {
      lrpair_df<-lrpair_df[lrpair_df$pair%in%celltype_pair,]
    }
    if (!is.null(lr_pair)) {
      lrpair_df<-lrpair_df[lrpair_df$lr%in%lr_pair,]
    }
    lrpair_df$pair=factor(lrpair_df$pair,levels =celltype_pair)
    lrpair_df<-lrpair_df[order(lrpair_df$pair,decreasing = T),]
    lrpair_df$lr=factor(lrpair_df$lr,levels = lr_pair)
    lrpair_df$score<-lrpair_df$score/max(lrpair_df$score)
    lrpair_df<-lrpair_df[order(lrpair_df$score),]
  }
  return(lrpair_df)
}

#' Scaling function
#'
#' @param mode
#'
#' @export
range01 <- function(x){(x-min(x))/(max(x)-min(x))}


