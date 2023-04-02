#' Visualization of spatrio results using heatmap
#'
#' @param input_list Spatrio module analysis results
#' @param group Modules to be plotted
#' @param main The title of the plot
#' @param cols Vector of colors used in heatmap
#' @param ann_cols Group color
#' @param name Prefix of module serial number
#'
#' @export
#'
#' @import magrittr
#' @import pheatmap
#' @import dplyr
spatrio_heatmap<- function(input_list,
                       group=NULL,
                       main=NA,
                       cols=c("#00578C","#8ED6EC","white","#FF7B67","#B00028"),
                       ann_cols=NA,
                       name="Module",...){
  col.heatmap <- colorRampPalette(cols)(1000)
  if(is.null(group)){
    group=1:length(input_list$group)
  }
  output = data.frame(feature="",G="")[-1,]
  for (i in group){
    output_sub=data.frame(feature=c(input_list$group[[i]]), G=paste(name,i, sep = ""))
    output<-rbind(output,output_sub)
  }
  output = output%>%magrittr::set_rownames(.$feature) %>% dplyr::select(-1)
  colnames(output)<-"Group"

  if (sum(is.na(ann_cols))>0){
    ann_colors<-NA
  } else{
    g<-paste(name,group, sep = "")
    ann_colors=list(Group=ann_cols)
    names(ann_colors$Group)<-g
  }
  print(ann_colors)
  p<-pheatmap::pheatmap(input_list$weighted_cor[rownames(output), rownames(output)],
                        clustering_method='ward.D2',
                        annotation_row=output,
                        annotation_names_row=F,
                        show_rownames=F,
                        show_colnames=F,
                        treeheight_row=10,
                        treeheight_col=10,
                        annotation_legend = T,
                        fontsize=8,
                        color=col.heatmap,
                        main=main,
                        legend_labels=F,
                        annotation_colors=ann_colors,
                        ...)
  return(p)

}

#' Visualization of spatrio results using dotplot
#'
#' @param df Dataframe collated from the results of spatrio
#' @param x Dataframe column name used as x axis data
#' @param y Dataframe column name used as y axis data
#' @param fill Dataframe column name used to color box
#' @param xlab Label of x axis
#' @param ylab Label of y axis
#' @param compare List of group for comparision
#' @param method A character string indicating which method to be used for comaprision
#' @param label A character string specifying label type of comaparision
#' @param vjust Vertical alignment adjustment of text
#' @param cols Group color
#' @param title Title of plot
#'
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @import ggpubr
spatrio_boxplot<-function(df,
                      x=NA,
                      y=NA,
                      fill=NA,
                      xlab="Cluster",
                      ylab=NULL,
                      compare=NULL,
                      method = "wilcox.test",
                      label = "p.signif",
                      vjust=0.2,
                      cols=NULL,
                      title=NULL,
                      x_level=NULL){
  my_comparisons<-compare
  plot_df=data.frame(x=df[,x],y=df[,y])
  if (!is.null(x_level)) {
    plot_df$x=factor(plot_df$x,levels = x_level)

  }
  library(ggpubr)
  p<-ggplot(plot_df, aes(x = x, y = y,fill=x))+
    geom_boxplot(width=0.5,position=position_dodge(0.9))+
    geom_point(size=0)+theme_bw(base_rect_size=1)+
    theme(panel.grid =element_blank())+
    FontSize(x.title = 15, y.title = 15,y.text=15,x.text=15)+
    theme(axis.text.x = element_text(color ="black"),axis.text.y  = element_text(color ="black")) +
    guides(color = guide_legend(override.aes = list(size = 8)))+
    theme(legend.position='none')+
    labs(x=xlab,y=ylab,title=title)
  if(!is.null(compare)){
    p<-p+stat_compare_means(comparisons = my_comparisons,method = method,label = label,vjust=vjust)
  }
  if(!is.null(cols)){
    p<-p+scale_fill_manual(values = cols)
  }
  return(p)
}

#' Visualization of distance
#'
#' @param object Input data in SeuratObject format
#' @param coord Matrix of spatial coordinate
#' @param root The cluster id used as root when calculating distance
#' @param level The level of cell clusters
#'
#' @export
#'
#' @import ggplot2
dist_plot<-function(object,
                   coord=object@images$image@coordinates[,c("x","y")],
                   root="Cell type 1",
                   level= c("Cell type 1","Cell type 2","Cell type 3")){
  dis_spot <- dist(object@images$image@coordinates[, coord]) %>% as.matrix
  dis_spot<-melt(dis_spot)
  dis_spot<-dis_spot[dis_spot$Var2%in%colnames(subset(object,idents = root)),]
  meta=data.frame(Var1 = colnames(object),type=object@active.ident)
  dis_spot<-left_join(dis_spot,meta,by="Var1")
  #dis_spot$type<-gsub("cluster","Cell type ",dis_spot$type)
  dis_spot<-dis_spot[dis_spot$type%in%level,]
  dis_spot$type=factor(dis_spot$type,levels = level)
  ggplot(dis_spot, aes(x = type, y = value,fill=type))+
    geom_boxplot(width=0.5,position=position_dodge(0.9),outlier.shape=NA)+
    labs(x="", y = "Distance")+
    theme_bw(base_rect_size=1)+
    theme(panel.grid =element_blank())+NoLegend()+
    FontSize(x.title = 15, y.title = 15,y.text=15,x.text=15)+
    theme(legend.text = element_text(size=15),legend.title = element_text(size=15))+
    theme(legend.key = element_rect(color = NA, fill = NA),legend.key.size = unit(1, "cm"))+
    theme(axis.text.x = element_text(color ="black",angle=45,hjust = 1),axis.text.y  = element_text(color ="black"))
}

#' Visualization of spatial clustering data
#'
#' @param object Input data in SeuratObject format
#' @param pt.size.factor Scale the size of the spots
#' @param image.alpha Transparency of background image
#' @param stroke The width of the border around the spots
#' @param highlight.clsuter The cluster name of cells to highlight
#' @param highlight.cell The cells to highlight
#' @param shape Adjust the shape of the spots
#' @param subset The cluster name of cells to splot
#' @param bg The color to fill the background
#' @param cols.highlight Colors to highlight the cells
#' @param label.color Color of the label text
#' @param ... Extra parameters used in sdplot. Parameters passed to \code{\link[Seurat]{SpatialDimPlot}}
#'
#' @export
#'
#' @import ggplot2
#' @import Seurat
sdplot<-function (object,
                  pt.size.factor=1,
                  image.alpha = 0,
                  stroke = NA,
                  highlight.clsuter =NULL,
                  highlight.cell =NULL,
                  shape=21,
                  subset=NULL,
                  ident=NA,
                  bg="white",
                  cols.highlight = c("#DE2D26", "grey80"),
                  label.color = "white",...
) {
  if (!is.na(ident)){
    Idents(object)<-object@meta.data[[ident]]
  }
  if (!is.null(subset)){
    object<-subset(object,ident=subset)
  }
  if (!is.null(highlight.clsuter)){
    cell.id = colnames(subset(object,ident=highlight.clsuter))
  } else if (!is.null(highlight.cell)) {
    cell.id = highlight.cell
  } else {
    cell.id = NULL
  }
  p<-SpatialDimPlot(object,
                    cells.highlight = cell.id,
                    stroke = stroke,
                    cols.highlight = cols.highlight,
                    label.color = label.color,
                    pt.size.factor=pt.size.factor,image.alpha=image.alpha,...)+
    theme(plot.background=element_rect(fill = bg),panel.grid =element_blank())+
    guides(fill=guide_legend(title=NULL))

  p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=shape)
  return(p)
}

#' Visualization of spatial expression data
#'
#' @param object Input data in SeuratObject format
#' @param stroke The width of the border around the spots
#' @param features Features to plot
#' @param shape Adjust the shape of the spots
#' @param option A character string of color map option.
#' @param ncol Number of columns
#' @param pt.size.factor Scale the size of the spots
#' @param ... Extra parameters used in sdplot. Parameters passed to \code{\link[Seurat]{SpatialFeaturePlot}}
#'
#' @export
#'
#' @import ggplot2
#' @import Seurat
sfplot<-function (object,
                  stroke = NA,
                  features=NA,
                  shape=21,
                  option="A",
                  ncol=NULL,
                  pt.size.factor=1,...) {
  if (length(features)>1) {
    p=list()
    for (i in features) {
      p1<-SpatialFeaturePlot(object,stroke = stroke,features=i,pt.size.factor=pt.size.factor,...)+scale_fill_viridis_c(option = option)
      p1$layers[[1]]$aes_params <- c(p1$layers[[1]]$aes_params, shape=shape)
      p[[i]]<-p1
    }
    p<-CombinePlots(p,ncol=ncol)
    return(p)
  }else{
    p<-SpatialFeaturePlot(object,stroke = stroke,features=features,pt.size.factor=pt.size.factor,...)+scale_fill_viridis_c(option = option)
    p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=shape)
    return(p)
  }
}

#' Visualization of significant gene-motif combos for celltypes
#'
#' @param df.sig Dataframe of significant gene-motif combos for celltypes
#' @param plot Select the type of gene-motif combos to plot. "sig","pos", or "neg"
#' @param return_data Whether to return data
#'
#' @export
#'
#' @import ggplot2
gene_motif_plot <- function(df.sig,plot='sig',return_data=FALSE){
  df.pos <- dplyr::filter(df.sig, corr > 0)
  sig_pos <- unique(df.pos$gene_motif)
  df.neg <- dplyr::filter(df.sig, corr < 0)
  sig_neg <- unique(df.neg$gene_motif)
  if (plot=="sig") {
    pearson <- cor.test(df.sig$chromvar, df.sig$rna, method="pearson", conf.level=0.95)
    max_chrom <- max(abs(df.sig$max_chrom)) + 1
    max_exp <- max(abs(df.sig$max_exp)) + 1
    p <- ggplot(df.sig, aes(x=chromvar, y=rna)) +
      geom_smooth(method="lm", color = "darkgray") +
      geom_point(aes(color=celltype)) +
      xlab("chromVAR activity (avg_logFC)") +
      ylab("Gene expression (avg_logFC)") +
      ggtitle("Significant gene-motif combos for celltypes",
              subtitle=paste0("Pearson r^2=",signif(pearson$estimate,2)," pval=",signif(pearson$p.val,2))) +
      xlim(c(-max_chrom,max_chrom)) +
      ylim(c(-max_exp,max_exp)) +
      theme_minimal() +
      geom_hline(yintercept=0,lty=3) +theme_bw(base_rect_size=1)+theme(panel.grid =element_blank())+
      geom_vline(xintercept=0,lty=3) +scale_color_manual(values =rna_col[levels(df.sig$celltype)])
  }else if(plot=="pos"){
    pearson <- cor.test(df.pos$chromvar, df.pos$rna, method="pearson", conf.level=0.95)
    max_chrom <- max(abs(df.pos$max_chrom)) + 1
    max_exp <- max(abs(df.pos$max_exp)) + 1
    p <- ggplot(df.pos, aes(x=chromvar, y=rna)) +
      geom_smooth(method="lm", color = "darkgray") +
      geom_point(aes(color=celltype)) +
      xlab("chromVAR activity (avg_logFC)") +
      ylab("Gene expression (avg_logFC)") +
      ggtitle("Positive correlation gene-motif combos for celltypes",
              subtitle=paste0("Pearson r^2=",signif(pearson$estimate,2)," pval=",signif(pearson$p.val,2))) +
      xlim(c(-max_chrom,max_chrom)) +
      ylim(c(-max_exp,max_exp)) +
      theme_minimal() +
      geom_hline(yintercept=0,lty=3) +theme_bw(base_rect_size=1)+theme(panel.grid =element_blank())+
      geom_vline(xintercept=0,lty=3) +scale_color_manual(values =rna_col[levels(df.sig$celltype)])
  }else if(plot=="neg"){
    pearson <- cor.test(df.neg$chromvar, df.neg$rna, method="pearson", conf.level=0.95)
    max_chrom <- max(abs(df.neg$max_chrom)) + 1
    max_exp <- max(abs(df.neg$max_exp)) + 1
    p <- ggplot(df.neg, aes(x=chromvar, y=rna)) +
      geom_smooth(method="lm", color = "darkgray") +
      geom_point(aes(color=celltype)) +
      xlab("chromVAR activity (avg_logFC)") +
      ylab("Gene expression (avg_logFC)") +
      ggtitle("Negative correlation gene-motif combos for celltypes",
              subtitle=paste0("Pearson r^2=",signif(pearson$estimate,2)," pval=",signif(pearson$p.val,2))) +
      xlim(c(-max_chrom,max_chrom)) +
      ylim(c(-max_exp,max_exp)) +
      theme_minimal() +
      geom_hline(yintercept=0,lty=3) +theme_bw(base_rect_size=1)+theme(panel.grid =element_blank())+
      geom_vline(xintercept=0,lty=3) +scale_color_manual(values =rna_col[levels(df.sig$celltype)])
  }

  if (return_data) {
    out=list()
    out[["sig"]]=df.sig
    out[["pos"]]=df.pos
    out[["neg"]]=df.neg
    return(out)
  }else{
    return(p)
  }
}

#' Visualization of selected gene-motif combos for celltypes
#'
#' @param object Input data in SeuratObject format
#' @param motif_gene_df Dataframe of motifs and genes
#' @param feature_gene Gene name
#'
#' @export
#'
#' @import ggplot2
#' @import Rmisc
gene_motif_lineplot<-function(object,feature_gene=NULL,motif_gene_df){
  feature_motif<-motif_gene_df$motif[motif_gene_df$gene_name==feature_gene]
  exp<-data.frame(gene=object@assays[["RNA"]]@data[feature_gene,],layer=Idents(object),group="gene")
  exp$gene<-scale(exp$gene)
  act<-data.frame(gene=object@assays[["chromvar"]]@data[feature_motif,],layer=Idents(object),group="motif")
  act$gene<-scale(act$gene)
  data<-rbind(exp,act)
  colnames(data)[1]<-"gene"
  tgc <- Rmisc::summarySE(data, measurevar="gene", groupvars=c("layer","group"))
  tgc$gene[tgc$group=="motif"]<-range01(tgc$gene[tgc$group=="motif"])
  tgc$gene[tgc$group=="gene"]<-range01(tgc$gene[tgc$group=="gene"])
  p<-ggplot(tgc, aes(x=layer, y=gene,group=group,color=group)) +
    geom_errorbar(aes(ymin=gene-se, ymax=gene+se), width=.1,size=1) +
    geom_line(size=1) +geom_point()+
    theme_bw(base_rect_size=1)+
    theme(panel.grid =element_blank())+
    FontSize(x.title = 15, y.title = 15,y.text=15,x.text=15)+
    theme(legend.text = element_text(size=15),legend.title = element_blank())+
    theme(legend.key = element_rect(color = NA, fill = NA),legend.key.size = unit(1, "cm"))+
    labs(x = "", y = "Exp/Act level")+guides(fill=guide_legend(title=NULL))+
    theme(axis.text.x = element_text(color ="black",angle = 45,hjust = 1),axis.text.y  = element_text(color ="black"))
  return(p)
}

#' Visualization of protein-protein communication
#'
#' @param df Dataframe of protein-protein conmmunication inference
#' @param color
#' @param mode
#'
#' @export
#'
#' @import ggplot2
adt_cci_plot<-function(df,color="#337EDE",mode="1"){
  if(mode=="1"){
    p<-ggplot(data = df)+
      geom_point(mapping = aes(x = receptor_gene_symbol,y = ligand_gene_symbol,size = score),colour=color)+
      labs(size = "LR score",y = "Ligand",x="Receptor")+
      theme_bw()+
      FontSize(x.title = 15, y.title = 15,y.text=15,x.text=15)+
      theme(axis.text.x = element_text(angle = 45,hjust = 1))
    return(p)
  }else if (mode=="2"){

  }
  p<-ggplot(data = df)+
    geom_point(mapping = aes(pair,y = lr,size = score),colour=color)+
    labs(size = "LR score",y = "Ligand",x="Receptor")+
    theme_bw()+
    FontSize(x.title = 15, y.title = 15,y.text=15,x.text=15)+
    theme(axis.text.x = element_text(angle = 45,hjust = 1))
  return(p)
}


#' Visualization of gene-gene communication
#'
#' @param object Input data in SpaTalk format
#' @param celltype_sender
#' @param celltype_receiver
#' @param top_lrpairs
#' @param color
#' @param return_data
#' @param lr_pair
#' @param celltype_pair
#'
#' @export
#'
#' @import ggplot2
gene_cci_plot <- function(object,
                          celltype_sender=NULL,
                          celltype_receiver=NULL,
                          top_lrpairs = 20,
                          color = "#337EDE",
                          return_data = FALSE,
                          lr_pair=NULL,
                          celltype_pair=NULL) {
  if (!is(object, "SpaTalk")) {
    stop("Invalid class for object: must be 'SpaTalk'!")
  }


  if (!is.null(lr_pair)&!is.null(celltype_pair)) {
    lrpair <- as.data.frame(object@lrpair)
    lrpair<-lrpair[lrpair$lr_co_ratio_pvalue<0.05,]
    lrpair$pair<-paste(lrpair$celltype_sender,lrpair$celltype_receiver,sep = "_")
    lrpair$lr<-paste(lrpair$ligand,lrpair$receptor,sep = "_")

    lrpair<-lrpair[lrpair$pair%in%celltype_pair,]
    lrpair$pair=factor(lrpair$pair,levels = celltype_pair)

    lrpair<-lrpair[lrpair$lr%in%lr_pair,]
    lrpair$lr=factor(lrpair$lr,levels = rev(lr_pair))

    p<-ggplot(data = lrpair)+
      geom_point(mapping = aes(pair,y = lr,size = score),color=color)+
      labs(size = "LR score",y = "Ligand",x="Receptor")+
      theme_bw()+
      theme(axis.text.x = element_text(color ="black",angle = 45,hjust = 1))+
      FontSize(x.title = 15, y.title = 15,y.text=15,x.text=15)

  }else{
    celltype<-unique(c(object@lrpair$celltype_sender,obj@lrpair$celltype_receiver))
    cell_pair <- object@cellpair
    cell_pair <- cell_pair[[paste0(celltype_sender, " -- ", celltype_receiver)]]
    if (is.null(cell_pair)) {
      stop("No LR pairs found from the celltype_sender to celltype_receiver!")
    }

    lrpair <- object@lrpair
    lrpair <- lrpair[lrpair$celltype_sender == celltype_sender & lrpair$celltype_receiver == celltype_receiver, ]
    lrpair <- lrpair[order(-lrpair$score), ]
    if (nrow(lrpair) > top_lrpairs) {
      lrpair <- lrpair[1:top_lrpairs, ]
    }
    lrpair$ligand=factor(lrpair$ligand,levels = rev(unique(lrpair$ligand)))
    lrpair$receptor=factor(lrpair$receptor,levels = unique(lrpair$receptor))
    p<-ggplot(data = lrpair)+
      geom_point(mapping = aes(x = receptor,y = ligand,size = score),colour=color)+
      labs(size = "LR score",y = "Ligand",x="Receptor")+
      theme_bw()+
      FontSize(x.title = 15, y.title = 15,y.text=15,x.text=15)+
      theme(axis.text.x = element_text(angle = 45,hjust = 1))
  }

  if (return_data) {
    return(lrpair)
  }else{
    return(p)
  }
}
