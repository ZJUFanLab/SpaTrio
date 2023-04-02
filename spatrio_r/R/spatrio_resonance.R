#' Spatially weighted correlation
#'
#' @param mat Expression matrix, n X p
#' @param weighted_mat Weight matrix, n X n
#' @param method Correlation method, pearson or spearman
#' @param na_zero Na to zero
#'
#' @return Weighted correlation matrix, p X p
#' @export
#'
#' @import magrittr
weighted_cor <- function(mat, weighted_mat, method='pearson', na_zero=T) {
  if (method=='spearman') mat <- apply(mat, 2, rank)
  mat <- scale(mat)
  mat[is.nan(mat)] <- NA
  cov <- (t(mat) %*% weighted_mat %*% mat)
  diag <- sqrt(diag(cov) %*% t(diag(cov)))
  cor <- cov/diag
  if (na_zero) cor[which(is.na(cor), arr.ind=T)] <- 0
  cor
}

#' Filter the genes with low correlation
#'
#' @param cor_mat
#' @param ave_cor_cut
#' @param min_n
#' @param max_n
#' @param na_diag
#'
#' @return
cor_filtering <- function (mat, min_cutoff = 0.5, min=5, max=100, na_diag=F) {
  if (na_diag) {
    diag(mat) <- NA
  }
  mat_temp <- mat
  cor_ave <- rowMeans(mat_temp, na.rm = T)
  cor_max <- max(mat_temp, na.rm = T)
  if ((cor_max < min_cutoff)|(nrow(mat_temp)<min)) {
    mat_temp <- matrix(NA, 1, 1)
  } else {
    cor_min <- min(cor_ave)
    cor_min_idx <- which.min(cor_ave)
    idx <- 1
    while ((cor_min < min_cutoff & idx <= (nrow(mat)-min))|((nrow(mat)-idx))>=max) {
      mat_temp <- mat_temp[-cor_min_idx, -cor_min_idx]
      cor_ave <- rowMeans(mat_temp, na.rm = T)
      cor_min <- min(cor_ave)
      cor_min_idx <- which.min(cor_ave)
      idx <- idx + 1
    }
  }
  return (mat_temp)
}

#' Smooth feature matrix using kNNs
#'
#' Method to smooth sparse scores (e.g. RNA expression, ADT expression, or motif scores) per cell per feature using cell K nearest neighbors (NNs).
#'
#'@param nn_matrix Matrix of K Nearest neighbor indices for each cell (row), where each column is one of K nearest neighbors per cell.
#'@param mat Matrix of gene/adt expression or motif activity score. Column names must match the rownames of nn_matrix.
#'@param nCores Integer specifying the number of cores to use. Default is 1.
#'@importFrom parallel mclapply
#'@import Matrix doParallel dplyr
#'@return A matrix of scores for each feature and each cell, smoothed over its K nearest neighbors.
#'@export
#'
#'@return
smoothScoresNN <- function(nn_matrix,
                           mat,
                           geneList = NULL,
                           subset=NULL,
                           nCores = 1)
{
  library(doParallel)
  library(FNN)
  stopifnot(all.equal(nrow(nn_matrix),ncol(mat)))

  if (is.null(rownames(nn_matrix)))
    stop("nn_matrix has to have matching cell id as rownames\n")
  if (!all.equal(rownames(nn_matrix), colnames(mat)))
    stop("Nearest-neighbor matrix and cell data matrix don't have matching cells barcodes ..\n")
  opts <- list()
  pb <- txtProgressBar(min = 0, max = ncol(mat), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  time_elapsed <- Sys.time()

  cl <- parallel::makeCluster(nCores)
  doSNOW::registerDoSNOW(cl)

  cat("Running in parallel using ", nCores, "cores ..\n")
  matL <- foreach::foreach(x=1:nrow(nn_matrix),.options.snow = opts,.packages = c("Matrix","data.table","dplyr")) %dopar% {
    smoothedScore <- data.table(Matrix::rowMeans(mat[, nn_matrix[x,]]))
    rownames(smoothedScore) <- rownames(mat)
    colnames(smoothedScore) <- rownames(nn_matrix)[x]
    smoothedScore
  }

  parallel::stopCluster(cl)

  close(pb)
  cat("Merging results ..\n")
  smoothedMat <- dplyr::bind_cols(matL) %>% data.matrix() %>% Matrix::Matrix(sparse=TRUE)
  rownames(smoothedMat) <- rownames(mat)

  time_elapsed <- Sys.time() - time_elapsed
  cat(paste("\nTime Elapsed: ", time_elapsed, units(time_elapsed),
            "\n"))

  return(smoothedMat)
}


#' Space resonance feature module
#'
#' @param object Input data in SeuratObject format
#' @param sigm
#' @param assay Which assay to use
#' @param feature_select
#' @param correlation
#' @param maxK
#' @param k
#' @param min_avg_con
#' @param min_avg_cor
#' @param min_featuer
#' @param max_featuer
#' @param min_pct_cutoff Filter features with low expression ratio
#' @param smooth
#' @param smooth_k
#' @param reduction Which low-dimensional representation to use for smoothing
#' @param loc_mat Coord matrix
#' @param scale
#' @param ...
#'
#' @return
#' @export
#'
#' @import Seurat
#' @import magrittr
#' @import FNN
#'
spatrio_resonance <- function(object,
                              sigma=NULL,
                              assay='RNA',
                              feature_select=NULL,
                              correlation='pearson',
                              maxK=8,
                              k=8,
                              min_avg_con=0.5,
                              min_avg_cor=0.5,
                              min_featuer=20,
                              max_featuer=100,
                              min_pct_cutoff=NULL,
                              smooth=NULL,
                              smooth_k=10,
                              reduction=NULL,
                              scale=T,
                              loc_mat=NULL){
  cat('Feature filtering...\n')
  if(!is.null(min_pct_cutoff)){
    cat('Filtering features...\n')
    feature_nz <- apply(object[[assay]]@data[feature_select,], 1, function(x) sum(x>0))
    features <- names(feature_nz)[feature_nz > ncol(object)*min_pct_cutoff]
    feature_select<-intersect(feature_select,features)
  }else{
    cat('Using all features...\n')
    feature_select<-rownames(object[[assay]]@data)
  }

  cat('Data smoothing\n')
  if(smooth=="reduction"){
    cat('Smoothing with reduction...\n')
    featureMat<-object[[assay]]@data
    cellKNN.mat<-object@reductions[[reduction]]@cell.embeddings
    cellkNN <- FNN::get.knn(cellKNN.mat,k = smooth_k)$nn.index
    rownames(cellkNN)<-rownames(cellKNN.mat)
    featureMat.smooth <- smoothScoresNN(nn_matrix=cellkNN,mat=featureMat,nCores=5)
    object[[assay]]@data<-as.sparse(featureMat.smooth)
  }
  if(smooth=="location"){
    cat('Smoothing with location...\n')
    featureMat<-object[[assay]]@data
    cellKNN.mat<-loc_mat
    cellkNN <- FNN::get.knn(cellKNN.mat,k = smooth_k)$nn.index
    rownames(cellkNN)<-rownames(cellKNN.mat)
    featureMat.smooth <- smoothScoresNN(nn_matrix=cellkNN,mat=featureMat,nCores=5)
    object[[assay]]@data<-as.sparse(featureMat.smooth)
  }
  if(smooth=="combined"){
    cat('Smoothing with reduction and location...\n')
    featureMat<-object[[assay]]@data
    cellKNN.mat<-object@reductions[[reduction]]@cell.embeddings
    cellkNN <- FNN::get.knn(cellKNN.mat,k = smooth_k)$nn.index
    rownames(cellkNN)<-rownames(cellKNN.mat)
    featureMat.smooth1 <- smoothScoresNN(nn_matrix=cellkNN,mat=featureMat,nCores=5)

    featureMat<-object[[assay]]@data
    cellKNN.mat<-loc_mat
    cellkNN <- FNN::get.knn(cellKNN.mat,k = smooth_k)$nn.index
    rownames(cellkNN)<-rownames(cellKNN.mat)
    featureMat.smooth2 <- smoothScoresNN(nn_matrix=cellkNN,mat=featureMat,nCores=5)

    featureMat.smooth<-featureMat.smooth1+featureMat.smooth2
    object[[assay]]@data<-as.sparse(featureMat.smooth/2)
  }

  dist_mat <- dist(object@meta.data[, c('xcoord', 'ycoord')]) %>% as.matrix
  kern_mat <- exp(-1*(dist_mat^2)/(2*sigma^2))
  DefaultAssay(object)<-assay
  cat("Calculating with",length(feature_select),"features...\n")
  if(scale){
    cat("Data scaling...\n")
    object<-ScaleData(object,features =feature_select )
    data_mat <- t(as.matrix(object[[assay]]@scale.data)[feature_select,])
  } else {
    data_mat <- t(as.matrix(object[[assay]]@data)[feature_select,])
  }

  cat('Spatially weighted correlation calculating...\n')
  wcor_mat <- weighted_cor(mat=data_mat, weighted_mat=kern_mat, method=correlation)

  wcor_dis <- as.dist(1-wcor_mat)
  cat('Consensus clustering...\n')
  c_cluster_res <- ConsensusClusterPlus::ConsensusClusterPlus(wcor_dis, maxK=maxK,
                                                              reps=20, distance="spearman", verbose=F,plot='png')
  cons_mat <- c_cluster_res[[k]]$consensusMatrix %>% data.frame
  colnames(cons_mat) <- rownames(wcor_mat)
  colnames(cons_mat) <- rownames(cons_mat)

  cat('Feature module calculating...\n')
  c_cluster_k <- c_cluster_res[[k]]
  con_k <- c_cluster_k$consensusMatrix %>%
    magrittr::set_rownames(names(c_cluster_k$consensusClass)) %>%
    magrittr::set_colnames(names(c_cluster_k$consensusClass))
  c_cluster_module <- list()
  for (i in 1:k) {
    features <- names(c_cluster_k$consensusClass)[c_cluster_k$consensusClass==i]
    cor_temp <- wcor_mat[features, features, drop=F]
    cor_temp <- cor_filtering(cor_temp, max=max_featuer, min=1, min_cutoff=min_avg_cor)
    features <- rownames(cor_temp)
    if (is.null(features)) features <- NA
    con_temp <- con_k[features, features, drop=F]
    con_temp <- cor_filtering(con_temp, max=max_featuer, min=1, min_cutoff=min_avg_con)
    features <- rownames(con_temp)
    if (is.null(features)) features <- NA
    c_cluster_module[[i]] <- features
  }
  names(c_cluster_module) <- paste0('k', 1:k)
  print(min_featuer)
  c_cluster_module <- c_cluster_module[sapply(c_cluster_module, function(x) length(x)>=min_featuer)]


  output <- list(group=c(), c_cluster=c(), rbf_kernel=c(), weighted_cor=c())
  output$group <- c_cluster_module
  output$c_cluster <- c_cluster_res
  output$rbf_kernel <- kern_mat
  output$weighted_cor <- wcor_mat

  return(output)

}


