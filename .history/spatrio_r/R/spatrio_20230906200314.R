#' Spatrio analysis using optimal transport
#'
#' @param input_path Path to load input data
#' @param output_path Path to save output data
#' @param py_path Path to python interpreter
#' @param spatrio_path Path to spatrio python package
#' @param marker_use Determines whether to select differential genes of each cell/spot type for subsequent analysis
#' @param top_marker_num Number of differential genes in each cell/spot type
#' @param hvg_use Determines whether to only reserve highly variable genes
#' @param alpha Alignment tuning parameter. Note: 0 <= alpha <= 1
#' @param dissimilarity Expression dissimilarity measure
#' @param k Number of neighbors to be used when constructing kNN graph
#' @param graph_mode Determines whether to use a connectivity graph or a distance graph
#' @param aware_spatial Determines whether to adjust the distance between spots according to areas (or other meta info)
#' @param aware_multi Determines whether to adjust the distance between cells according to types (or other meta info)
#' @param aware_power Type aware parameter. The greater the parameter, the greater the distance between different areas/types of spots/cells in the graph
#' @param The maximum number of cells allocated in a spot
#' @param random Determines whether to randomly assign cell coordinates or assign coordinates based on pearson correlation coefficient
#' @param verbose
#'
#' @export
#'

spatrio <- function(
    input_path = getwd(),
    ref_path = getwd(),
    output_path = getwd(),
    py_path = NULL,
    spatrio_path = NULL,
    marker_use = TRUE,
    top_marker_num = 100,
    hvg_use = FALSE,
    alpha = 0.1,
    dissimilarity = "scaled_euc",
    k = 10,
    graph_mode = "connectivity",
    aware_spatial = TRUE,
    aware_multi = TRUE,
    aware_power = 2,
    top_num = 5,
    random = FALSE,
    verbose = TRUE,
    ...
) {
  if (!dir.exists(paths = input_path)) {
    stop("Requested input directory does not exist")
  }
  if (!dir.exists(paths = output_path)) {
    stop("Requested output directory does not exist")
  }
  if (!dir.exists(paths = spatrio_path)) {
    stop("Requested SpaTrio directory does not exist")
  }
  # find python
  if (!file.exists("/home/yph/anaconda3/envs/spatrio/bin/python")) {
    stop("Please enter the correct python path")
  }else{
    print(paste("Using the Python interpreter with a path of ",py_path,sep = ""))
  }

  cmd <- paste0(
    py_path," ",
    spatrio_path,'/run.py',
    ' --input_path ',input_path,
    ' --output_path ',output_path,
    ' --marker_use ',marker_use,
    ' --top_marker_num ',top_marker_num,
    ' --hvg_use ',hvg_use,
    ' --alpha ',alpha,
    ' --dissimilarity ',dissimilarity,
    ' --k ',k,
    ' --graph_mode ',graph_mode,
    ' --aware_spatial ',aware_spatial,
    ' --aware_multi ',aware_multi,
    ' --aware_power ',aware_power,
    ' --top_num ',top_num,
    ' --random ',random
  )
  if (!is.null(ref_path)) {
  cmd <- paste(cmd, ' --ref_path ', ref_path)
}

  # call SpaTrio
  system(
    command = cmd,
    wait = TRUE,
    ignore.stderr = !verbose,
    ignore.stdout = !verbose
  )
}
