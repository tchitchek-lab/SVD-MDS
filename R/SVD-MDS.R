#' @title Construction of a SVD-MDS representation
#' 
#' @description
#' This function computes a SVD-MDS representation based on a distance matrix. 
#'
#' @param dist a numeric matrix with all pairwise distances between objects of the representation
#' @param k a numeric value specifying the desired number of dimensions in the resulting Reference Map representation
#' @param metric a character indicating the distance metric to use ("euclidean" or "manhattan")
#' @param max_it a numeric defining the maximal number of steps the algorithm can perform
#' @param stress_sd_th a numeric defining the threshold for the standard deviation of Kruskal Stress
#' @param stack_length a numeric defining the length of the Kruskal Stress stack (used to compute the standard deviation of the Kruskal Stress)
#' @param verbose a boolean enabling the display of debug information at each step of the algorithm
#' 
#' @return a list of 3 elements containing the position of the objects ('points' element), the Kruskal Stress ('stress' element), and the Entourage Score ('entourage' element)
#' 
#' @export
SVDMDS = function(dist, k = 2, metric = "euclidean", max_it = 6*10^6, stress_sd_th = 10^-4, stack_length = 500, verbose = TRUE) {
	
	init = "svd"
	
	if (class(dist) != "dist" && is.null(nrow(dist))) {
        stop("dist must be a dist object or a distance matrix")
		}
    if (!is.matrix(init) && !(init %in% c("rand","svd"))) {
        stop("init must be a matrix, a character rand or a character svd")
		}
    if (!(metric %in% c("euclidean","manhattan"))) {
        stop("metric must be a character euclidean or a character manhattan")
		}
    if (k < 1) {
		stop("k must be >= 1")
		}
    if (max_it < 1) {
		stop("max_it must be >= 1")
		}
    if (stress_sd_th < 0) {
		stop("stress_sd_th must be >= 0")
		}
    if (stack_length < 2) {
		stop("stack_length must be >= 2")
		}
    if (!(verbose %in% c(TRUE,FALSE))) {
		stop("verbose must be a logical")
		}
    
    dist = as.matrix(dist)
    if (nrow(dist) != ncol(dist)) {
		stop("dist must be a dist object or a distance matrix")
		}
		
    manhattan = if (metric == "euclidean") {0} else {1}
    if (is.matrix(init)) {
        if (verbose) {print("Initialization with the given init matrix")}
        positions = init    
    } else if (init == "rand") {
        if (verbose) {print("Initialization with random positions")}
        np        = dim(dist)[1]
        positions = matrix(stats::runif(np * k, min = -1, max = 1), nrow = np, ncol = k)
    } else if (init == "svd") {
        if (verbose) {print("Initialization with a SVD on the dist matrix, may takes times")}
        svd       = svd(dist)
        positions = t(svd$u)[,c(1:k)]
    }
    
    stress_best    = Inf
    ref_point      = 0
    positions_best = matrix(positions, ncol = ncol(positions))
    if (verbose) {print("Computation of the SVD-MDS")}
    if (any(is.na(dist))) {
        print("NA value(s) detected in the input distance matrix, computations will be slower")
        invisible(.Call('SVDMDS_NA_core', PACKAGE = 'SVD-MDS', positions_best, stress_best, positions, dist, manhattan, ref_point, max_it, stack_length, stress_sd_th, verbose))
    }else {
        invisible(.Call('SVDMDS_core', PACKAGE = 'SVD-MDS', positions_best, stress_best, positions, dist, manhattan, ref_point, max_it, stack_length, stress_sd_th, verbose))
    }
	
    res           = c()
    res$points    = positions_best
    res$stress    = stress_best
    dist_res      = if (manhattan) {distManhattan(positions_best)} else {distEuclidean(positions_best)}
    res$entourage = computeEntourageScore(dist, dist_res)
    
    return(res)
}
