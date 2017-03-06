#' function to automatically learn the structure of data by either using L1-graph or the spanning-tree formulization
#' @param X the input data DxN
#' @param y the initial cluster assignment
#' @param C0 the initialization of centroids
#' @param G graph matrix with side information where cannot-link pair is 0
#' @param maxiter maximum number of iteraction
#' @param eps relative objective difference
#' @param gstruct graph structure to learn, either L1-graph or the span-tree
#' @param lambda regularization parameter for inverse graph embedding
#' @param gamma regularization parameter for k-means (the prefix of 'param' is used to avoid name collision with gamma)
#' @param sigma bandwidth parameter
#' @param nn number of nearest neighbors
#' @param verbose emit results from iteraction
#' @return a list of X, C, W, P, objs
#' X is the input data
#' C is the centers for principal graph
#' W is the pricipal graph matrix
#' P is the cluster assignment matrix
#' objs is the objective value for the function
#' @export
principal_graph_large <- function(X, y,
	maxiter = 10, eps = 1e-5,
	gstruct = c('l1-graph', 'span-tree'),
	lambda = 1,
	gamma = 0.5,
	sigma = 0.01,
	nn = 5,
	ncenter = NULL,
	verbose = T) { # [y_center, C, W, P, objs] =

	D <- nrow(X); N <- ncol(X)

	if(verbose){
		message("running k-means clustering")
	}

	K <- ncenter

	set.seed(20170126)
	kmean_res <- kmeans(t(X), K, nstart = 10) #Hartigan-Wong algorithm
	PI <- kmean_res$cluster
	centers <- kmean_res$centers

	y_set <- list(ncenter)
	for(i in 1:N)
		y_set[[PI[i]]] <- c( y_set[PI[i]][[1]], y[i])

	y_center <- t(rep(0, ncenter))
for(i in 1:length(y_center)) {
		tbl <- table(y_set[[1]])
		max_val <- max(tbl); max_idx <- which.max(tbl)
	}

	C0 <- t(centers)
	G <- c()

	if(gstruct == 'l1-graph') {
		C0 <- t(centers)
		nC0 <- ncol(C0)

		if(nn < nC0) {
			G_list <- get_knn(C0, nn)
		}
		else
			G <- matrix(1, nrow = nC0, ncol = nC0) - eye(nC0,nC0)
	}

	res <- principal_graph(X, C0, G_list$G,
	    	maxiter = maxiter, eps = eps,
	    	gstruct = gstruct,
	    	lambda = lambda,
	    	gamma = gamma,
	    	sigma = sigma,
	    	nn = nn,
	    	verbose = verbose)

	return(res)
}

