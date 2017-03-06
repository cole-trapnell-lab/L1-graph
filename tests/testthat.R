library(testthat)
library(L1Graph)
library(ggplot2)
library(Matrix)
library(monocle)
library(destiny)
library(R.matlab)
library(simplePPT)

test_check("L1Graph")

run_our_method <- function(filename, maxiter = 20,
                           eps = 1e-5,
                           gstruct = "l1-graph",
                           lambda = 1.0,
                           gamma = 0.5,
                           sigma = 0.01,
                           nn = 5,
                           verbose = TRUE) {

  if(filename == 'Circle'){
    X <- readMat(paste(data_dir, filename, '.mat', sep =''))$X
    X <- t(X)
    sigma <- 0.1; nn = 10
  }
  else if(filename == 'two_moon') {
    X <- readMat(paste(data_dir, filename, '.mat', sep =''))$X
    X <- t(X)
    lambda <- 0.1; gamma = 3
  }
  else if(filename == 'tree_300') {
    X <- readMat(paste(data_dir, filename, '.mat', sep =''))$X
    X <- t(X)
    gamma = 10
  }
  else if(filename == 'Spiral_SUN') {
    X <- readMat(paste(data_dir, filename, '.mat', sep =''))$X
    X <- t(X)
    nn = 10
  }
  else if(filename == 'three_clusters') {
    X <- readMat(paste(data_dir, filename, '.mat', sep =''))$X
    X <- t(X)
    lambda = 0.1
  }
  else if(filename == 'DistortedSShape') {
    X <- readMat(paste(data_dir, filename, '.mat', sep =''))$X
    X <- t(X)
    lambda = 0.1
  }
  else
    stop('unexpected data settings for ', name)

  D <- nrow(X); N <- ncol(X)
  Z <- X
  C0 <- Z
  Nz <- ncol(C0)

  if(nn < N)
    G <- get_knn(C0, nn)
  else
    G <- matrix(1, nrow = Nz, ncol = Nz) - eye(Nz,Nz);

  # ptm <- proc.time()
  res <- principal_graph(X, C0, G$G, maxiter = maxiter,
                         eps = eps, gstruct = gstruct,
                         lambda = lambda, gamma = gamma,
                         sigma = sigma, nn = 5, verbose = verbose)
  # proc.time() - ptm
  return(list(C = res$C, W = res$W, P = res$P, objs = res$objs, X = X))
}

###################################################################################################################################################################################
# test this on  all the examples from Qi Mao
###################################################################################################################################################################################
return_myself <- function(data) {
  return(t(data))
}

create_cds <- function(data) {
  pd <- new("AnnotatedDataFrame", data = data.frame(cell = 1:ncol(data), row.names = paste('cell', 1:ncol(data), sep = '')))
  fd <- new("AnnotatedDataFrame", data = data.frame(gene = 1:nrow(data), row.names = paste('gene', 1:nrow(data), sep = '')))

  data <- as.data.frame(data)
  dimnames(data) <- list(paste('gene', 1:nrow(data), sep = ''), paste('cell', 1:ncol(data), sep = ''))
  new_cds <- newCellDataSet(as(as.matrix(data), "sparseMatrix"),
                            phenoData = pd,
                            featureData = fd,
                            expressionFamily=gaussianff(),
                            lowerDetectionLimit=1)

  new_cds <- estimateSizeFactors(new_cds)
  # new_cds <- estimateDispersions(new_cds)
  new_cds
}

maxiter = 20
eps = 1e-5
gstruct = "l1-graph"
lambda = 1.0
gamma = 0.5
sigma = 0.01
nn = 5

filenames = c('Circle','two_moon','tree_300','Spiral_SUN','three_clusters','DistortedSShape')
data_dir <- '~/Dropbox (Personal)/Projects/SPL+L1graph/tree-l1-code/toy/'
res_list <- list()
data_len <- rep(0, length(filenames))
for(i  in 1:length(filenames)){
  print(filenames[i])
  res1 <- run_our_method(filenames[i])
  new_cds <- create_cds(data = res_list[[1]]$X)


  DistortedSShape <- reduceDimension(Circle, norm_method = 'none', pseudo_expr = 0, scaling = F,
                                     reduction_method = 'L1-graph', maxiter = maxiter, initial_method = return_myself,
                                     eps = eps, lambda = 1, gamma = gamma,
                                     sigma = 0.1, nn = 10, verbose = T)
  DistortedSShape <- orderCells(DistortedSShape)
  plot_cell_trajectory(DistortedSShape, color_by = 'cell')

  # print(qplot(res1$C[1, ], res1$C[2, ], color = 'red') + geom_point(aes(res1$X[1, ], res1$X[2, ]), color = 'black') + ggtitle('l1-graph'))

  res2 <- run_our_method(filenames[i], gstruct = 'span-tree')
  # print(qplot(res2$C[1, ], res2$C[2, ], color = 'red') + geom_point(aes(res2$X[1, ], res2$X[2, ]), color = 'black') + ggtitle('spanning-tree'))

  print(dim(res1$C))
  res_list <- c(res_list, list(res1, res2))
}

data_len <- unlist(lapply(res_list, function(x) nrow(t(x$C))))
data_len2 <- unlist(lapply(res_list, function(x) nrow(t(x$X))))
C_data <- do.call(rbind.data.frame, lapply(res_list, function(x) t(x$C)))
C_data$name <- c(rep(filenames[6], data_len[1] * 2), rep(filenames[5], data_len[3] * 2),
                 rep(filenames[4], data_len[5] * 2), rep(filenames[3], data_len[7] * 2),
                 rep(filenames[2], data_len[9] * 2), rep(filenames[1], data_len[11] * 2))
C_data$method <- c(rep(c('L1-graph', 'L1-span-tree'), each = data_len[1]), rep(c('L1-graph', 'L1-span-tree'), each = data_len[3]),
                 rep(c('L1-graph', 'L1-span-tree'), each = data_len[5]), rep(c('L1-graph', 'L1-span-tree'), each = data_len[7]),
                 rep(c('L1-graph', 'L1-span-tree'), each = data_len[9]), rep(c('L1-graph', 'L1-span-tree'), each = data_len[11]))

X_data <- do.call(rbind.data.frame, lapply(res_list, function(x) t(x$X)))
X_data$name <- c(rep(filenames[6], data_len[1] * 2), rep(filenames[5], data_len[3] * 2),
                 rep(filenames[4], data_len[5] * 2), rep(filenames[3], data_len[7] * 2),
                 rep(filenames[2], data_len[9] * 2), rep(filenames[1], data_len[11] * 2))
X_data$method <- c(rep(c('L1-graph', 'L1-span-tree'), each = data_len[1]), rep(c('L1-graph', 'L1-span-tree'), each = data_len[3]),
                   rep(c('L1-graph', 'L1-span-tree'), each = data_len[5]), rep(c('L1-graph', 'L1-span-tree'), each = data_len[7]),
                   rep(c('L1-graph', 'L1-span-tree'), each = data_len[9]), rep(c('L1-graph', 'L1-span-tree'), each = data_len[11]))

XC_data <- rbind(X_data, C_data)
XC_data$Type <- rep(c('raw', 'principal_point'), each = nrow(X_data))

qplot(V1, V2, color = Type, data = XC_data) + facet_wrap(~name + method, nrow = 6)

XC_data_s1 <- subset(XC_data, name %in% filenames[1:3])
XC_data_s2 <- subset(XC_data, name %in% filenames[4:6])
qplot(V1, V2, color = Type, data = XC_data_s1) + facet_wrap(~name + method, nrow = 3, scales = 'free')
qplot(V1, V2, color = Type, data = XC_data_s2) + facet_wrap(~name + method, nrow = 3, scales = 'free')


#create cds to ensure we can learn the circle and parallel trajectory:
filename = filenames{i};

maxiter = 20
eps = 1e-5
gstruct = "l1-graph"
lambda = 1.0
gamma = 0.5
sigma = 0.01
nn = 5

###################################################################################################################################################################################
# run with monocle 2
###################################################################################################################################################################################
library(devtools)
load_all('/Users/xqiu/Dropbox (Personal)/Projects/monocle-dev')

filenames = c('Circle','two_moon','tree_300','Spiral_SUN','three_clusters','DistortedSShape')
data_dir <- '~/Dropbox (Personal)/Projects/SPL+L1graph/tree-l1-code/toy/'

DistortedSShape <- create_cds(data = res_list[[1]]$X)
return_myself <- function(data) {
  return(t(data))
}

DistortedSShape <- reduceDimension(Circle, norm_method = 'none', pseudo_expr = 0, scaling = F,
                          reduction_method = 'L1-graph', maxiter = maxiter, initial_method = return_myself,
                          eps = eps, lambda = 1, gamma = gamma,
                          sigma = 0.1, nn = 10, verbose = T)
DistortedSShape <- orderCells(DistortedSShape)
plot_cell_trajectory(DistortedSShape, color_by = 'cell')

###################################################################################################################################################################################
# test this on the large example data from Qi's paper:
###################################################################################################################################################################################

experiment_large <- readMat('/Users/xqiu/Dropbox (Personal)/Projects/SPL+L1graph/tree-l1-code/matlab.mat')

newX <- experiment_large$newX
C0 <- experiment_large$C0
G <- experiment_large$G

maxiter <- 20
eps <- 1.0000e-05
gstruct <- 'span-tree'
gamma <- 1 #smoothness
sigma <- 0.01000 #
lambda <- 1 #L1 g
nn <- 10
verbose = T

  experiment_large_res <- principal_graph(t(newX), C0, G, maxiter = maxiter,
                                  eps = eps, gstruct = 'span-tree',
                                  lambda = lambda, gamma = gamma,
                                  sigma = sigma, nn = 5, verbose = T)

print(qplot(experiment_large_res$C[1, ], experiment_large_res$C[2, ], color = 'red') + geom_point(aes(newX[, 1], newX[, 3]), color = 'black') + ggtitle('spanning-tree'))

###################################################################################################################################################################################
# run principal tree (simplePPT)
###################################################################################################################################################################################
load('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/RData/analysis_score_ordering_MAR_seq.RData')

dpt_res <- run_new_dpt(valid_subset_GSE72857_cds[, ], normalize = F)
dm_res <- DM(log2(exprs(valid_subset_GSE72857_cds) + 1))
mar_seq_pt_res <- principal_tree(t(dpt_res$dm@eigenvectors[, 1:2]), MU = NULL, lambda = params.lambda, bandwidth = params.bandwidth, maxIter = 10, verbose = T)

###################################################################################################################################################################################
# run principal graph (L1 graph)
###################################################################################################################################################################################
X <- t(dm_res)
D <- nrow(X); N <- ncol(X)
Z <- X
C0 <- Z
Nz <- ncol(C0)

G <- get_knn(C0, nn)

# choose appropriate lambda, gamma and sigma

# arbitray
mar_seq_pg_res <- principal_graph(t(dm_res), C0, G$G, maxiter = maxiter,
                                  eps = eps, gstruct = 'l1-graph',
                                  lambda = lambda, gamma = gamma,
                                  sigma = sigma, nn = 5, verbose = T)

print(qplot(mar_seq_pg_res$C[1, ], mar_seq_pg_res$C[2, ], color = 'red') + geom_point(aes(X[1, ], X[2, ]), color = 'black') + ggtitle('L1-graph'))

# spanning tree
mar_seq_pg_res_span_tree <- principal_graph(t(dm_res), C0, G$G, maxiter = maxiter,
                                  eps = eps, gstruct = 'span-tree',
                                  lambda = lambda, gamma = gamma,
                                  sigma = sigma, nn = 5, verbose = T)

print(qplot(mar_seq_pg_res_span_tree$C[1, ], mar_seq_pg_res_span_tree$C[2, ], color = 'red') + geom_point(aes(X[1, ], X[2, ]), color = 'black') + ggtitle('spanning-tree'))

###################################################################################################################################################################################
# run principal graph (L1 graph) for the muscle dataset
###################################################################################################################################################################################
load('/Users/xqiu/Dropbox (Personal)/Projects/DDRTree_fstree/DDRTree_fstree/RData/fig_si2.RData')

pca_res <- PCA(as.matrix(log(exprs(HSMM_myo[HSMM_myo_ordering_genes, ]) + 1)))
qplot(pca_res[, 1], pca_res[, 2], color = pData(HSMM_myo)$Time)
HSMM_pca_res <- principal_tree(t(pca_res[, 1:2]), MU = NULL, lambda = params.lambda, bandwidth = params.bandwidth, maxIter = 10, verbose = T)
pca_HSMM_myo <- reduceDimension(HSMM_myo, reduction_method = 'simplePPT', initial_method = PCA, verbose = T)
pca_HSMM_myo <- orderCells(pca_HSMM_myo)

maxiter <- 20
eps <- 1.0000e-05
gstruct <- 'span-tree'
gamma <- 0.1 #smoothness
sigma <- 0.01000 #
lambda <- 1 #L1 g
nn <- 5
verbose = T

X <- t(dm_res[, 1:2])
# D <- nrow(X); N <- ncol(X)
# Z <- X
C0 <- X
Nz <- ncol(C0)

G <- get_knn(C0, nn)

# choose appropriate lambda, gamma and sigma

# spanning tree
HSMM_seq_pg_res_span_tree <- principal_graph(t(dm_res[, 1:2]), C0, G$G, maxiter = maxiter,
                                            eps = eps, gstruct = 'span-tree',
                                            lambda = lambda, gamma = gamma,
                                            sigma = sigma, nn = 5, verbose = T)

print(qplot(HSMM_seq_pg_res_span_tree$C[1, ], HSMM_seq_pg_res_span_tree$C[2, ], color = I('black')) + geom_point(aes(dm_res[, 1], dm_res[, 2], color = as.character(pData(HSMM_myo)$Time))) + ggtitle('spanning-tree'))

HSMM_myo_l1_span_tree <- reduceDimension(HSMM_myo, reduction_method = 'L1-span-tree', maxiter = maxiter, initial_method = DM,
                        eps = eps, lambda = lambda, gamma = gamma,
                        sigma = sigma, nn = 5, verbose = T)
HSMM_myo_l1_span_tree <- orderCells(HSMM_myo_l1_span_tree)
pdf(paste(SI_fig_dir, 'dm_l1_span_tree.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(HSMM_myo_l1_span_tree, color_by = 'Time', cell_size = 0.2, show_branch_points = F) + nm_theme()
dev.off()

HSMM_myo_l1_graph <- reduceDimension(HSMM_myo, reduction_method = 'L1-graph', C0 = HSMM_seq_pg_res_span_tree$C, maxiter = 2, initial_method = DM,
                            eps = eps,
                            lambda = 1, gamma = 0.5,
                            sigma = 0.015, nn = 5, verbose = T)
HSMM_myo_l1_graph <- orderCells(HSMM_myo_l1_graph)

pdf(paste(SI_fig_dir, 'dm_l1_graph.pdf', sep = ''), height = 1, width = 1)
plot_cell_trajectory(HSMM_myo_l1_graph, color_by = 'Time', cell_size = 0.2, show_branch_points = F) + nm_theme()
dev.off()

# arbitray
G <- get_knn(HSMM_seq_pg_res_span_tree$C, nn)
HSMM_seq_pg_res <- principal_graph(t(dm_res[, 1:2]), HSMM_seq_pg_res_span_tree$C, G$G, maxiter = 2,
                                   eps = eps, gstruct = 'l1-graph',
                                   lambda = 1, gamma = 0.5,
                                   sigma = 0.01, nn = 5, verbose = T)

print(qplot(HSMM_seq_pg_res$C[1, ], HSMM_seq_pg_res$C[2, ], color = I('black')) + geom_point(aes(dm_res[, 1], dm_res[, 2], color = as.character(pData(HSMM_myo)$Time))) + ggtitle('l1-graph'))

