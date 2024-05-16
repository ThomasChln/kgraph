
###############################################################################
#' sparse_encode
#'
#' Sparse encoding method by closest neighbors. 
#' Three methods are available:
#'   - hard encoding: each patient's closest neighbors are set to 1, others
#'     are set to 0
#'   - soft encoding: each patient's closest neighbors distances are
#'     transformed by the exponential norm, others are set to 0
#'   - epsilon encoding: each patient's neighbors closer than the mean of
#'     the distance matrix are transformed by the exponential norm and
#'     others are set to 0.
#'
#' @param m_data      Numeric matrix
#' @param dist_method Distance method passed to qb_dist
#' @param encoding    Encoding method: one of hard, soft, or epsilon
#' @param sigma       Parameter for the exponential norm transform.
#'                    Default is mean of std. dev. of distance matrix columns
#' @param n_neighbors Number of neighbors (ignored in epsilon encoding)
#' @param scale_obs   Scale by observations
#' @return Projected matrix
#' @export
sparse_encode <- function(m_data, dist_method = 'norm_inprod',
  encoding = c('epsilon', 'hard', 'soft'), sigma,
  n_neighbors = floor(nrow(m_data) / 10), scale_obs = TRUE) {

  if (scale_obs) m_data <- t(scale(t(m_data)))
  m_data <- as.matrix(m_data)
  stopifnot(is.matrix(m_data) && is.numeric(m_data))

  encoding <- match.arg(encoding)
  n_neighbors <- as.numeric(n_neighbors)
  if(!is.na(n_neighbors) && (n_neighbors < 1 || n_neighbors > nrow(m_data))) {
    stop('n_neighbors out of bounds')
  }

  m_dist <- dist_matrix(m_data, dist_method)

  if (missing(sigma)) sigma <- stddev_mean(m_dist)
  if (sigma <= 0) stop('sigma must be strictly positive')

  sparse_encode_dist(m_dist, dist_method, encoding, sigma, n_neighbors)
}

# sparse encode a distance matrix
# set diagonal to NA and iterate on columns (distance matrix)
# for epsilon encoding, apply norm_inprod and substract by mean
sparse_encode_dist <- function(m_dist, dist_method, encoding, ...) {

  diag(m_dist) <- NA
  sp_enc <- apply(m_dist, 2, sparse_encode_vec, encoding, ...)

  if (encoding == 'epsilon') {
    sp_enc <- norm_inprod(sp_enc)
    sp_enc <- matrix(pmax(sp_enc - mean(sp_enc), 0), nrow(m_dist),
      dimnames = list(rownames(m_dist), rownames(m_dist)))
  }

  sp_enc
}

# sparse encode a vector
# set farthest neighbors to 0 and apply exp_transform for epsilon and soft
sparse_encode_vec <- function(vec, encoding, sigma, n_neighbors) {
  switch(encoding, epsilon = {
      vec[which(!is.na(vec))] <- exp_transform(na.omit(vec), sigma)
      vec[which(is.na(vec))] <- 0
      vec
    }, {
      neighbors <- utils::head(order(vec), n_neighbors)
      vec[-neighbors] <- 0
      vec[neighbors] <- switch(encoding, hard = 1,
        soft = exp_transform(vec[neighbors], sigma))
      vec
    })
}

# normed exponential transform
exp_transform <- function(vec, sigma) {
  vec <- exp( - vec / sigma)
  vec / sum(vec)
}

##############################################################################
#' dist_matrix
#'
#' Dispatch of amap::Dist, lsa::cosine, and norm_inprod methods.
#'
#' @param data    Rectangular numeric matrix [Observations, Features]
#' @param method  Methods accepted by amap::Dist or cosine
#' @param n_cores Number of cores
#' @return Dissimarility symmetric matrix
#'
#' @export
dist_matrix <- function(data, method = 'euclidean', n_cores = 1) {

  amap_methods <- c("euclidean", "maximum", "manhattan", "canberra",
    "binary", "pearson", "correlation", "spearman", "kendall",
    "abspearson", "abscorrelation")
  methods <- c(amap_methods, 'cosine', 'norm_inprod')

  if (!method %in% methods) stop('Invalid distance method')

  switch(method, cosine = 1 - lsa::cosine(t(data)),
    norm_inprod = -norm_inprod(data), {
      as.matrix(amap::Dist(data, method, nbproc = n_cores))
    })
}

###############################################################################
#' stddev_mean
#'
#' Get mean of standard deviations of matrix columns
#'
#' @param m Numeric matrix
#' @return Mean of standard deviations of matrix columns
#'
#' @export
stddev_mean <- function(m) {
  m <- as.matrix(m)
  stopifnot(is.matrix(m) && is.numeric(m))
  mean(sqrt(apply(m, 2, var)))
}


###############################################################################
#' norm_inprod 
#'
#' Normalized inner product with transposed input matrix
#'
#' @param m Numeric matrix
#' @return Numeric matrix
#'
#' @export
norm_inprod <- function(m) {
  m_inprod <- inprod(m - rowMeans(m))
  m_inprod <- apply(m_inprod, 1, function(i) i / as.vector(sqrt(i %*% i)))
  rownames(m_inprod) <- rownames(m)

  m_inprod
}

# inner product
inprod <- function(m) m %*% t(m)

