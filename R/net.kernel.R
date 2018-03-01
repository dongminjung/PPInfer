# kernel matrix
net.kernel <- function (g, decay = 0.5) 
{
  L <- as.matrix(laplacian_matrix(g, normalized = TRUE))
  colname <- colnames(L)
  rowname <- rownames(L)
  I <- diag(1, dim(L))
  K <- chol2inv(chol(I + decay * L))
  colnames(K) <- colname
  rownames(K) <- rowname
  return(K)
}