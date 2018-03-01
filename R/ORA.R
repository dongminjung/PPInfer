# Over-representation Analysis
ORA <- function(pathways, gene.id, minSize = 1, maxSize = Inf,
                p.adjust.methods = NULL)
{
  K <- length(gene.id)
  M <- lengths(pathways)
  index <- (minSize < M) & (M < maxSize)
  M <- M[index]
  index <- seq(1,length(pathways))[index]
  m <- length(index)
  geneNames <- unique(unlist(pathways))
  N <- length(geneNames)
  
  pvalue <- 0
  x <- 0
  pb <- txtProgressBar(min = 0, max = m, style = 3)
  for(i in 1:m)
  {
    x[i] <- length(intersect(gene.id, pathways[[index[i]]]))
    ct <- matrix(c(x[i], K-x[i], M[i]-x[i], N-M[i]-(K-x[i])), 2, 2)
    pvalue[i] <- fisher.test(ct, alternative = 'greater')$p.value
    Sys.sleep(0.01)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  padj <- p.adjust(pvalue, method = p.adjust.methods)
  Count <- x
  Size <- M
  Category <- names(M)
  data.frame(Category, Count, Size, pvalue, padj, row.names = NULL)
}