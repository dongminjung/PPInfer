# enrichment network
enrich.net <- function(x, gene.set, node.id, node.name = node.id, pvalue,
                       n = 50, numChar = NULL, pvalue.cutoff = 0.05,
                       edge.cutoff = 0.05, degree.cutoff = 0,
                       edge.width = function(x) {10*x^2},
                       node.size = function(x) {2.5*log10(x)},
                       group = FALSE, group.color = c('red', 'green'),
                       group.shape = c('circle', 'square'),
                       legend.parameter = list('topright'),
                       show.legend = TRUE, ...)
{
  x <- data.frame(x, group)
  colnames(x)[length(colnames(x))] <- 'Group'
  x <- x[x[,pvalue] < pvalue.cutoff, ]
  x <- x[order(x[,pvalue]),]
  n <- min(nrow(x),n)
  if (n == 0)
  {
    stop("no enriched term found...")
  }
  x <- x[1:n,]
    
  index <- match(x[, node.id], names(gene.set))
  geneSets <- list()
  for(i in 1:n)
  {
    geneSets[[i]] <- gene.set[[index[i]]]
  }
  names(geneSets) <- x[,node.name]
    
  if (is.null(numChar))
  {
    numChar <- max(nchar(as.character(x[,node.name])))
  }
  else
  {
    if (length(unique(substr(x[,node.name], 1, numChar))) < nrow(x))
    {
      numChar <- max(nchar(as.character(x[,node.name])))
      message("Note : numChar is too small.", "\n")
    }
  }
  x[,node.name] <- paste(substr(x[,node.name], 1, numChar), 
                          ifelse(nchar(as.character(x[,node.name])) > numChar, 
                                 "...", ""), sep = "")
    
  w <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:n)
  {
    for (j in i:n)
    {
      u <- unlist(geneSets[i])
      v <- unlist(geneSets[j])
      w[i, j] = length(intersect(u, v))/length(unique(c(u, v)))
    }
  }
    
  list.edges <- stack(data.frame(w))
  list.edges <- cbind(list.edges[,1], rep(x[,node.name], n),
              rep(x[,node.name], each = n))
  list.edges <- list.edges[list.edges[, 2] != list.edges[, 3], ]
  list.edges <- list.edges[!is.na(list.edges[, 1]), ]
  g <- graph.data.frame(list.edges[, -1], directed = FALSE)
    
  E(g)$width = edge.width(as.numeric(list.edges[, 1]))
  V(g)$size <- node.size(lengths(geneSets))
    
  g <- delete.edges(g, E(g)[as.numeric(list.edges[, 1]) < edge.cutoff])
  index.deg <- igraph::degree(g) >= degree.cutoff
  g <- delete.vertices(g, V(g)[!index.deg])
  x <- x[index.deg,]
  index <- index[index.deg]
  if (length(V(g)) == 0)
  {
    stop("no categories greater than degree.cutoff...")
  }
  n <- min(nrow(x),n)
  x <- x[1:n,]
  
  group.level <- sort(unique(group))
  pvalues <- x[,pvalue]
  for(i in 1:length(group.level))
  {
    index <- x[,'Group'] == group.level[i]
    V(g)$shape[index] <- group.shape[i]
      
    group.pvalues <- pvalues[index]
    if(length(group.pvalues) > 0)
    {
      if(max(group.pvalues) == min(group.pvalues))
      {
        V(g)$color[index] <- adjustcolor(group.color[i], alpha.f = 0.5)
      }
      else
      {
        V(g)$color[index] <- sapply(1-(group.pvalues-min(group.pvalues))/
                                      (max(group.pvalues)-min(group.pvalues)),
                                      function(x) {adjustcolor(group.color[i],
                                                               alpha.f = x)})
      }
    }
  }
    
  plot(g, ...)
  if(show.legend)
  {
    legend.parameter$legend <- group.level
    legend.parameter$text.col <- group.color
    do.call(legend, legend.parameter)
  }
  return(g)
}
