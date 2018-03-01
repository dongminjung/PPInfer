\name{enrich.net}
\alias{enrich.net}

\title{
Visualize network for the functional enrichment analysis
}
\description{
The connection between nodes depends on the proportion of overlapping genes between two categories. 
}
\usage{
enrich.net(x, gene.set, node.id, node.name = node.id, pvalue,
           n = 50, numChar = NULL, pvalue.cutoff = 0.05,
           edge.cutoff = 0.05, degree.cutoff = 0,
           edge.width = function(x) {10*x^2},
           node.size = function(x) {2.5*log10(x)},
           group = FALSE, group.color = c('red', 'green'),
           group.shape = c('circle', 'square'),
           legend.parameter = list('topright'),
           show.legend = TRUE, ...)
}

\arguments{
  \item{x}{
a result with category and p-value of gene sets
}
  \item{gene.set}{
gene sets which is already used for functional enrichment
}
  \item{node.id}{
name of gene sets
}
  \item{node.name}{
label of nodes in the network (default: node.id)
}
  \item{pvalue}{
pvalues for categories
}
  \item{n}{
number of top categories (default: 50)
}
  \item{numChar}{
the maximal number of characters of the label of gene sets
}
  \item{pvalue.cutoff}{
nodes with p-values which are greater than pvalue.cutoff are removed (default: 0.05)
}
  \item{edge.cutoff}{
edges with the proportion which is less than edge.cutoff are removed (default: 0.05)
}
  \item{degree.cutoff}{
nodes with the degrees which are less than degree.cutoff are removed (default: 0)
}
  \item{edge.width}{
width of edges
}
  \item{node.size}{
size of nodes
}
  \item{group}{
variable for group
}
  \item{group.color}{
color for group (default: red and green for 2 groups)
}
  \item{group.shape}{
shape for group (default: circle and square for 2 groups)
}
  \item{legend.parameter}{
list of parametres for the legend
}
  \item{show.legend}{
show the legend (default: TRUE)
}
  \item{...}{
additional parameters for the igraph
}  
}


\value{
plot for the network. The size of nodes is proportional to the size of gene sets. The more significant categories are, the less transparent their nodes are.
}



\references{
Yu G, Wang L, Yan G and He Q (2015). "DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis." Bioinformatics, 31(4), pp. 608-609.
}

\author{
Dongmin Jung, Xijin Ge
}



\seealso{
igraph
}

\examples{
data(examplePathways)
data(exampleRanks)
set.seed(1)
result.GSEA <- fgsea(examplePathways, exampleRanks, nperm = 1000)
enrich.net(result.GSEA, examplePathways, node.id = 'pathway',
           pvalue = 'pval', edge.cutoff = 0.6, degree.cutoff = 1,
           n = 50, vertex.label.cex = 0.75, show.legend = FALSE,
           edge.width = function(x) {5*sqrt(x)},
           layout = igraph::layout.kamada.kawai)
}



