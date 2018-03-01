\name{GSEA.barplot}
\alias{GSEA.barplot}

\title{
Visualize the gene set enrichment analysis
}
\description{
For the functional enrichment analysis, we can visualize the result from the gene set enrichment analysis.
}
\usage{
GSEA.barplot(object, category, score, pvalue, top = 10,
             sort = NULL, decreasing = FALSE, numChar = NULL, 
             title = NULL, transparency = 0.5, plot = TRUE)
}

\arguments{
  \item{object}{
a table with category, enrichment score and p-value of gene sets
}
  \item{category}{
name of gene sets
}
  \item{score}{
enrichment score
}
  \item{pvalue}{
p-value of gene sets
}
  \item{top}{
the number of top categories (default: 10)
}
  \item{sort}{
a variable used for sorting data
}
  \item{decreasing}{
logical indicating whether ascending or descending order (default: FALSE)
}
  \item{numChar}{
the maximal number of characters of the name of gene sets
}
\item{title}{
title for the plot
}
  \item{transparency}{
transparency (default: 0.5)
}
  \item{plot}{
return plot when plot is true, otherwise return table (default: TRUE)
}
}


\value{
GSEA barplot
}



\references{
Yu G, Wang L, Yan G and He Q (2015). "DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis." Bioinformatics, 31(4), pp. 608-609.
}

\author{
Dongmin Jung, Xijin Ge
}



\seealso{
ggplot2
}

\examples{
data(examplePathways)
data(exampleRanks)
set.seed(1)
result.GSEA <- fgsea(examplePathways, exampleRanks, nperm = 1000)
GSEA.barplot(result.GSEA, category = 'pathway', score = 'NES',
             pvalue = 'pval', sort = 'NES', decreasing = TRUE)
}



