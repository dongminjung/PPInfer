\name{ORA.barplot}
\alias{ORA.barplot}

\title{
Visualize the over-representation analysis
}
\description{
For the functional enrichment analysis, we can visualize the result from the over-representation analysis.
}
\usage{
ORA.barplot(object, category, size, count, pvalue, top = 10,
            sort = NULL, decreasing = FALSE, p.adjust.methods = NULL,
            numChar = NULL, title = NULL, transparency = 0.5,
            plot = TRUE)
}

\arguments{
  \item{object}{
a table with category, size, count and p-value of gene sets
}
  \item{category}{
name of gene sets
}
  \item{size}{
size of gene sets
}
  \item{count}{
count of gene sets
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
  \item{p.adjust.methods}{
a correction method
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
ORA barplot
}



\references{
Yu G, Wang L, Yan G and He Q (2015). "DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis." Bioinformatics, 31(4), pp. 608-609.
}

\author{
Dongmin Jung, Xijin Ge
}



\seealso{
p.adjust, ggplot2
}

\examples{
data(examplePathways)
data(exampleRanks)
geneNames <- names(exampleRanks)
set.seed(1)
gene.id <- sample(geneNames, 100)
result.ORA <- ORA(examplePathways, gene.id)
ORA.barplot(result.ORA, category = "Category", size = "Size",
            count = "Count", pvalue = "pvalue", sort = "pvalue")
}



