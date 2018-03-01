\name{ORA}
\alias{ORA}

\title{
Over-representation Analysis
}
\description{
the result from the over-representation analysis
}
\usage{
ORA(pathways, gene.id, minSize = 1, maxSize = Inf,
    p.adjust.methods = NULL)
}

\arguments{
  \item{pathways}{
list of gene sets
}
  \item{gene.id}{
set of genes
}
  \item{minSize}{
Minimal size of a gene set
}
  \item{maxSize}{
Maximal size of a gene set
}
  \item{p.adjust.methods}{
a correction method
}
}


\value{
ORA result
}


\author{
Dongmin Jung, Xijin Ge
}



\seealso{
fisher.test
}

\examples{
data(examplePathways)
data(exampleRanks)
geneNames <- names(exampleRanks)
set.seed(1)
gene.id <- sample(geneNames, 100)
ORA(examplePathways, gene.id)
}



