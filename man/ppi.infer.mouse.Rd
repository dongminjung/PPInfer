\name{ppi.infer.mouse}
\alias{ppi.infer.mouse}

\title{
Inferring functionally related proteins using protein networks for mouse
}
\description{
This function is designed for mouse protein-protein interaction from STRING database. Default format is 'mgi'. The number of proteins is 10 in default. Note that the number of proteins used as a target may be different from the number of proteins in the input since mapping between formats is not always one-to-one in getBM.
}
\usage{
ppi.infer.mouse(target, kernel, top = 10, classifier = net.infer,
                input = "mgi_symbol", output = "mgi_symbol", ...)
}

\arguments{
  \item{target}{
set of interesting proteins or target class
}
  \item{kernel}{
the regularized Laplacian matrix for a graph
}
  \item{top}{
number of top proteins most closely related to target class (default: 10)
}
  \item{classifier}{
net.infer or net.infer.ST (default: net.infer)
}
  \item{input}{
  input format
}
  \item{output}{
  output format
}
  \item{...}{
additional parameters for the chosen classifier
}
}


\value{
\item{list}{
list of a target class used in the model
}

\item{error}{
training error
}

\item{CVerror}{
cross validation error, (when cross > 0 in net.infer)
}

\item{top}{
top proteins
}

\item{score}{
decision values for top proteins
}
}

\author{
Dongmin Jung, Xijin Ge
}




\seealso{
net.infer, net.infer.ST, getBM
}
\examples{
string.db.10090 <- STRINGdb$new(version = '10', species = 10090,
                                score_threshold = 999)
string.db.10090.graph <- string.db.10090$get_graph()
K.10090 <- net.kernel(string.db.10090.graph)
rownames(K.10090) <- substring(rownames(K.10090), 7)
colnames(K.10090) <- substring(colnames(K.10090), 7)
target <- colnames(K.10090)[1:100]
infer.mouse <- ppi.infer.mouse(target, K.10090, input="ensembl_peptide_id")
}

