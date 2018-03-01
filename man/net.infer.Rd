\name{net.infer}
\alias{net.infer}

\title{
Inferring functionally related proteins using networks
}
\description{
Proteins can be classified by using networks to identify functionally closely related proteins.
}
\usage{
net.infer(target, kernel, top = NULL, cross = 0,
          C = 1, nu = 0.2, epsilon = 0.1, cache1 = 40,
          tol1 = 0.001, shrinking1 = TRUE, cache2 = 40,
          tol2 = 0.001, shrinking2 = TRUE)
}

\arguments{
  \item{target}{
set of interesting proteins or target class
}
  \item{kernel}{
the regularized Laplacian matrix for a graph
}
  \item{top}{
number of top proteins most closely related to target class  (default: all proteins except for target and pseudo-absence class)
}
  \item{cross}{
if a integer value k>0 is specified, a k-fold cross validation on the training data is performed to assess the quality of the model
}
  \item{C}{
cost of constraints violation for SVM (default: 1)
}
  \item{nu}{
The nu parameter for OCSVM (default: 0.2)
}
  \item{epsilon}{
epsilon in the insensitive-loss function for OCSVM (default: 0.1)
}
  \item{cache1}{
cache memory in MB for OCSVM (default: 40)
}
  \item{tol1}{
tolerance of termination criterion for OCSVM (default: 0.001)
}
  \item{shrinking1}{
option whether to use the shrinking-heuristics for OCSVM (default: TRUE)
}
  \item{cache2}{
cache memory in MB for SVM (default: 40)
}
  \item{tol2}{
tolerance of termination criterion for SVM (default: 0.001)
}
  \item{shrinking2}{
option whether to use the shrinking-heuristics for SVM (default: TRUE)
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
cross validation error, (when cross > 0)
}

\item{top}{
top proteins
}

\item{score}{
decision values for top proteins
}
}

\references{
Senay, S. D. et al. (2013). Novel three-step pseudo-absence selection technique for improved species distribution modelling. PLOS ONE. 8(8), e71218.
}
\author{
Dongmin Jung, Xijin Ge
}



\seealso{
ksvm
}
\examples{
# example 1
\dontrun{
string.db.9606 <- STRINGdb$new(version = '10',species = 9606,
                               score_threshold = 999)
string.db.9606.graph <- string.db.9606$get_graph()
K.9606 <- net.kernel(string.db.9606.graph)
rownames(K.9606) <- substring(rownames(K.9606), 6)
colnames(K.9606) <- substring(colnames(K.9606), 6)
target <- colnames(K.9606)[1:100]
infer <- net.infer(target, K.9606, 10)
}

# example 2
data(litG)
litG <- igraph.from.graphNEL(litG)
sg <- decompose(litG, min.vertices = 50)
sg <- sg[[1]]
K <- net.kernel(sg)
litG.infer <- net.infer(names(V(sg))[1:10], K, top=20)
}


