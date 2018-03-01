\name{net.infer.ST}
\alias{net.infer.ST}

\title{
Inferring functionally related proteins with self training
}
\description{This function is the self-training version of net.infer. The function net.infer is the special case of net.infer.ST where a single iteration is conducted.}
\usage{
net.infer.ST(target, kernel, top = NULL, C = 1, nu = 0.2,
            epsilon = 0.1, cache1 = 40, tol1 = 0.001, shrinking1 = TRUE,
            cache2 = 40, tol2 = 0.001, shrinking2 = TRUE, thrConf = 0.9,
            maxIts = 10, percFull = 1, verbose = FALSE) 
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
  \item{thrConf}{
A number between 0 and 1, indicating the required classification confidence for an unlabelled case to be added to the labelled data set with the label predicted predicted by the classification algorithm (default: 0.9)
}
  \item{maxIts}{
The maximum number of iterations of the self-training process (default: 10)
}
  \item{percFull}{
A number between 0 and 1. If the percentage of labelled cases reaches this value the self-training process is stoped (default: 1)
}
  \item{verbose}{
A boolean indicating the verbosity level of the function. (default: FALSE)
}
}


\value{
\item{list}{
list of a target class used in the model
}

\item{error}{
training error
}

\item{top}{
top proteins
}

\item{score}{
decision values for top proteins
}
}


\author{Dongmin Jung, Xijin Ge}


\seealso{self.train}
\examples{
data(litG)
litG <- igraph.from.graphNEL(litG)
sg <- decompose(litG, min.vertices = 50)
sg <- sg[[1]]
K <- net.kernel(sg)
litG.infer.ST <- net.infer.ST(names(V(sg))[1:10], K, top=20)
}


