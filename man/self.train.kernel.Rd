\name{self.train.kernel}
\alias{self.train.kernel}

\title{ Self training for a kernel matrix }
\description{This function can be used for classification of semi-supervised data by using the kernel support vector machine.}
\usage{
self.train.kernel(K, y, type = 'response', C = 1, cache = 40,
                  tol = 0.001, shrinking = TRUE, thrConf = 0.9,
                  maxIts = 10, percFull = 1, verbose = FALSE) 
}

\arguments{
  \item{K}{
kernel matrix
}
  \item{y}{
lable vector
}
  \item{type}{
one of response, probabilities ,votes, decision indicating the type of output (default: response)
}
  \item{C}{
cost of constraints violation for SVM (default: 1)
}
  \item{cache}{
cache memory in MB for SVM (default: 40)
}
  \item{tol}{
tolerance of termination criterion for SVM (default: 0.001)
}
  \item{shrinking}{
option whether to use the shrinking-heuristics for OCSVM (default: TRUE)
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
A boolean indicating the verbosity level of the function (default: FALSE)
}
}

\value{
prediction from the SVM
}

\author{Dongmin Jung, Xijin Ge}


\references{
Torgo, L. (2016) Data Mining using R: learning with case studies, second edition, Chapman & Hall/CRC.
}
\examples{
data(litG)
litG <- igraph.from.graphNEL(litG)
sg <- decompose(litG, min.vertices = 50)
sg <- sg[[1]]
K <- net.kernel(sg)
y <- rep(NA, length(V(sg)))
y[1:10] <- 1
y[11:20] <- 0
y <- factor(y)
self.train.kernel(K, y)
}


