\name{net.kernel}
\alias{net.kernel}

\title{ Kernel matrix for a graph }
\description{This function gives the regularized Laplacian matrix for a graph.}
\usage{
net.kernel(g, decay = 0.5)
}

\arguments{

\item{g}{
graph
}
\item{decay}{
decaying constant (default: 0.5)
}
}

\value{
the regularized Laplacian matrix
}

\author{Dongmin Jung, Xijin Ge}


\seealso{laplacian_matrix}
\examples{
# example 1
\dontrun{
string.db.9606 <- STRINGdb$new(version = '10', species = 9606,
                               score_threshold = 999)
string.db.9606.graph <- string.db.9606$get_graph()
K.9606 <- net.kernel(string.db.9606.graph)
}

# example 2
data(litG)
litG <- igraph.from.graphNEL(litG)
sg <- decompose(litG, min.vertices=50)
sg <- sg[[1]]
K <- net.kernel(sg)
}


