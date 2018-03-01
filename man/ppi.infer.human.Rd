\name{ppi.infer.human}
\alias{ppi.infer.human}

\title{
Inferring functionally related proteins using protein networks for human
}
\description{
This function is designed for human protein-protein interaction from STRING database. Default format is 'hgnc'. The number of proteins is 10 in default. Note that the number of proteins used as a target may be different from the number of proteins in the input since mapping between formats is not always one-to-one in getBM.
}
\usage{
ppi.infer.human(target, kernel, top = 10, classifier = net.infer,
                input = "hgnc_symbol", output = "hgnc_symbol", ...)
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
# example 1
string.db.9606 <- STRINGdb$new(version = '10', species = 9606,
                               score_threshold = 999)
string.db.9606.graph <- string.db.9606$get_graph()
K.9606 <- net.kernel(string.db.9606.graph)
rownames(K.9606) <- substring(rownames(K.9606), 6)
colnames(K.9606) <- substring(colnames(K.9606), 6)
target <- colnames(K.9606)[1:100]
infer.human <- ppi.infer.human(target, K.9606, input = "ensembl_peptide_id")

\dontrun{
# example 2
library(graph)
data(apopGraph)
target <- nodes(apopGraph)
apoptosis.infer <- ppi.infer.human(target, K.9606, 100)

# example 3
library(KEGGgraph)
library(KEGG.db)
pName <- "p53 signaling pathway"
pId <- mget(pName, KEGGPATHNAME2ID)[[1]]
getKGMLurl(pId, organism = "hsa")
p53 <- system.file("extdata/hsa04115.xml", package="KEGGgraph")
p53graph <- parseKGML2Graph(p53,expandGenes=TRUE)

entrez <- translateKEGGID2GeneID(nodes(p53graph))
ensembl <- useMart("ensembl")
human.ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
target <- getBM(attributes=c('entrezgene', 'hgnc_symbol'),
                filter = 'entrezgene', values = entrez,
                mart = human.ensembl)[,2]
p53.infer <- ppi.infer.human(target, K.9606, 100)
}
}

