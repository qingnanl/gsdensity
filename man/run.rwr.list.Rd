% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsdensity_functions.R
\name{run.rwr.list}
\alias{run.rwr.list}
\title{result in a matrix (number of rows = number of cells; number of columns = number of gene sets)
reflecting the probability each cell is labeled during the
propagation (relevance to the gene set); same idea as run.rwr but with multiple
gene sets}
\usage{
\method{run}{rwr.list}(el, gene_set_list, cells)
}
\arguments{
\item{el:}{edge list; output of 'compute.nn.edges'}

\item{gene_set:}{a vector of genes of interest}

\item{cells:}{name of cells; usually the same as 'colnames(object)'}
}
\description{
result in a matrix (number of rows = number of cells; number of columns = number of gene sets)
reflecting the probability each cell is labeled during the
propagation (relevance to the gene set); same idea as run.rwr but with multiple
gene sets
}
\examples{
rs <- sample(length(gene.set.list), 10, replace = F)
cv.df <- run.rwr.list(el = el, gene_set_list = gene.set.list[rs], cells = cells)
}
