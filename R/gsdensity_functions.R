
#' 1. compute MCA embeddings
#'

#' @param object a seurat object
#' @param dims.use which mca dimensions to use; default is the first 10 dimensions
#' @param genes.use which genes to use; default is all genes in the object
#' @return returns a dataframe with cells as rows and mca coordinates as columns
#' @export
#' @import Seurat dplyr Matrix CelliD
#'
#'
compute.mca <- function(object, dims.use = 1:10, genes.use = rownames(object)){

  genes.use.1 <- intersect(genes.use, rownames(object))
  object <- object %>%
    NormalizeData() %>%
    RunMCA(features = genes.use.1)
  genes.use.update <- intersect(genes.use.1,
                                rownames(object@reductions$mca@feature.loadings))
  coembed <- rbind(object@reductions$mca@feature.loadings[genes.use.update, dims.use],
                   object@reductions$mca@cell.embeddings[, dims.use])
  return(coembed)
}

#' 2. compute density of gene sets of interest
#' 2.1 compute grid point coordinates

#' @param coembed the result from compute.mca
#' @param genes.use which genes to use; no default;
#' can use genes based on the gene set selection or use rownames(object)
#' @param n.grids number of grid points used for gene set density estimation;
#' larger number is more accurate and slower;
#' default is 100 (recommended to test 100 first)
#' @return grid coordinates
#' @export
#' @import anticlust
#'
compute.grid.coords <- function(coembed, genes.use, n.grids = 100){
  coembed <- scale(coembed[genes.use, ])
  cl <- balanced_clustering(coembed, K = n.grids)
  centroid.coords <- aggregate(coembed, list(cl), mean)[, -1]
  return(centroid.coords)
}

#' 2.2 compute KL-divergence
#' (some are adapted from https://github.com/alexisvdb/singleCellHaystack/)

#' @param coembed the result from compute.mca
#' @param genes.use which genes to use; no default;
#' can use genes based on the gene set selection or use rownames(object)
#' @param n.grids number of grid points used for gene set density estimation;
#' larger number is more accurate and slower;
#' default is 100 (recommended to test 100 first)
#' 'coembed', 'genes.use', 'n.grids' are passed to 'compute.grid.coords()'
#' @param gene.set.list a list of gene sets;
#' e.g., gene.set.list <- list(gene.set.a = c("A", "B", "C"),
#'                             gene.set.b = c("a", "b", "c"))
#' @param gene.set.cutoff gene sets with length less than this cutoff will
#' not be used; the length is after the intersection of the gene set and
#' genes.use
#' @param n.times to evaluate how likely the gene set density is not caused
#' by randomness, size-matched gene sets will be used to compute the background
#' density distribution; This simulation will be done n.times; default is 100
#' @return kl-divergence between given gene set and random gene sets
#' @export
#' @import
#' @examples
#' library(SeuratData)
#' library(Seurat)
#' data('pbmc3k')
#' res <- compute.kld(coembed = ce,
#'                    genes.use = intersect(rownames(ce), rownames(pbmc3k)),
#'                    n.grids = 100,
#'                    gene.set.list = gene.set.list[1:5],
#'                    gene.set.cutoff = 3,
#'                    n.times = 100)
#'
#'
compute.kld <- function(coembed, genes.use,
                        n.grids = 100,
                        gene.set.list,
                        gene.set.cutoff = 3,
                        n.times = 100){
  coembed.scale <- scale(coembed[genes.use, ])
  grid.co <- compute.grid.coords(coembed = coembed,
                                 genes.use = genes.use,
                                 n.grids = n.grids)
  # for background
  dist.to.grid <- vectorized_pdist(A = as.matrix(coembed.scale), B = as.matrix(grid.co))
  bandwidth <- median(apply(dist.to.grid,1,min))
  dist.to.grid.norm <- dist.to.grid / bandwidth
  density.contributions <-
    exp(-dist.to.grid.norm * dist.to.grid.norm / 2)
  Q <- compute.db(density.df = density.contributions)

  # compute kld for each gene set and randomly sampled gene sets of the same size
  gene.set.names <- rep(NA, length(gene.set.list))
  klds <- rep(0, length(gene.set.list))
  len.gl <- rep(0, length(gene.set.list))
  rklds.avg <- rep(0, length(gene.set.list))
  rklds.sd <- rep(0, length(gene.set.list))
  pvalues <- rep(0, length(gene.set.list))
  for (i in 1:length(gene.set.list)){
    gene.set.name <- names(gene.set.list)[i]
    gene.set <- intersect(genes.use, gene.set.list[[i]])
    #                print(gene.set)
    len.gene.set <- length(gene.set)
    if (len.gene.set < gene.set.cutoff){
      next
    }
    gene.set.names[i] <- gene.set.name
    len.gl[i] <- len.gene.set
    P <- compute.db(density.df = density.contributions[gene.set, ])
    kld <- sum(P * log(P/Q))
    klds[i] <- log(kld)
    # sample
    rkld <- sapply(1:n.times, function(x)sample.kld(density.df = density.contributions,
                                                    ref = Q,
                                                    len.gene.set = len.gene.set))
    rkld.avg <- mean(log(rkld))
    rkld.sd <- sd(log(rkld))
    rklds.avg[i] <- rkld.avg
    rklds.sd[i] <- rkld.sd
    pvalue <- pnorm(log(kld), rkld.avg, rkld.sd, lower.tail = FALSE)
    pvalues[i] <- pvalue
  }
  out <- as.data.frame(cbind(gene.set.names,
                             klds,
                             len.gl,
                             rklds.avg,
                             rklds.sd,
                             pvalues))
  out <- out[complete.cases(out), ]
  colnames(out) <- c("gene.set",
                     "kld",
                     "gene.set.length",
                     "rkld.mean",
                     "rkld.sd",
                     "p.value")
  out$p.adj <- p.adjust(p = out$p.value, method = "fdr")

  return(out)
}

#' from an excellent post: https://www.r-bloggers.com/2013/05/pairwise-distances-in-r/
#' enhanced the speed
#' this function is called by 'compute.kld' to quickly compute the distance between
#' genes to grid points
#' 
#' @param A matrix
#' @param B matrix
#' @return returns pairwise-distances
#' @export
vectorized_pdist <- function(A,B){
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))

  m = nrow(A)
  n = nrow(B)

  tmp = matrix(rep(an, n), nrow=m)
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  sqrt( tmp - 2 * tcrossprod(A,B) )
}

#' this function is called by 'compute.kld' to aggregate the density contribution
#' of each gene to each grid point, and then normalize the densities of grid points
#' to 1.
#' @param density.df an intermediate object in 'compute.kld'
#' @return distribution
#' @export

compute.db <- function(density.df){
  Q <- apply(density.df, 2, sum)
  pseudo <- 1e-300
  Q <- Q + pseudo
  Q <- Q / sum(Q)
  return(Q)
}

#' this function is called by 'compute.kld' to calculate the kl-divergence between
#' sampled (background) gene set and the ref (all) gene set
#'
#' @param density.df density.df
#' @param ref ref
#' @param len.gene.set len.gene.set
#' @return returns random klds
#' @export
sample.kld <- function(density.df, ref, len.gene.set){
  idx <- sample(nrow(density.df), len.gene.set, replace = F)
  rP <- compute.db(density.df = density.df[idx, ])
  rkld.vec <- sum(rP * log(rP/ref))
  return(rkld.vec)
}

#' 3. compute nearest neighbor graph for genes and cells
#' This graph will be used for fetching the most relevant cells of a gene set
#'

#' @param coembed the result from compute.mca
#' @param nn.use the number of nearest neighbors for building the graph; default 300
#' @return nearest neighbor graph (edges)
#' @export
#' @import RANN igraph

#'
compute.nn.edges <- function(coembed, nn.use = 300){
  nbrs <- nn2(coembed, k = nn.use)
  el_idx <- el_nn_search(nn2_out = nbrs)
  el_nn <- cbind(rownames(coembed)[el_idx[, 1]], rownames(coembed)[el_idx[, 2]])
  return(el_nn)
}

#' this function is called by 'compute.nn.edges' to convert nearest neighbor
#' identity matrix to edge list
#' @param nn2_out nn2_out
#' @return returns edge list
#' @export
el_nn_search <- function(nn2_out){
  n.df <- nn2_out$nn.idx
  n.df <- cbind(1:nrow(n.df), n.df)
  el <- cbind(n.df[, 1], c(n.df[, -1]))
  return(el)
}

#' 4. compute label propagation from gene set to cells

#' this function is to form a 'seed matrix' used by the dRWR function (dnet R package);
#' the seed matrix is specifying which nodes are the sources for label propagation
#' @param gene_set gene_set
#' @param graph.use graph.use
#' @return returns seed matrix
#' @export
seed.mat <- function(gene_set, graph.use){
  gs <- intersect(gene_set, names(V(graph.use)))
  ss <- data.frame(n = names(V(graph.use)))
  rownames(ss) <- ss$n
  ss$weight <- ifelse(rownames(ss) %in% gs, 1, 0)
  ss$n <- NULL
  return(ss)
}

#' this function is used when more than one 'seed sets' will be used (when there
#' are multiple gene sets of interest)
#' @param gene_set_list gene_set_list
#' @param graph.use graph.use
#' @return returns seed matrix
#' @export
seed.mat.list <- function(gene_set_list, graph.use){
  sml <- future.apply::future_lapply(names(gene_set_list),
                                     function(x) seed.mat(gene_set = gene_set_list[[x]],
                                                          graph.use = graph.use))
  sm <- do.call(cbind, sml)
  colnames(sm) <- names(gene_set_list)
  sm <- as.matrix(sm)
  return(sm)
}

#' 4.1 To calculate the label propagation probability for a gene set among cells;
#' result in a vector (length = number of cells) reflecting the probability
#' each cell is labeled during the propagation (relevance to the gene set)
#'
#'

#' @param el edge list; output of 'compute.nn.edges'
#' @param gene_set a vector of genes of interest
#' @param cells name of cells; usually the same as 'colnames(object)'
#' @param restart the probability of the propagation to restart
#' @return cell vector (representing gene set activity)
#' @export
#' @import future future.apply dnet


run.rwr <- function(el, gene_set, cells, restart = 0.75){
  g <- graph_from_edgelist(as.matrix(el), directed = F)
  g <- simplify(g)
  ss <- seed.mat(gene_set = gene_set, graph.use = g)
  rwr <- dRWR(g, setSeeds = ss, parallel = F, verbose = F, restart = restart)
  rownames(rwr) <- rownames(ss)
  cell_vec <- rwr[cells, 1]
  cell_vec <- cell_vec / sum(cell_vec)
  return(cell_vec)
}

#' result in a matrix (number of rows = number of cells; number of columns = number of gene sets)
#' reflecting the probability each cell is labeled during the
#' propagation (relevance to the gene set); same idea as run.rwr but with multiple
#' gene sets
#'

#' @param el edge list; output of 'compute.nn.edges'
#' @param gene_set_list a list of gene sets
#' @param cells name of cells; usually the same as 'colnames(object)'
#' @param restart the probability of the propagation to restart
#' @return activity of pathways in cells
#' @export


run.rwr.list <- function(el, gene_set_list, cells, restart = 0.75){
  g <- graph_from_edgelist(as.matrix(el), directed = F)
  g <- simplify(g)
  ss <- seed.mat.list(gene_set_list = gene_set_list, graph.use = g)
  rwrl <- dRWR(g, setSeeds = ss, parallel = T, verbose = F, restart = restart)
  rownames(rwrl) <- rownames(ss)
  cell_rwr <- rwrl[cells, ]
  colnames(cell_rwr) <- colnames(ss)
  cell_rwr <- cell_rwr[, colSums(is.na(cell_rwr)) == 0]
  cell_rwr_norm <- cell_rwr %*% diag(1/colSums(cell_rwr))

  rownames(cell_rwr_norm) <- rownames(cell_rwr)
  colnames(cell_rwr_norm) <- colnames(cell_rwr)
  cell_rwr_norm <- as.data.frame(cell_rwr_norm)
  return(cell_rwr_norm)
}

#' 4.2. binarize the label propagation probability in the cell population;
#' result in a binarized vector of cells with 'nagative' and 'positive' labels;
#' 'positive' means that the cells are relevant to the gene set
#'
#'
#' @param cell_vec output of 'run.rwr'
#' @return cell label of 'negative' or 'positive' for a given pathway
#' @export
#' @import multimode

#'
compute.cell.label <- function(cell_vec){
  nm <- names(cell_vec)
  m <- locmodes(cell_vec, mod0 = 2)
  split.value <- m$locations[2]
  cell.label <- ifelse(cell_vec < split.value, "negative", "positive")
  return(cell.label)
}

#' similar to compute.cell.label; used when working with multiple gene sets
#' @param cell_df output of 'run.rwr.list'
#' @return cell labels of 'negative' or 'positive' for given pathways
#' @export
#' @import multimode

compute.cell.label.df <- function(cell_df){
  cell.labels <- future_apply(cell_df,
                              MARGIN = 2,
                              function(x) {compute.cell.label(x)
                              })
  return(cell.labels)
}


#' 5. compute the specificity of gene set when cell partition information is available;
#' the information could be clustering, sample origins, or other conditions
#' inspired by https://github.com/FloWuenne/scFunctions/blob/0d9ea609fa72210a151f7270e61bdee008e8fc88/R/calculate_rrs.R
#'
#' this function is called by compute.spec.single to calculate the similarity
#' between two vectors
#' @param x x
#' @param y y
#' @return returns jsd_distance
#' @export
#'
compute.jsd <- function(x, y){
  input_df <- rbind(x, y)
  jsd_divergence <- suppressMessages(philentropy::JSD(input_df))
  jsd_distance <- 1-sqrt(jsd_divergence)
  return(jsd_distance)
}

#' This is to calculate the similarity between:
#' 1. the label propagation probability of cells for gene sets and
#' 2. the identify of cells in a certain partition
#' This is called by 'compute.spec'; can also run by itself
#' @param vec cell partition vector (usually a column name in object@meta.data)
#' @param positive the positive label, e.g. "disease" or "cluster_1"
#' @param cell_df the output of run.rwr.list
#' @return specificity of a pathway activity and other levels of cell annotations (e.g., cell type)
#' @export
#' @import infotheo philentropy

#'

compute.spec.single <- function(vec, positive, cell_df){
  num <- ifelse(vec == positive, 1, 0)
  num <- num / sum(num)
  spec.single <- future_apply(cell_df,
                              MARGIN = 2,
                              function(x) {compute.jsd(x = x, y = num)})
  names(spec.single) <- colnames(cell_df)
  return(spec.single)
}

#' This is to calculate the similarity between:
#' 1. the label propagation probability of cells for gene sets and
#' 2. the identify of cells in partitions
#' @param cell_df the output of run.rwr.list
#' @param metadata a data frame with cell information (each row is a cell;
#' usually object@meta.data)
#' @param cell_group cell partition vector (usually a column name
#' @return specificity of a pathway activity and other levels of cell annotations (e.g., cell type)
#' in object@meta.data)
#' @export


compute.spec <- function(cell_df, metadata, cell_group){
  gene.sets <- colnames(cell_df)
  cells <- rownames(cell_df)
  metadata <- metadata[cells, ]
  cell_groups <- unique(as.character(metadata[[cell_group]]))

  jsd <- future.apply::future_lapply(cell_groups,
                                     function(x)
                                     {compute.spec.single(vec = metadata[[cell_group]],
                                                          positive = x,
                                                          cell_df = cell_df)})
  jsd.df <- do.call(cbind, jsd)
  colnames(jsd.df) <- cell_groups
  return(jsd.df)
}

#' 6. find gene sets with spatial relevance
#'
#'
#' This function is to calculate how likely the cells relevant to a gene set is
#' randomly distributed spatially
#'

#' @param spatial.coords a data frame with each row as a cell and each column
#' as a spatial coordinate (usually 2: x and y)
#' @param weight_vec output of run.rwr
#' @param n split the spatial map for local density estimation;
#' n is the number of split for each dimension; for n = 10, the spatial map is
#' split to n * n = 100 grids for the density estimation
#' @param n.times the weight_vec is shuffled several times (n.times) to generate
#' a background distribution (shuffled weights vs. equal weights) for statistical
#' significance estimation (p.value); larger n.times will be more time-consuming and
#' more accurate
#' @return spatial kl-divergence 
#' @export
#' @import MASS
#'
compute.spatial.kld <- function(spatial.coords, weight_vec, n = 10, n.times = 20){
  bg.weight <- rep(1/nrow(spatial.coords), nrow(spatial.coords))
  bg.dens <- kde2d.weighted(x = spatial.coords[, 1],
                            y = spatial.coords[, 2],
                            w = bg.weight,
                            n = n)
  Q <- c(bg.dens$z) + 1e-300
  dens <- kde2d.weighted(x = spatial.coords[, 1],
                         y = spatial.coords[, 2],
                         w = weight_vec,
                         n = n)
  P <- c(dens$z)
  spatial.kld <- sum(P * log(P/Q))
  rspatial.kld <- sapply(1:n.times, function(x)sample.spatial.kld(weight_vec = weight_vec,
                                                                  ref = Q,
                                                                  n = n,
                                                                  spatial.coords = spatial.coords))
  rspatial.kld.avg <- mean(log(rspatial.kld))
  rspatial.kld.sd <- sd(log(rspatial.kld))
  pvalue <- pnorm(log(spatial.kld), rspatial.kld.avg, rspatial.kld.sd, lower.tail = FALSE)
  spatial.kld.vec <- c(log(spatial.kld), rspatial.kld.avg, rspatial.kld.sd, pvalue)
  names(spatial.kld.vec) <- c("spatial.kld", "rspatial.kld", "rspatial.kld.sd", "p.value")
  return(spatial.kld.vec)
}

#' This function is to calculate how likely the cells relevant to
#' multiple gene sets are randomly distributed spatially
#'

#' @param spatial.coords a data frame with each row as a cell and each column
#' as a spatial coordinate (usually 2: x and y)
#' @param weight_df output of run.rwr.list
#' @param n split the spatial map for local density estimation;
#' n is the number of split for each dimension; for n = 10, the spatial map is
#' split to n * n = 100 grids for the density estimation
#' @param n.times the same as n.times in function 'compute.spatial.kld'
#' @return spatial kl-divergence for multiple gene sets
#' @export
#'
compute.spatial.kld.df <- function(spatial.coords, weight_df, n = 10, n.times = 20){
  weight_df <- weight_df[rownames(spatial.coords), ] # make sure the order is the same
  klds <- future_apply(weight_df,
                       MARGIN = 2, future.seed=TRUE,
                       function(x) {compute.spatial.kld(spatial.coords = spatial.coords,
                                                        weight_vec = x,
                                                        n = n, n.times = n.times)})
  klds.df <- as.data.frame(t(klds))
  #klds.df <- do.call(rbind, klds)
  #rownames(klds.df) <- colnames(weight_df)
  klds.df$p.adj <- p.adjust(p = klds.df$p.value, method = "fdr")
  return(klds.df)
}

#' this function is called by 'compute.spatial.kld' to calculate the kl-divergence between
#' cell-weighted with shuffled weight vector and the ref (all cells, unweighted)
#'
#' @param weight_vec weight_vec
#' @param spatial.coords spatial.coords
#' @param n n
#' @param ref ref
#' @return returns randomly sampled spatial klds for gene sets
#' @export
sample.spatial.kld <- function(weight_vec, spatial.coords, n, ref){
  rweight_vec <- sample(weight_vec)
  rdens <- kde2d.weighted(x = spatial.coords[, 1],
                          y = spatial.coords[, 2],
                          w = rweight_vec,
                          n = n)
  rP <- c(rdens$z)
  rspatial.kld <- sum(rP * log(rP/ref))
  return(rspatial.kld)
}

#' based on https://stat.ethz.ch/pipermail/r-help/2006-June/107405.html
#' this is called by compute.spatial.kld to calculate the kernel density estimation
#' in 2d space with each data point weighted.
#'
#' @param x x
#' @param y y
#' @param w w
#' @param h h
#' @param n n
#' @param lims lims
#' @return weighted kde2d estimation
#' @export
kde2d.weighted <- function (x, y, w, h, n, lims = c(range(x), range(y))) {
  nx <- length(x)
  if (length(y) != nx)
    stop("data vectors must be the same length")
  gx <- seq(lims[1], lims[2], length = n) # gridpoints x
  gy <- seq(lims[3], lims[4], length = n) # gridpoints y
  if (missing(h))
    h <- c(bandwidth.nrd(x), bandwidth.nrd(y));
  if (missing(w))
    w <- numeric(nx)+1;
  h <- h/4
  ax <- outer(gx, x, "-")/h[1] # distance of each point to each grid point in x-direction
  ay <- outer(gy, y, "-")/h[2] # distance of each point to each grid point in y-direction
  z <- (matrix(rep(w,n), nrow=n, ncol=nx, byrow=TRUE)*matrix(dnorm(ax), n, nx)) %*% t(matrix(dnorm(ay), n, nx))/(sum(w) * h[1] * h[2]) # z is the density
  return(list(x = gx, y = gy, z = z))
}

