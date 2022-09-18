# gsdensity: Density-based Gene Set Evaluation


### Introduction

Investigation of gene sets among cell populations is a common step in single-cell transcriptomic studies. Several hurdles do exist in this approach, such as: 

1. the noise/sparsity of the detected transcripts per cell; 
2. the uncertainty/bias in clustering and annotation; and 
3. the necessity of manual annotation (for factorization/pattern approaches). 

Here, we present 'gsdensity', a computational tool dedicated for 'pathway centric' analysis of single-cell data. Given a gene set and a cell-by-gene matrix, we want to ask: is this gene set somehow enriched by a subpopulation of the cells? If so, can we fetch the cells for downstream analysis?

The motivation is that in many cases, to analyze single cell data, there is not a standard way of finding out which cell population is the most 'interesting' and worth a deeper analysis. However, scientists of the very field always have insights regarding 'it will be interesting to look at this pathway in the data'. This is especially true when clustering and annotation is non-trivial for the dataset.

The core idea of gsdensity is that we use multiple correspondence analsysis (MCA) to co-embed cells and genes (thanks to the CelliD package) and then investigate the 'density' of the gene sets of interest in the MCA space. The position of each gene in the MCA space reflects the relationship between that gene and cells and other genes. Thus, when a gene set shows a very different 'density' in the MCA space compared with the background, it indicates that there does exist some 'affinity' between the gene set and some cell subpopulation. We can then fetch the most relevant cells for a gene set using the label propagation algorithm since the genes and cells can be placed in the same graph (nearest neighbor graph built with the MCA space), with the genes in the gene set being 'seeds'. Thus, gsdensity is not affected by any of the three hurdles brought up above:

1. The co-embedding of cells and genes will alleviate of the sparsity of single-cell data
2. This gsdensity method does not require any clustering information; however, when clustering/annotation is available, we can find highly specific gene sets for a cluster of interest (please check tutorials)
3. We can use well-annotated gene sets directly as the input; there is no need for manual annotation of patterns.

### Installation

Have been tested for windows and linux systems with no problem. Running on a linux server (with parallel computing) is strongly recommended when there are many gene sets to be tested. 

Installation not successful in macos with M1 chip, due to a dependency of the CelliD package (https://github.com/LTLA/scuttle/issues/14#issuecomment-989653173). Will keep testing and update when this is solved. (9/17: please check Update section #2)

```
# First, install dependencies: Seurat, CelliD, dnet, supraHex, Rgraphviz, infotheo, anticlust, multimode, philentropy; then:

# install.packages("remotes")

#Turn off warning-error-conversion, because the tiniest warning stops installation
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

#install from github
remotes::install_github("https://github.com/qingnanl/gsdensity")

```

### Usage

We included two examples showing how to use gsdensity with scRNA-seq data and spatial transcriptomics data. 

In the former one we introduced the basic workflow and main functions of gsdensity, including a quick comparison between gsdensity and gsea.

The latter one is mainly focused on how to associate gene sets with the spatial information of cells or other partition of cells.

[Using gsdensity to analyze scRNA-seq data](http://htmlpreview.github.io/?https://github.com/qingnanl/gsdensity/blob/master/vignette/pbmc3k_example.html)

[Using gsdensity to analyze spatial transcriptomics data](http://htmlpreview.github.io/?https://github.com/qingnanl/gsdensity/blob/master/vignette/spatial_example_10x_visium.html)

### Quick start

```
library(gsdensity)
library(ggplot2) # for plotting
library(reshape2)
library(msigdbr) # for gathering gene sets
library(Seurat)
library(SeuratData)
# library(future) # for parallel computing
# library(future.apply) # for parallel computing

# use GO_BP gene sets 
# Conver the format to 'list'

mdb_c5 <- msigdbr(species = "Homo sapiens", category = "C5")
mdb_c5_bp <- mdb_c5[mdb_c5$gs_subcat == "GO:BP", ]
gene.set.list <- list()
for (gene.set.name in unique(mdb_c5_bp$gs_name)){
        gene.set.list[[gene.set.name]] <- mdb_c5_bp[mdb_c5_bp$gs_name %in% gene.set.name, ]$gene_symbol
}

# from SeuratData
data("pbmc3k")
# Normalize the data
pbmc3k <- NormalizeData(pbmc3k)

# compute cell/gene embeddings
ce <- compute.mca(object = pbmc3k)
# find gene sets with differential density 
res <- compute.kld(coembed = ce, 
                   genes.use = rownames(pbmc3k), 
                   n.grids = 100, 
                   gene.set.list = gene.set.list,
                   gene.set.cutoff = 3,
                   n.times = 100)                   
gene.set.deviated <- res[res$p.adj < 0.05, ]$gene.set

# build nearest neighbor graph
cells <- colnames(pbmc3k)
el <- compute.nn.edges(coembed = ce, nn.use = 300)

# get label propagation probability for each cell of a gene set 'GOBP_B_CELL_ACTIVATION'
cv <- run.rwr(el = el, gene_set = gene.set.list[["GOBP_B_CELL_ACTIVATION"]], cells = cells)
# get a 'positive' or 'negative' label for each cell of this gene set
cl <- compute.cell.label(cv)

# Compute the UMAP coordinates for visualization
pbmc3k <- pbmc3k %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA() %>%
        FindNeighbors() %>%
        RunUMAP(dims = 1:20)
# plot the cell annotation
DimPlot(pbmc3k,
        group.by = "seurat_annotations",
        label = T,
        raster = T)
```
![image](https://user-images.githubusercontent.com/53788946/186203473-4b47d5c1-6543-4f48-a858-65bc3e2b2b49.png)


```
# add the label propagation probability and binarized label to meta data
pbmc3k@meta.data$b_cell_activition <- cv[colnames(pbmc3k)]
pbmc3k@meta.data$b_cell_activition_bin <- cl[colnames(pbmc3k)]

# plot the probabilities
FeaturePlot(pbmc3k,
        features = "b_cell_activition",
        raster = T)
```
![image](https://user-images.githubusercontent.com/53788946/186203736-fa9b03f8-716f-4275-908a-438e1b4b1799.png)

```
# plot the 'positive' or 'negative' labels        
DimPlot(pbmc3k,
        group.by = "b_cell_activition_bin",
        raster = T)        

```
![image](https://user-images.githubusercontent.com/53788946/186203804-ad29d828-3980-4b00-afe3-2cc6cddf7779.png)


### Updates

1. (9/15/22)When such error appears when running run.rwr():
```
Error: as(<dsCMatrix>, “dgCMatrix”) is deprecated since Matrix 1.5-0; do as(., “generalMatrix”) instead
```
please set:
```
igraph::igraph.options(sparsematrices = FALSE)
```
then run the code again

2. (9/17/22) For MacOS with M1 chip, the tool can be installed with the intel build R, but not for the arm64 version of R. This has been tested with an M1-2020 computer.
