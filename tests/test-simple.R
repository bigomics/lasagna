library(data.table)
library(ggplot2)
library(igraph)
library(devtools)
#load_all("~/Playground/playbase")
#load_all("~/Projects/WGCNAplus")
load_all()

## Imports
gset.rankcor <- playbase::gset.rankcor
mat2gmt <- playbase::mat2gmt
ai.ask <- playbase::ai.ask
ai.create_image_gemini <- playbase::ai.create_image_gemini
lasagna.multisolve <- playbase::lasagna.multisolve

data <- playbase::mofa.exampledata("brca")
#data <- playbase::mofa.exampledata("cll")

ntop=1000;nc=20;add.sink=TRUE;intra=TRUE;use.gmt=FALSE;use.graphite=0;
fully_connect=FALSE;add.revpheno=TRUE;condition.edges=1

names(data)
names(data$X)
##names(data$X) <- substring(names(data$X),1,2)

obj <- create_model(data, pheno="pheno", ntop=1000, nc=10,
  add.sink=TRUE, intra=FALSE, fully_connect=FALSE, add.revpheno=TRUE,
  condition.edges=1)
names(obj)

## color by WGCNA clustering
if(0) {

  names(data$X)
  wgcna <- computeWGCNA_multiomics(data$X, data$samples,
    power=12, minmodsize = 3, minKME=0.1, mergeCutHeight = 0.5)
  names(wgcna)
  names(wgcna$layers)
  
  wgcna$me.genes
  table(wgcna$me.colors)
  head(wgcna$me.colors)
  
  V(obj$graph)$color <- wgcna$me.colors[V(obj$graph)$name]
  ##V(obj$graph)$color <- substring(wgcna$me.colors[V(obj$graph)$name],3,99)
  V(obj$graph)$color[is.na(V(obj$graph)$color)] <- "red"
  table(V(obj$graph)$color)

}

wgcna$me.genes

## solve the graph for a certain phenotype
colnames(obj$Y)
pheno = "activated=act"
pheno = "condition=Her2"
graph <- lasagna::solve(obj, pheno, min_rho=0.01, max_edges=1000,
  value="rho", sp.weight=1, prune=FALSE) 
graph

## prune graph for cleaner plotting
pdf("multipartite-moxbrca.pdf", w=14, h=8)
par(mfrow=c(1,1), mar=c(1,1,1,1)*0)
mp <- lasagna::plot_multipartite(
  graph,
  min.rho = 0.1,
  ntop = 40,
  xdist = 1,
  color.var = "color",
  labpos = c(2,2,4,4),
  cex.label = 0.8,
  vx.cex = 1.1,
  edge.cex = 1.5,
  edge.alpha = 0.4,
  edge.sign = "both",
  edge.type = "inter",
  edge.gamma = 2,
  yheight = 0.99,
  normalize.edges = 1,
  strip.prefix = TRUE,
  prune = 0
) 
dev.off()

names(mp)
mp$graph
table(V(mp$graph)$color)


## interactive multipartite
M <- mp$layout
vis <- plot_visgraph(mp$graph, layers=NULL, ntop=-1,
  min_rho=0.2, ecex=3, vcex=3, labcex=1, egamma=2,
  color.var="color", mst=0, layout=M, physics=0) 
vis


## MST - FR layout (NEED RETHINK)
require(visNetwork)
mst <- mp$graph
##ew <- 1 / pmax(E(mp$graph)$weight,0)
ew <- 1 / (1e-8 + abs(E(mp$graph)$weight)**2)
mst <- igraph::mst(mp$graph, weights=ew)
vis <- plot_visgraph(mst, layers=NULL, ntop=-1,
  min_rho=0.1, ecex=3, egamma=2, vcex=2, labcex=1, 
  mst=1, color.var="color", layout=NULL, physics=TRUE) 
vis

## hierarchical layout
vis %>% visHierarchicalLayout()

##-------------------------------------------------
## 3D lasagna
##-------------------------------------------------
##source("~/Playground/playbase/dev/include.R", chdir=TRUE)
load_all("..")

X=obj$X
xpos <- layout_multipartite_3d(graph, obj$X, clust='svd')
xpos <- layout_multipartite_3d(graph, obj$X, clust='tsne')
xpos <- layout_multipartite_3d(graph, obj$X, clust='umap')

plot_3d(graph, layout=xpos, draw_edges=TRUE,
  color.by = "color", min_rho=0.0, sign_rho="pos",
  cex=1.5, cex.gamma=0.5, num_edges=100, znames=NULL) 

  
##-------------------------------------------------
##-------------------------------------------------
##-------------------------------------------------
