#' Title
#'
#' @param Thresh.Adj.Mat A thresholded square adjecency matrix
#' @param threshold the threshold that thresholded the adjecency matrix
#' @param nnodes.mod number of nodes in the smallest possible module, default is 5
#'
#' @return the node roles based on the modified module size, and a the network with the outmodule vertices removed
#' @export
#'
#'
#'
Module.reducer <- function(Thresh.Adj.Mat, threshold, nnodes.mod = 5)
{

  #revised nole role labelling based on z-score and
  #participation coefficient

  require(igraph)
  require(brainGraph)
  require(WGCNA)
  #Read in a thresholded, no singleton square correlation matrix

  threshmat           <- Thresh.Adj.Mat
  rownames(threshmat) <- colnames(threshmat)
  oli.Threshold       <- threshold

  #Build a weighted graph and an unweighted graph
  weight.g        <- graph_from_adjacency_matrix(as.matrix(threshmat),mode="lower",weighted=TRUE,diag=FALSE)
  unweight.adjmat <- signumAdjacencyFunction(threshmat,oli.Threshold)
  unweight.g      <- graph_from_adjacency_matrix(as.matrix(unweight.adjmat),mode="lower",diag=FALSE)

  #Apply clustering algorithm on the weighted graph.  Should give same results as unweighted graph
  oliclus       <- cluster_fast_greedy(weight.g,weights = abs(E(weight.g)$weight))
  oliclus.membs <- as.factor(oliclus$membership) #make the group numbers factors
  oli.sim.mod   <- modularity(oliclus) #retain modularity
  V(weight.g)$module <- oliclus.membs #set the module ids for each vertex

  #Now we need to separate out the network so that labels are correctly
  #applied to species that are members of modules of nnodes.mod or more species.
  species          <- as.data.frame(matrix(nrow=nrow(unweight.adjmat),ncol=2))
  species[,1]      <- colnames(unweight.adjmat)
  species[,2]      <- oliclus.membs
  species.table    <- table(species$V2)
  species.goodmods <- as.numeric(names(species.table[which(species.table >= nnodes.mod)]))
  species.badmods  <- as.numeric(names(species.table[which(species.table < nnodes.mod)]))
  species.goodmods <- species[which(species$V2 %in% species.goodmods),]
  species.badmods2 <- species[which(species$V2 %in% species.badmods),]

  smaller.network  <- weight.g - vertices(species.badmods2$V1)

  return(smaller.network) #This should have the roles of the species now accurately defined
}
