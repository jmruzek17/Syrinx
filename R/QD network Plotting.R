
#' Title
#'
#' @param Tall.Mat A thresholded ADJ mat that has been melted into tall format (3 columns, origin, destination and weight of the edge)
#' @param MAINTEXT A title for the plot
#'
#' @return The output from \code{\link{plot}} An ugly, but representitive network
#' @export
#'
#'
QD_Network_Plot <- function(Tall.Mat, MAINTEXT = ""){

  ##Start with a tall correlation matrix, with 3 columns, source, target and weight

  mDXX.cor <- Tall.Mat

  ## remove the diagonal
  mDXX.cor$value[mDXX.cor$value==1] <- 0
  DXX      <- mDXX.cor[which(mDXX.cor$value!=0),]
  DXX.Th   <- DXX

  ##Rename columns, this is important
  DXX.Tha <- DXX.Th[]
  colnames(DXX.Tha) <- c("Source", "Target", "Data")

  DlinksXX <- aggregate(DXX.Tha[,3], DXX.Tha[,-3], sum)             #Edges are links
  DlinksXX <- DlinksXX[order(DlinksXX$Source, DlinksXX$Target),] #Should be set up as Source, Weight, Attributes
  colnames(DlinksXX)[3] <- c("weight")
  rownames(DlinksXX)     <- NULL

  NodesXX   <- as.matrix(DXX.Tha[,2])        # Nodes are species
  NodesXX   <- cbind(NodesXX, "Dump","Trump") # I could put node attribuites here, like guild, but for now...
  NodesXX   <- unique(NodesXX)                # get rid of duplicate nodes

  #Create a network, and remove symetric connections
  DnetXX <- graph_from_data_frame(d=DlinksXX, vertices=NodesXX, directed=F)
  DnetXX <- igraph::simplify(DnetXX, remove.multiple = T, remove.loops = T)

  # Plot, there is a bunch of layout options, fr is the best
  l <- layout_with_dh(DnetXX)
  #l <- layout_in_circle(DnetXX)
  V(DnetXX)$size <- 2.5
  V(DnetXX)$frame.color <- "black"
  V(DnetXX)$color <- "black"
  V(DnetXX)$label <- ""
  E(DnetXX)$arrow.mode <- 0
  E(DnetXX)$color <- "grey"
  Kanye   <- plot(DnetXX, rescale = T, layout = l*0.0003, main = MAINTEXT)

  cfg.XX  <- cluster_fast_greedy((DnetXX), weights = abs((E(DnetXX)$weight)))
  Yeezus  <- plot(cfg.XX, y = DnetXX, rescale = T, layout = l, main = MAINTEXT)

  return(Kanye)

}





#' Title
#'
#' @param Network  A Network from iGraph package
#' @param MAINTEXT A title for the plot
#'
#' @return The output from \code{\link{plot}} An ugly, but representitive network
#' @export
#'

QD_Network_Plot_with_Network <- function(Network, MAINTEXT = ""){

  DnetXX  <- Network

  DnetXX <- igraph::simplify(DnetXX, remove.multiple = T, remove.loops = T)
  # Plot, there is a bunch of layout options, fr is the best
  l <- layout_with_dh(DnetXX)
  #l <- layout_in_circle(DnetXX)
  V(DnetXX)$size <- 2.5
  V(DnetXX)$frame.color <- "black"
  V(DnetXX)$color <- "black"
  V(DnetXX)$label <- ""
  E(DnetXX)$arrow.mode <- 0
  E(DnetXX)$color <- "grey"
  Kanye   <- plot(DnetXX, rescale = T, layout = l*0.0003, main = MAINTEXT)

  cfg.XX  <- cluster_fast_greedy((DnetXX), weights = abs((E(DnetXX)$weight)))
  Yeezus  <- plot(cfg.XX, y = DnetXX, rescale = T, layout = l, main = MAINTEXT)

  return(Kanye)

}
