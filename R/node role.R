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
node.role <- function(Thresh.Adj.Mat, threshold, nnodes.mod = 5)
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

  #Now we calculate the within module z-scores and participation coefficients
  #at the network level first, because the P coefficient needs information
  #about the whole network, including modules with < 5 species.  This code is divided
  #into two sections.  First, we need to define within and out of module connections
  #and then count, for each species, how many inbound/out of module connections they have
  #Second, we then need to calculate for each species their Z and pcoef scores and then
  #assign the network roles for that species.

  #Build the species matrix and calculate in/out of module degree/strength

  #Edge type designation accounting for module connection types#
  edge.df      <- data.frame(as_edgelist(weight.g,names=TRUE))
  edge.df[,1]  <- as.character(edge.df[,1])
  edge.df[,2]  <- as.character(edge.df[,2])
  edge.df2     <- data.frame(edge.df[,1],seq(1,nrow(edge.df),1))
  colnames(edge.df2) <- c("V1","V2")
  edge.df3     <- data.frame(cbind(V(weight.g)$name,V(weight.g)$module))
  colnames(edge.df3) <- c("V1","V2")
  edge.df4     <- merge(edge.df2,edge.df3,by="V1")
  edge.df4     <- edge.df4[order(edge.df4$V2.x),]
  edge.df5     <- data.frame(edge.df[,2],seq(1,nrow(edge.df),1))
  colnames(edge.df5) <- c("V1","V2")
  edge.df6     <- data.frame(cbind(V(weight.g)$name,V(weight.g)$module))
  colnames(edge.df6) <- c("V1","V2")
  edge.df7     <- merge(edge.df5,edge.df6,by="V1")
  edge.df7     <- edge.df7[order(edge.df7$V2.x),]
  edge.dflist  <- data.frame(cbind(edge.df4,edge.df7))
  colnames(edge.dflist)  <- c("spp1","index","module1","spp2","index","module2")
  edge.dflist$modcontype <- ifelse(edge.dflist$module1==edge.dflist$module2,"InModule","OutModule")
  E(weight.g)$connectiontype <- edge.dflist$modcontype

  #Now build the species matrix
  spp.mat    <- igraph::degree(weight.g)
  spp.strmat <- igraph::strength(weight.g)
  mod.conmat <- as.data.frame(matrix(nrow=length(spp.mat),ncol=7))

  for(p in 1:nrow(mod.conmat)){
    mod.conmat[p,1] <- names(spp.mat)[p]
    mod.conmat[p,2] <- as.numeric(igraph::degree(weight.g)[p])

    Snoop    <- paste0(mod.conmat[p,1])
    iso.sps1 <- edge.dflist[(edge.dflist$spp1 == Snoop),]
    iso.sps2 <- edge.dflist[(edge.dflist$spp2 == Snoop),]

    iso.sps  <- rbind(iso.sps1, iso.sps2)

    mod.conmat[p,3] <- nrow(iso.sps[which(iso.sps$modcontype == "OutModule"),])
    mod.conmat[p,4] <- nrow(iso.sps[which(iso.sps$modcontype == "InModule"),])

    #
    #
    spp.edge        <- table(E(weight.g)$connectiontype)

  }

# for(p in 1:nrow(mod.conmat)){
#    mod.conmat[p,1] <- names(spp.mat)[p]
#    mod.conmat[p,2] <- as.numeric(igraph::degree(weight.g)[p])
#    spp.edge        <-table(E(weight.g)$connectiontype)#
#
#    if(length(spp.edge)==1){
#      if(names(spp.edge)=="InModule"){
#        mod.conmat[p,3] <- spp.edge
#      }else {
#        mod.conmat[p,4] <- spp.edge}}
#    else{
#      for(pp in 1:length(spp.edge)){
#        if(names(spp.edge[pp])=="InModule"){
#          mod.conmat[p,3] <- spp.edge[pp]}else {mod.conmat[p,4] <- spp.edge[pp]
#        }
#      }
#   }
#  }

  for(p in 1:nrow(mod.conmat)){
    mod.conmat[p,5] <- as.numeric(igraph::strength(weight.g))[p]
    spp.edge        <- table(E(weight.g)$connectiontype)

    if(length(spp.edge)==1){
      if(names(spp.edge)=="InModule"){
        mod.conmat[p,6]      <- sum(E(weight.g)[from(p)]$weight)
      }else {mod.conmat[p,7] <- sum(E(weight.g)[from(p)]$weight)}}
    else{
      for(pp in 1:length(spp.edge)){
        if(names(spp.edge[pp])=="InModule"){
          df <- data.frame(E(weight.g)[from(p)]$connectiontype,E(weight.g)[from(p)]$weight)
          mod.conmat[p,6] <- sum(df[which(df[,1]=="InModule"),2])
          }
        else {
          df <-data.frame(E(weight.g)[from(p)]$connectiontype,E(weight.g)[from(p)]$weight)
          mod.conmat[p,7] <- sum(df[which(df[,1]=="OutModule"),2])
        }
      }
    }
  }

  mod.conmat[is.na(mod.conmat)] <- 0
  netintramod.deg        <- mod.conmat
  netintramod.deg$module <- oliclus.membs

  #calculation of node within module z-scores to define hub vs. nonhub
  netintramod.deg$z.score <- within_module_deg_z_score(unweight.g, as.numeric(oliclus.membs))

  #calculation of node degree participation coefficient
  netintramod.deg$pcoef   <- part_coeff(unweight.g, as.numeric(oliclus.membs))


  #separate the species based on valid module size
  species.goodmods <- species[which(species$V2 %in% species.goodmods$V2),]
  species.badmods2 <- species[which(species$V2 %in% species.badmods),]
  goodintramod.deg <- netintramod.deg[netintramod.deg[,1] %in% species.goodmods$V1,]
  badintramod.deg  <- netintramod.deg[netintramod.deg[,1] %in% species.badmods2$V1,]

  #labels for species
  goodintramod.deg$role <-ifelse(goodintramod.deg$z.score > 2.5 & goodintramod.deg$pcoef > 0.62,"NetHub",
                             ifelse(goodintramod.deg$z.score > 2.5 & goodintramod.deg$pcoef < 0.62, "ModHub",
                                ifelse(goodintramod.deg$z.score < 2.5 & goodintramod.deg$pcoef < 0.62, "Peripheral",
                                   ifelse(goodintramod.deg$z.score < 2.5 & goodintramod.deg$pcoef > 0.62,"Connector","Nominal")
                                      )
                                    )
                                  )

  #Now we repeat, but for species who fall into modules of < 5 species, we need to
  #correctly label them at the network level, because it is possible for
  #these small species modules to be connected with other modules, thus
  #increasing their participation coefficient and making them possibly
  #connectors.

  badintramod.deg$role <- ifelse(badintramod.deg$z.score < 2.5 & badintramod.deg$pcoef < 0.62,
                                 "Peripheral","Connector")

  #now we combine the data, label the columns, and reorder it so that the
  #species are in the same order as the original network matrix

  goodbadintramod.mat           <- rbind(goodintramod.deg,badintramod.deg)
  colnames(goodbadintramod.mat) <- c("spp", "kTotal","kIn","kOut","Strength",
                                     "StrengthIn","StrengthOut","module","Z.score","pcoef",
                                     "Role")

  netintramod.deg2 <- goodbadintramod.mat[order(match(rownames(goodbadintramod.mat),
                                                      rownames(netintramod.deg))), ]


  return(netintramod.deg2) #This should have the roles of the species now accurately defined
}
