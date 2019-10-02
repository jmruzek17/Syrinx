


#' Title
#'
#' @param Thresholded.Matrix A thresholded matrix, species as columns and samples as rows.
#' @param abun.min  The minimal local abundance a specialist species must average in across all sites it occurs. Default is 0.02
#' @param upper.site.limit The maximal number of sites a specialist can occur in. Must be a proportion, defaults to 0.1
#' @param lower.site.limit The minimal number of sites a genralist can occur in. Must be a proportion, defaults to 0.5
#'
#' @return A matrix object with columns "No Sites", :Mean Abiundance" The mean abudance of that species across all sites it occurs,
#'         "Log abundance", the log of the total abundance of that species across all sites it occurs, and finally Specialist, generalist
#'         and nominal deliniation. A Generalist occurs in equal or more sites than the upper.site.limit, a specialist occurs in
#'         fewer or equal to sites than the lower.site.limit, and also must have a mean abundace greater than or equal to the
#'         abun.min.
#'         Nominal is all else.
#'
#' @export
#'
Spec.Gen <- function(Thresholded.Matrix, abun.min = 0.02, upper.site.limit = 0.5, lower.site.limit = 0.1){

  Thresh.mat  <- Thresholded.Matrix

  #Specialist/generalist classification, sensu Barberan et al. 2012

  sps.Pres <- as.matrix(colSums(decostand(Thresh.mat, method = "pa"))) #as.matrix(colSums(d2))
  nperm    <- length(Thresh.mat[1,])
  Maynard  <- matrix(ncol = 1, nrow = nperm)
  rownames(Maynard) <- colnames(Thresh.mat)

  for(i in 1:nperm){  ## Calculates the mean abundace of each species
    sps1 <- Thresh.mat[,i]
    sps2 <- sps1[sps1!=0]
    Maynard[i,] <- mean(sps2)
  }

  log_abun <- as.matrix((log(colSums(Thresh.mat))))
  Killer.Mike <- cbind(sps.Pres, Maynard, log_abun)
  colnames(Killer.Mike) <- c( "No.Sites", "Mean.abun", "Log.Abun.")
  Killer.Mike <- as.data.frame(Killer.Mike)
  #abun.min <- (0.02) ## species must be locally abundant, which barbaran classified as >2% of sequences (in our case density)
  upper.site.limit2 <- ((length(Thresh.mat[,1]))* upper.site.limit) # Barberan used 52% as his cut off, I rounded down to 50% of sites
  lower.site.limit2 <- ((length(Thresh.mat[,1]))* lower.site.limit)  ## Barberan used 10 sites,which was 6% of his sites, I round up to 10%
  Killer.Mike$Gen.Specialist <- ifelse(Killer.Mike$No.Sites >= upper.site.limit2 ,"Generalist",
                                       ifelse(Killer.Mike$Mean.abun >= abun.min & Killer.Mike$No.Sites <= lower.site.limit2 , "Specialist",
                                              "Nominal"))

  return(Killer.Mike)

}
