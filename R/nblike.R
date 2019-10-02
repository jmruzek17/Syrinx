#' Title
#'
#' @param theta Yo no se?
#' @param n IDK, Will wrote this one
#'
#' @return Something that I thnk has to do with degree distros
#' @export
#'
#'
nblike <- function(theta, n)
{
  occi <- colSums(n > 0)
  ni=colSums(n[,which(colSums(n)>0)])
  k<- theta[1]
  m <-nrow(n)
  propi <- 1-(1+ni/(m*k))^-k
  logl <- sum(occi*log(propi)+(m-occi)*log(1-propi))
  return(-logl)
}
