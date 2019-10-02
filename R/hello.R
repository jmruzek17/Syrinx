

#' Title
#'
#' @param x the name of the person to say hi to
#'
#' @return The output from \code{\link{print}}
#' @export
#'
#' @examples
#' hello("Big Money")
#'
hello <- function(x) {
  print(paste("Hello", x, ", world!"))
}
