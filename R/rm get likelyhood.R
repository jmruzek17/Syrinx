
#' Title
#'
#' @param rand.mat A Square correlation matrix, AKA adj. mat
#' @param nr.thresholds 101
#' @param unfold.method Gaussian
#' @param bandwidth "nrd0"
#' @param nr.fit.points 101
#' @param dist.method "LL"
#' @param nr.breaks 101
#' @param discard.outliers TRUE
#' @param discard.zeros FALSE
#' @param min.mat.dim Minimum number of speices for calculations to proceed. Lower for datasets with low diversity
#' @param max.ev.spacing 3
#' @param interval NULL
#' @param interactive FALSE
#' @param smooth.par 0.5
#' @param plot.comp FALSE
#' @param save.fit FALSE
#' @param plot.spacing FALSE
#' @param wait.seconds Not important
#'
#' @return The output from \code{\link{cbind}}
#' @return The output from \code{\link{print}}
#' @export
#'

rm.get.likelyhood <- function (rand.mat, nr.thresholds = 101, unfold.method = "gaussian",
                               bandwidth = "nrd0", nr.fit.points = 101, dist.method = "LL",
                               nr.breaks = 101, discard.outliers = TRUE, discard.zeros = FALSE,
                               min.mat.dim = 40, max.ev.spacing = 3, interval = NULL, interactive = F,
                               smooth.par = 0.5, plot.comp = F, save.fit = FALSE, plot.spacing = FALSE,
                               wait.seconds = 0)
{
  if (!is.matrix(rand.mat))
    stop("\n\n rm.get.threshold: 'rand.mat' must be a matrix.\n\n")
  if (nrow(rand.mat) != ncol(rand.mat))
    stop("\n\n rm.get.threshold: 'rand.mat' must be a quadratic matrix.\n\n")
  if (!isSymmetric(rand.mat))
    stop("\n\n rm.get.threshold: 'rand.mat' must be a symmetric matrix.\n\n")
  if (max.ev.spacing <= sqrt(2/pi))
    stop("\n\n rm.get.threshold: 'max.ev.spacing' should not be lower than the abscissa of the peak of the Wigner-Dyson distribution (which is approximately 0.8).\n\n")
  if (max.ev.spacing > 5)
    cat("\n   WARNING: parameter 'max.ev.spacing' is quite big. Might provoke errors.\n\n")
  if (nrow(rand.mat) < 100)
    cat("\n   WARNING: matrix 'rand.mat' seems to be quite small for this approach but let's give it a try ...\n\n")
  if (!unfold.method %in% c("gaussian", "spline"))
    stop("\n\n  rm.get.threshold: parameter 'unfold.method' must be 'gaussian' or 'spline'.\n\n")
  if (!dist.method %in% c("LL", "KLD"))
    stop("\n\n  rm.get.threshold: parameter 'dist' must be 'LL' (log likelihood) or 'KLD' (Kullback-Leibler).\n\n")
  min.cell = min(abs(rand.mat[upper.tri(rand.mat, diag = F)]))
  max.cell = max(abs(rand.mat[upper.tri(rand.mat, diag = F)]))
  if (!is.null(interval)) {
    if (!is.numeric(interval) | (length(interval) != 2))
      stop("\n\n rm.get.threshold: 'interval' must be a two-component numeric vector.\n\n")
    request = ((min(interval) >= min.cell) || (min(interval) ==
                                                 0)) && (max(interval) <= max.cell)
    if (!request) {
      cat(paste("\n  Range of the absolute values of the matrix elements:",
                signif(min.cell, 5), " to ", signif(max.cell,
                                                    5), "\n\n"))
      stop("\n  rm.get.threshold: parameter 'interval' must be inside the range of the absolute values of the matrix elements.\n\n")
    }
    thresholds = seq(min(interval), max(interval), len = nr.thresholds)
  }
  else {
    thresholds = seq(min.cell, max.cell, len = nr.thresholds)
  }
  N = nrow(rand.mat)
  cat(paste("\n ", N, "times", N, "symmetric matrix read.\n\n"))
  results = list()
  results[["unfold.method"]] = unfold.method
  results[["dist.method"]] = dist.method
  results[["tested.thresholds"]] = numeric(0)
  results[["dist.Wigner"]] = numeric(0)
  results[["dist.Expon"]] = numeric(0)
  results[["nr.zeros"]] = integer(0)
  results[["nr.uniq.ev"]] = integer(0)
  results[["max.ev.mult"]] = integer(0)
  results[["nr.spacings"]] = integer(0)
  results[["nr.small.spacings"]] = integer(0)
  results[["perc.small.spacings"]] = integer(0)
  results[["eff.dimension"]] = integer(0)
  if (plot.comp)
    results[["comparison.plots"]] = character(0)
  if (save.fit)
    results[["cumfit.plots"]] = character(0)
  if (plot.spacing)
    results[["space.plots"]] = character(0)
  if (discard.zeros)
    results[["rm.dimension"]] = integer(0)
  if (discard.outliers)
    results[["nr.outliers.removed"]] = integer(0)
  results[["p.ks"]] = integer(0)
  results[["sse.exp"]] = integer(0)
  for (i in 1:nr.thresholds) {
    loop.start.time <- Sys.time()
    thres = thresholds[i]
    cat(paste(" ---> Loop =", i, "  threshold =", signif(thres,
                                                         3), "\n"))
    diagon = diag(rand.mat)
    rand.mat[which(abs(rand.mat) < abs(thres), arr.ind = T)] = 0
    diag(rand.mat) = diagon
    eff.mat = rm.discard.zeros(rand.mat, silent = T)
    cat(paste("      Effective matrix size =", nrow(eff.mat),
              "\n"))
    if (nrow(eff.mat) < min.mat.dim) {
      cat(paste("\n  Remaining number of non-zero rows & columns is below",
                min.mat.dim, "...\n  Breaking loop.\n\n"))
      break
    }
    if (discard.zeros)
      rand.mat = eff.mat
    if (save.fit)
      fn.fit = paste("RMT.Fit", i, "png", sep = ".")
    else fn.fit <- NULL
    res <- rm.ev.unfold(rand.mat, unfold.method = unfold.method,
                        bandwidth = bandwidth, nr.fit.points = nr.fit.points,
                        discard.outliers = discard.outliers, fn = fn.fit,
                        pop.up = FALSE, silent = TRUE)
    ev.spacing = res$ev.spacing
    ev = res$eigenvalues
    l1 = length(ev.spacing)
    if (!is.null(max.ev.spacing))
      ev.spacing = ev.spacing[ev.spacing <= max.ev.spacing]
    l2 = length(ev.spacing)
    cat(paste("      Number of large spacings not considered ( larger than",
              max.ev.spacing, ") :", l1 - l2, "\n"))
    epsilon = max.ev.spacing/1000
    nr.small.spacings = sum(ev.spacing < epsilon)
    perc.small.spacings = nr.small.spacings/length(ev.spacing) *
      100
    cat(paste("      Percentage of small spacings ( <", epsilon,
              ") =", round(perc.small.spacings, 2), "\n"))
    p.val.ks.test = ks.test(unique(ev.spacing), "pexp", 1)$p.value
    sse.exp = rm.sse(ev.spacing)
    results[["eff.dimension"]][i] = nrow(eff.mat)
    results[["nr.small.spacings"]][i] = nr.small.spacings
    results[["perc.small.spacings"]][i] = perc.small.spacings
    results[["nr.spacings"]][i] = length(ev.spacing)
    results[["cumfit.plots"]][i] = fn.fit
    results[["tested.thresholds"]][i] = thres
    results[["nr.zeros"]][i] = sum(rand.mat == 0)
    results[["nr.uniq.ev"]][i] = length(unique(ev))
    results[["max.ev.mult"]][i] = max(table(ev))
    if (discard.zeros)
      results[["rm.dimension"]][i] = nrow(rand.mat)
    if (discard.outliers)
      results[["nr.outliers.removed"]][i] = res$nr.outliers.removed
    results[["p.ks"]][i] = p.val.ks.test
    results[["sse.exp"]][i] = sse.exp
    if (plot.spacing) {
      fn = paste("RMT.Spaceplot", i, "png", sep = ".")
      rm.spacing.scatter(ev.spacing, pop.up = F, fn = fn)
      results[["space.plots"]][i] = fn
    }
    dres <- rm.get.distance(ev.spacing, dist.method = dist.method,
                            nr.breaks = nr.breaks)
    results[["dist.Wigner"]][i] = dres$dist.Wigner
    results[["dist.Expon"]][i] = dres$dist.Expon
    if (plot.comp) {
      fn = paste("RMT.Spacing", i, "png", sep = ".")
      rm.spacing.distribution(ev.spacing, threshold = thres,
                              dist.Wigner = dres$dist.Wigner, dist.Expon = dres$dist.Expon,
                              fn = fn)
      results[["comparison.plots"]][i] = fn
    }
    loop.end.time <- Sys.time()
    loop.time = as.numeric(loop.end.time - loop.start.time)
    if (loop.time < wait.seconds)
      Sys.sleep(wait.seconds - loop.time)
  }
  thresholds = results[["tested.thresholds"]]
  nr.unique.ev = results[["nr.uniq.ev"]]
  max.ev.mult = results[["max.ev.mult"]]
  nr.zeros = results[["nr.zeros"]]
  nr.spacings = results[["nr.spacings"]]
  dist.Wigner = results[["dist.Wigner"]]
  dist.Expon = results[["dist.Expon"]]
  sum.sq.err = results[["sse.exp"]]
  if (discard.zeros)
    mat.dim = results[["rm.dimension"]]
  fn = paste("RMT.numzeros", "png", sep = ".")
  png(fn)
  mtxt = "Number of zero matrix-elements vs. threshold"
  plot(thresholds, nr.zeros, col = "blue", main = mtxt, font.main = 1,
       xlab = "threshold", ylab = "nr zeros")
  dev.off()
  results[["number.zeros.plot"]] = fn
  fn = paste("RMT.num.uniq.ev", "png", sep = ".")
  png(fn)
  mtxt = "Number of unique eigenvalues vs. threshold"
  plot(thresholds, nr.unique.ev, col = "blue", main = mtxt,
       font.main = 1, xlab = "threshold", ylab = "nr unique ev")
  dev.off()
  results[["number.uniq.ev.plot"]] = fn
  fn = paste("RMT.max.ev.mult", "png", sep = ".")
  png(fn)
  mtxt = "Maximum eigenvalue multiplicity vs. threshold"
  plot(thresholds, max.ev.mult, col = "blue", main = mtxt,
       font.main = 1, xlab = "threshold", ylab = "max. ev multiplicity")
  dev.off()
  results[["max.ev.mult.plot"]] = fn
  if (discard.zeros) {
    fn = paste("RMT.mat.dimension", "png", sep = ".")
    png(fn)
    mtxt = "Dimension of non-zero matrix vs. threshold"
    plot(thresholds, mat.dim, col = "blue", main = mtxt,
         font.main = 1, xlab = "threshold", ylab = "matrix dimension")
    dev.off()
    results[["mat.dimension.plot"]] = fn
  }
  fn = paste("RMT.num.ev.spacings", "png", sep = ".")
  png(fn)
  mtxt = "Number of ev spacings vs. threshold"
  plot(thresholds, nr.spacings, col = "blue", main = mtxt,
       font.main = 1, xlab = "threshold", ylab = "nr. ev spacings")
  dev.off()
  results[["num.ev.spacings.plot"]] = fn
  fn = paste("RMT.Dist.vs.Thres", "png", sep = ".")
  if (dist.method == "LL")
    main.res <- rm.likelihood.plot(thresholds, log.le = dist.Expon,
                                   log.lw = dist.Wigner, smooth.par = smooth.par, fn = fn,
                                   interactive = F)
  if (dist.method == "KLD")
    main.res <- rm.distance.plot(thresholds, dist.Expon = dist.Expon,
                                 dist.Wigner = dist.Wigner, smooth.par = smooth.par,
                                 fn = fn, interactive = F)

  likelyhood <- cbind(thresholds, dist.Expon, dist.Wigner )

  return(likelyhood)
}

