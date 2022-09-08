#' EQF plot
#'
#' Depict the EQF plot by the result of permutation process to detect
#' the QTL hotspot.
#'
#' @param result list. The data list of the output from LOD.QTLdetect()
#' or EQF.permu().
#' @param plot.all logical. If being TURE, output one figure of the
#' EQF values over the bins.
#' @param plot.cr logical. If being TURE, output the figures of the
#' EQF values over the bins of each chromosome.
#'
#' @return
#'
#' One or several EQF plots.
#'
#' @export
#'
#' @references
#'
#' Wu, P.-Y., M.-.H. Yang, and C.-H. KAO 2021 A Statistical Framework
#' for QTL Hotspot Detection. G3: Genes, Genomes, Genetics: jkab056.
#'
#' @seealso
#' \code{\link[QTLEMM]{LOD.QTLdetect}}
#' \code{\link[QTLEMM]{EQF.permu}}
#'
#' @examples
#' # load the example data
#' load(system.file("extdata", "LODexample.RDATA", package = "QTLEMM"))
#'
#' # run and result
#' EQF.plot(LOD.QTLdetect.result)
#' EQF.plot(EQF.permu.result)
EQF.plot <- function(result, plot.all = TRUE, plot.cr = TRUE){

  dat <- result
  name0 <- names(dat)
  if(length(name0) == 6){
    datatest <- name0 != c("detect.QTL.number", "QTL.matrix", "EQF.matrix",
                           "linkage.QTL.number", "LOD.threshole", "bin")
  } else if (length(name0) == 9){
    datatest <- name0 != c("EQF.matrix", "bin", "LOD.threshole", "cluster.number", "cluster.id",
                           "cluster.matrix", "permu.matrix.cluster", "permu.matrix.Q", "EQF.threshold")
  } else {
    stop("Input data error, please input the original output data of LOD.QTLdetect or EQF.permu.", call. = FALSE)
  }

  if(TRUE %in% (datatest)){
    stop("Input data error, please input the original output data of LOD.QTLdetect or EQF.permu.", call. = FALSE)
  }

  if(!plot.all[1] %in% c(0,1) | length(plot.all) > 1){plot.all <- TRUE}
  if(!plot.cr[1] %in% c(0,1) | length(plot.cr) > 1){plot.cr <- TRUE}

  EQF <- dat$EQF.matrix
  clumatrix <- dat$cluster.matrix
  thre <- dat$LOD.threshole
  eqfthre <- dat$EQF.threshold
  bin <- dat$bin
  nc <- nrow(bin)
  lcr <- bin[, 2]
  ncr <- c()
  for(i in 1:nc){
    ncr[i] <- sum(bin[1:i, 2])
  }
  cr0 <- c()
  for(i in 1:nc){
    cr0 <- c(cr0, rep(i, bin[i, 2]))
  }
  eqf.all <- apply(EQF, 2, sum)

  x0 <- 1:length(cr0)+(cr0-1)*100
  xn <- ncr+(0:(nc-1))*100-lcr/2

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  if(plot.all & nc>1){
    graphics::par(mfrow = c(1, 1))
    graphics::par(mai = c(1, 1, 1, 1))
    yli <- max(eqf.all)*1.2
    ma <- paste("LOD thresholds =", thre)
    if(length(clumatrix) > 0){
      ma <- paste("LOD thresholds =", thre, "  # of group =", nrow(clumatrix))
      graphics::par(mai = c(1, 1, 1, 1.5))
    }
    plot(x0, eqf.all, type = "h", ylab = "EQF", xlab = "chromosome", main = ma,
         xaxt = "n", ylim = c(0, yli), yaxt = "n",  cex.main = 2, cex.lab = 1.2, axes = FALSE)
    graphics::axis(side = 1, pos = -yli/35, at = xn, labels = 1:nc, cex.axis = 1.2, tick = FALSE)
    lse <- 4000/max(x0)
    graphics::segments(x0, rep(-yli/35, length(x0)), x0, rep(-yli/70, length(x0)), lwd = lse)
    graphics::segments(-100, 0, -100, yli)
    graphics::segments(-100, yli, max(x0)+100, yli)
    graphics::segments(max(x0)+100, 0, max(x0)+100, yli)
    if(length(eqfthre) > 0){
      for(j in 1:nrow(eqfthre)){
        graphics::axis(2, eqfthre[j, 1], las = 2)
        graphics::axis(4, paste(rownames(eqfthre)[j], " (", eqfthre[j, 2], ")", sep = ""), at = eqfthre[j], las = 2)
      }
      graphics::abline(h = eqfthre[, 1], col = "red")
    } else {graphics::axis(2, seq(0, yli, 20), las = 2)}
  }

  if(plot.cr | (nc == 1 & plot.all)){
    if(nc >= 4){
      graphics::par(mai = c(0.8, 0.8, 0.8, 0.8))
      if(length(clumatrix) > 0){
        graphics::par(mai = c(0.8, 0.8, 0.8, 1.2))
      }
    } else {
      graphics::par(mai = c(1, 1, 1, 1.5))
      if(length(clumatrix) > 0){
        graphics::par(mai = c(1, 1, 1, 1))
      }
    }
    for(i in 1:nc){
      eqf <- eqf.all[cr0 == i]
      plot(eqf, type = "h", ylab = "EQF", xlab = "position(bin)", main = paste("cr", i, "   LOD thresholds = ", thre),
           ylim = c(0, max(eqf.all)*1.2), yaxt = "n", cex.main = 2, cex.lab = 1.2)
      if(length(eqfthre)>0){
        for(k in 1:nrow(eqfthre)){
          graphics::axis(2, eqfthre[k, 1], las = 2.5)
          graphics::axis(4, paste(rownames(eqfthre)[k], " (", eqfthre[k, 2], ")", sep = ""), at = eqfthre[k], las = 2.5)
        }
        graphics::abline(h = eqfthre[, 1], col = "red")
      } else {graphics::axis(2, seq(0, yli, 20), las = 2.5)}
    }
  }
}


