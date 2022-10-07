#' QTL search by IM with Selective Genotyping
#'
#' Expectation-maximization algorithm for QTL interval mapping to search
#' the possible position of QTL in all chromosome. This function can
#' handle the genotype data witch is selective genotyping.
#'
#' @param marker matrix. A k*2 matrix contains the marker information,
#' where the row dimension k is the number of markers in the chromosomes.
#' The first column labels the chromosomes where the markers are located,
#' and the second column labels the positions of QTLs (in morgan (M) or
#' centimorgan (cM)). Note that chromosome and position must be divided
#' in order.
#' @param geno matrix. A n*k matrix contains the k markers of the n
#' individuals. The marker genotypes of P1 homozygote (MM),
#' heterozygote (Mm) and P2 homozygote (mm) are coded as 2, 1 and 0,
#' respectively, and NA for missing value.
#' @param y vector. A vector that contains the phenotype values of
#' individuals with genotyped.
#' @param yu vector. A vector that contains the phenotype value
#' of the individuals without genotyped.
#' @param sele.g character. If sele.g="n", it will consider that the
#' data is not a selective genotyping data but the complete genotyping
#' data. If sele.g="p", it will consider that the data is a selective
#' genotyping data, and use the proposed model (Lee 2014) to analyze.
#' If sele.g="t", it will consider that the data is a selective
#' genotyping data, and use the truncated model (Lee 2014) to analyze.
#' If sele.g="f, it will consider that the data is a selective
#' genotyping data, and use the frequency-based model (Lee 2014) to
#' analyze. Note that the yu must be input when sele.g="p" of "f".
#' @param tL numeric. The lower truncation point of phenotype value
#' when sele.g="t". Note that when sele.g="t" and tL=NULL, the yu
#' must be input and the function will consider the minimum of yu
#' as the lower truncation point.
#' @param tR numeric. The upper truncation point of phenotype value
#' when sele.g="t". Note that when sele.g="t" and tR=NULL, the yu
#' must be input and the function will consider the maximum of yu
#' as the upper truncation point.
#' @param method character. method="EM" means the interval mapping method
#' by Lander and Botstein (1989) is used in the analysis, while
#' method="REG" means  the approximate regression interval mapping method
#' by Haley and Knott (1992) is considered in the analysis.
#' @param type character. The population type of the dataset. Include
#' backcross (type="BC"), advanced intercross population (type="AI"), and
#' recombinant inbred population (type="RI").
#' @param D.matrix matrix. The design matrix of the IM model. If
#' D.matrix=NULL, the design matrix will be the constructed using the
#' Cockerham’s model. In BC population, it is a 2*1 matrix which contains
#' 0.5, -0.5 for additive effect. In RI or AI population, it is a 3*2 matrix
#' which contains  1, 0, -1 for additive effect and 0.5, -0.5, 0.5 for
#' dominant effect.
#' @param ng integer. The generation number of the population type. For
#' example, the BC1 population is type="BC" with ng=1; the AI F3
#' population is type="AI" with ng=3.
#' @param cM logical. Specify the unit of marker position. cM=TRUE for
#' centi-Morgan. Or cM=FALSE for Morgan.
#' @param speed numeric. The walking speed of the QTL search (in cM).
#' @param conv numeric. The convergent criterion of EM algorithm.
#' The E and M steps will be iterated until a convergent criterion
#' is satisfied.
#' @param d.eff logical. Specify whether the dominant effect will be
#' considered in the parameter estimation or not for AI or RI population.
#' @param LRT.thre logical or numeric. If being TRUE, the LRT threshold
#' will be computed based on the Gaussian stochastic process
#' (Kao and Ho 2012). Or users can input a numerical value as the LRT
#' threshold to assessing the significance of QTL detection.
#' @param simu integer. To decide how many simulation samples will be used
#' to compute the LRT threshold using the Gaussian process.
#' @param alpha numeric. The type I error rate for the LRT threshold.
#' @param detect logical. Whether the significant QTL whose LRT statistic
#' is larger than the LRT threshold will be shown in the output dataset or
#' not.
#' @param QTLdist numeric. The minimum distance (cM) among different
#' linked significant QTL.
#' @param plot.all logical. If being TRUE, output the profile of LRT
#' statistics for the genome in one figure.
#' @param plot.chr logical. If being TRUE, output the profile of LRT
#' statistics for the chromosomes.
#' @param console logical. To decide whether the process of algorithm will
#' be shown in the R console or not.
#'
#' @return
#' \item{effect}{The estimated effects and LRT statistics of all positions.}
#' \item{LRT.threshold}{The LRT threshold value computed for the data using the
#' Gaussian stochastic process (Kuo 2011; Kao and Ho 2012).}
#' \item{detect.QTL}{The positions, effects and LRT statistics of the detected
#' QTL significant using the obtained LRT threshold value.}
#' \item{model}{The model of selective genotyping data in this analyze.}
#'
#' Graphical outputs including LOD value and effect of each position.
#'
#' @export
#'
#' @references
#'
#' KAO, C.-H. and Z.-B. ZENG 1997 General formulas for obtaining the maximum
#' likelihood estimates and the asymptotic variance-covariance matrix in QTL
#' mapping when using the EM algorithm. Biometrics 53, 653-665.
#'
#' KAO, C.-H., Z.-B. ZENG and R. D. TEASDALE 1999 Multiple interval mapping
#' for Quantitative Trait Loci. Genetics 152: 1203-1216.
#'
#' H.-I LEE, H.-A. HO and C.-H. KAO 2014 A new simple method for improving
#' QTL mapping under selective genotyping. Genetics 198: 1685-1698.
#'
#' KAO, C.-H. and H.-A. Ho 2012 A score-statistic approach for determining
#' threshold values in QTL mapping. Frontiers in Bioscience. E4, 2670-2682.
#'
#' @seealso
#' \code{\link[QTLEMM]{EM.MIM2}}
#' \code{\link[QTLEMM]{IM.search}}
#' \code{\link[QTLEMM]{LRTthre}}
#'
#' @examples
#' # load the example data
#' load(system.file("extdata", "exampledata.RDATA", package = "QTLEMM"))
#'
#' # make the seletive genotyping data
#' ys <- y[y > quantile(y)[4] | y < quantile(y)[2]]
#' yu <- y[y >= quantile(y)[2] & y <= quantile(y)[4]]
#' geno.s <- geno[y > quantile(y)[4] | y < quantile(y)[2],]
#'
#' # run and result
#' result <- IM.search2(marker, geno.s, ys, yu, sele.g = "p", type = "RI", ng = 2,
#' speed = 7.5, conv = 10^-3, LRT.thre = 10)
#' result$detect.QTL
IM.search2 <- function(marker, geno, y, yu = NULL, sele.g = "n", tL = NULL, tR = NULL, method = "EM",
                       type = "RI", D.matrix = NULL, ng = 2, cM = TRUE, speed = 1, conv = 10^-5,
                       d.eff = FALSE, LRT.thre = TRUE, simu = 1000, alpha = 0.05, detect = TRUE,
                       QTLdist = 15, plot.all = TRUE, plot.chr = TRUE, console = TRUE){

  if(is.null(marker) | is.null(geno) | is.null(y)){
    stop("Input data is missing, please cheak and fix", call. = FALSE)
  }

  genotest <- table(geno)
  datatry <- try(geno*geno, silent=TRUE)
  if(class(datatry)[1] == "try-error" | FALSE%in%(names(genotest)%in%c(0, 1, 2))  | !is.matrix(geno)){
    stop("Genotype data error, please cheak your genotype data.", call. = FALSE)
  }

  marker <- as.matrix(marker)
  markertest <- c(ncol(marker) != 2, NA%in%marker, marker[,1] != sort(marker[,1]), nrow(marker) != ncol(geno))
  datatry <- try(marker*marker, silent=TRUE)
  if(class(datatry)[1] == "try-error" | TRUE%in%markertest){
    stop("Marker data error, or the number of marker does not match the genetype data.", call. = FALSE)
  }

  Y <- y
  Y[is.na(Y)] <- mean(Y,na.rm = TRUE)

  if(!is.null(yu)){
    yu <- yu[!is.na(yu)]
    datatry <- try(yu%*%yu, silent=TRUE)
    if(class(datatry)[1] == "try-error"){
      stop("yu data error, please check your yu data.", call. = FALSE)
    }
  }

  datatry <- try(y%*%geno, silent=TRUE)
  if(class(datatry)[1] == "try-error"){
    stop("Phenotype data error, or the number of individual does not match the genetype data.", call. = FALSE)
  }

  if(!sele.g[1] %in% c("n", "t", "p", "f") | length(sele.g)!=1){
    stop("Parameter sele.g error, please check and fix.", call. = FALSE)
  }

  if(sele.g == "t"){
    lrtest <- c(tL[1] < min(c(y,yu)), tR[1] > max(c(y,yu)), tR[1] < tL[1],
                length(tL) > 1, length(tR) > 1)
    datatry <- try(tL[1]*tR[1], silent=TRUE)
    if(class(datatry)[1] == "try-error" | TRUE %in% lrtest){
      stop("Parameter tL or tR error, please check and fix.", call. = FALSE)
    }
    if(is.null(tL) | is.null(tR)){
      if(is.null(yu)){
        stop("yu data error, the yu data must be input for truncated model when parameter tL or tR is not set.", call. = FALSE)
      }
    }
  }

  if(!method[1] %in% c("EM","REG") | length(method) > 1){method <- "EM"}

  if(!type[1] %in% c("AI","RI","BC") | length(type) > 1){
    stop("Parameter type error, please input AI, RI, or BC.", call. = FALSE)
  }

  if(!is.null(D.matrix)){
    if(is.matrix(D.matrix)){
      if(type == "BC"){
        if(FALSE %in% (dim(D.matrix) == c(2,1))){
          stop("Design matrix error, please cheak and fix", call. = FALSE)
        }
      } else {
        if(FALSE %in% (dim(D.matrix) == c(3,2))){
          stop("Design matrix error, please cheak and fix", call. = FALSE)
        }
      }
    } else {
      stop("Design matrix error, please cheak and fix", call. = FALSE)
    }
  }

  if(!is.numeric(ng) | length(ng) > 1 | min(ng) < 1){
    stop("Parameter ng error, please input a positive integer.", call. = FALSE)
  }
  ng <- round(ng)

  if(!cM[1] %in% c(0,1) | length(cM > 1)){cM <- TRUE}

  if(!is.numeric(speed) | length(speed) > 1 | min(speed) < 0){
    stop("Parameter speed error, please input a positive number.", call. = FALSE)
  }

  if(!is.numeric(conv) | length(conv) > 1 | min(conv) < 0){
    stop("Parameter conv error, please input a positive number.", call. = FALSE)
  }

  if(!d.eff[1] %in% c(0,1) | length(d.eff) > 1){d.eff <- TRUE}

  if(min(LRT.thre) < 0 | length(LRT.thre) > 1 | is.character(LRT.thre)){LRT.thre <- TRUE}

  if(!is.numeric(simu) | length(simu) > 1 | min(simu) < 50 | max(simu) > 10^8){
    simu = 1000
  }

  if(!is.numeric(alpha) | length(alpha) > 1 | min(alpha) < 0 | max(alpha) > 1){
    stop("Parameter alpha error, please input a positive number between 0 and 1.", call. = FALSE)
  }

  if(!detect[1] %in% c(0,1) | length(detect > 1)){detect <- TRUE}

  if(!is.numeric(QTLdist) | length(QTLdist) > 1 | min(QTLdist) < speed*2){
    stop("Parameter QTLdist error, please input a bigger positive number.", call. = FALSE)
  }

  if(!plot.all[1] %in% c(0,1) | length(plot.all) > 1){plot.all <- TRUE}
  if(!plot.chr[1] %in% c(0,1) | length(plot.chr) > 1){plot.chr <- TRUE}
  if(!console[1] %in% c(0,1) | length(console) > 1){console <- TRUE}

  if(!cM){
    marker[, 2] <- marker[, 2]*100
  }

  if(is.null(D.matrix)){
    if(type == "BC"){
      D.matrix <- matrix(c(0.5, -0.5), 2, 1)
      row.names(D.matrix) <- c(2, 1)
      colnames(D.matrix) <- "a1"
    } else {
      a1 <- c(1, 0, -1)
      d1 <- c(-0.5, 0.5, -0.5)
      D.matrix <- cbind(a1, d1)
      row.names(D.matrix) <- c(2, 1, 0)
    }
  }
  if(!d.eff & type != "BC"){
    D.matrix <- matrix(D.matrix[, 1], 3, 1)
    row.names(D.matrix) <- c(2, 1, 0)
    colnames(D.matrix) <- "a1"
  }

  if(method == "EM"){
    meth <- function(QTL, marker, geno, D.matrix, y, yu, tL, tR, type, ng, sele.g, conv){
      EM <- EM.MIM2(QTL, marker, geno, D.matrix, y = y, yu = yu, tL = tL, tR = tR, type = type,
                    ng = ng, sele.g = sele.g, conv = conv, console = FALSE)
      eff <- as.numeric(EM$E.vector)
      mu0 <- as.numeric(EM$beta)
      sigma <- sqrt(as.numeric(EM$variance))
      LRT <- EM$LRT
      R2 <- EM$R2
      model <- EM$model
      result <- list(eff, mu0, sigma, LRT, R2, model)
      return(result)
    }
    if(sele.g == "p" | sele.g == "f"){
      Ya <- c(y, yu)
    } else {
      Ya <- y
    }
  } else if (method=="REG"){
    if(sele.g == "p" | sele.g == "f"){
      Ya <- c(y, yu)
    } else {
      Ya <- y
    }

    mixprop <- function(QTL, nu, marker, geno, model, cM = TRUE, type = "RI", ng = 2, cp.matrix = NULL){
      nQTL <- nrow(QTL)
      if(nQTL > 1){QTL <- QTL[order(QTL[, 1], QTL[, 2]),]}
      marker <- marker[order(marker[, 1], marker[, 2]),]

      Q3 <- 3^nQTL
      Q2 <- 2^nQTL
      if(is.null(cp.matrix)){
        cp.matrix <- Q.make(QTL, marker, geno, cM = cM, type = type, ng = ng)[[(nQTL+1)]]
      }
      QTL.freq <- cp.matrix

      if(cM){
        QTL[, 2] <- QTL[, 2]/100
        marker[, 2] <- marker[, 2]/100
      }

      if(model == 2){
        mix.prop <- QTL.freq
        popu.freq <- NULL
        freq.u <- NULL
      } else {
        N <- nrow(QTL.freq)+nu
        freq.s <- colSums(QTL.freq)/N

        if(nQTL == 1){
          if(type == "BC"){
            popu.freq <- c(1-0.5^ng, 0.5^ng)
            Qn <- 2
          } else if (type == "RI"){
            hetro <- 0.5^(ng-1)
            popu.freq <- c((1-hetro)/2, hetro, (1-hetro)/2)
            Qn <- 3
          } else {
            popu.freq <- c(0.25, 0.5, 0.25)
            Qn <- 3
          }
        } else {
          M <- c(1, 0)
          GAM <- matrix(0, Q2, nQTL)
          Pg <- c(2, 1, 0)
          PG <- matrix(0, Q3, nQTL)
          rmn <- rep(0, nQTL-1)
          for(n0 in 1:nQTL){
            GAM[, n0] <- rep(rep(M, 2^(n0-1)), each <- 2^(nQTL-n0))
            PG[, n0] <- rep(rep(Pg, 3^(n0-1)), each <- 3^(nQTL-n0))
            if(n0 == nQTL){break}
            if(QTL[n0, 1]==QTL[(n0+1), 1]){
              rmd <- QTL[(n0+1), 2]-QTL[n0, 2]
              rmn[n0] <- (1-exp(-2*rmd))/2
            }else{
              rmn[n0] <- 0.5
            }
          }

          if(type == "AI" & ng > 2){
            k0 <- seq(0, (ng-1), 2)
            rmn0 <- matrix(0, (nQTL-1), length(k0))
            for(i in 1:length(k0)){
              rmn0[, i] <- choose((ng-1), k0[i])*rmn^k0[i]*(1-rmn)^(ng-1-k0[i])
            }
            rmn <- 1-apply(rmn0, 1, sum)
          }

          gamf <- GAM.f <- matrix(0, Q2, nQTL-1)
          GAM.freq <- rep(1, Q2)
          for(nr in 1:(nQTL-1)){
            gamf[, nr] <- GAM[, nr] == GAM[, nr+1]
            for(ngam in 1:Q2){
              GAM.f[ngam, nr] <- ifelse(gamf[ngam, nr] == 1, 1-rmn[nr], rmn[nr])
            }
            GAM.freq <- GAM.freq*GAM.f[, nr]
          }

          if(type == "BC"){
            popu.freq <- GAM.freq/2
            if(ng > 1){
              G.matrix <- matrix(0, Q2, Q2)
              G.matrix[1, 1] <- 1
              G.matrix[Q2,] <- GAM.freq/2
              for(i in 2:(Q2-1)){
                p0 <- apply(matrix(t(GAM[, which(GAM[i,] != 0)]) == GAM[i, which(GAM[i,] != 0)],
                                   nrow = Q2, byrow = TRUE), 1, sum) == sum(GAM[i,] != 0)
                p1 <- popu.freq[p0]/sum(popu.freq[p0])
                G.matrix[i, p0] <- p1
              }
              for(i in 1:(ng-1)){
                popu.freq <- crossprod(G.matrix, popu.freq)
              }
            }
            Qn <- Q2
          }else{
            Gfreq <- rep(0, (Q2)*(Q2))
            Gtype <- matrix(0, (Q2)*(Q2), nQTL)
            wg <- 0
            for(ga1 in 1:Q2){
              for(ga2 in 1:Q2){
                wg <- wg+1
                Gfreq[wg] <- (GAM.freq[ga1]/2)*(GAM.freq[ga2]/2)
                Gtype[wg,] <- GAM[ga1,]+GAM[ga2,]
              }
            }
            popu.freq <- rep(0, Q3)
            ind <- rep(0, (Q2)*(Q2))
            for(pfi in 1:Q3){
              for(gti in 1:((Q2)*(Q2))){
                if(sum(Gtype[gti,] == PG[pfi,]) == nQTL){ind[gti] <- pfi}
              }
              popu.freq[pfi] <- sum(Gfreq[ind == pfi])
            }

            if(type == "RI" & ng > 2){
              G.matrix <- matrix(0, Q3, Q3)
              G.matrix[1, 1] <- 1
              G.matrix[Q3,Q3] <- 1
              for(i in 2:(Q3-1)){
                p0 <- apply(matrix(t(PG[, which(PG[i,] != 1)]) == PG[i, which(PG[i,] != 1)],
                                   nrow <- Q3, byrow = TRUE), 1, sum) == sum(PG[i,] != 1)
                p1 <- popu.freq[p0]/sum(popu.freq[p0])
                G.matrix[i, p0] <- p1
              }
              for(i in 1:(ng-1)){
                popu.freq <- crossprod(G.matrix, popu.freq)
              }
            }
            Qn <- Q3
          }
        }

        freq.u <- as.vector(popu.freq)-freq.s
        if(model == 3){
          freq.u <- as.vector(popu.freq)
          names(freq.u) <- names(freq.s)
        }
        freq.u[freq.u < 0] <- 0
        Freq.u <- matrix(rep(freq.u/sum(freq.u), each <- nu), nu, Qn)
        mix.prop <- rbind(QTL.freq, Freq.u)
      }
      re <- list(popu.freq = popu.freq, mix.prop = mix.prop, QTL = QTL, marker = marker, freq.u = freq.u)
      return(re)
    }

    meth <- function(QTL, marker, geno, D.matrix, y, yu, tL, tR, type, ng,sele.g, conv){
      mp <- switch(sele.g,
                   p = mixprop(QTL, length(yu), marker, geno, model = 1, cM = cM, type = type, ng = ng)[[2]],
                   t = mixprop(QTL, length(yu), marker, geno, model = 2, cM = cM, type = type, ng = ng)[[2]],
                   f = mixprop(QTL, length(yu), marker, geno, model = 3, cM = cM, type = type, ng = ng)[[2]],
                   n = Q.make(QTL, marker, geno, type = type, ng = ng)$cp.matrix)
      X <- mp%*%D.matrix
      fit <- stats::lm(Ya~X)
      eff <- as.numeric(fit$coefficients[-1])
      mu0 <- as.numeric(fit$coefficients[1])
      ms <- stats::anova(fit)$`Mean Sq`
      sigma <- ms[2]^0.5
      R2 <- summary(fit)$r.squared

      L0 <- c()
      L1 <- c()
      for(k in 1:nrow(mp)){
        L00 <- c()
        L01 <- c()
        for(m in 1:nrow(D.matrix)){
          L00[m] <- mp[k, m]*stats::dnorm((Ya[k]-mu0)/sigma)
          L01[m] <- mp[k, m]*stats::dnorm((Ya[k]-(mu0+D.matrix[m,]%*%eff))/sigma)
        }
        L0[k] <- sum(L00)
        L1[k] <- sum(L01)
      }
      LRT <- -2*sum(log(L0[!is.na(L0) & !is.na(L1)]/L1[!is.na(L0) & !is.na(L1)]))

      result <- list(eff, mu0, sigma, LRT, R2, model = "regression interval mapping model")
      return(result)
    }
  }

  if(console){cat("chr", "cM", "LRT", "\n", sep = "\t")}
  effect <- c()
  cr0 <- unique(marker[, 1])
  for(i in cr0){
    cr <- marker[marker[, 1] == i,]
    for(j in seq(ceiling(floor(min(cr[,2]))+speed),(max(cr[, 2])), speed)){
      QTL <- matrix(c(i, j), 1, 2)

      result <- meth(QTL, marker, geno, D.matrix, y, yu, tL, tR, type, ng, sele.g, conv)
      eff <- result[[1]]
      mu0 <- result[[2]]
      sigma <- result[[3]]
      LRT <- result[[4]]
      R2 <- result[[5]]
      model <- result[[6]]

      eff0 <- c(i, j, eff, LRT, R2)
      effect <- rbind(effect, eff0)
      if(console){cat(i, j, LRT, "\n", sep = "\t")}
    }
  }
  row.names(effect) <- 1:nrow(effect)
  colnames(effect) <- c("chr", "cM", colnames(D.matrix), "LRT", "R2")
  effect <- data.frame(effect)
  effect[effect[, ncol(effect)] == Inf & !is.na(effect[, ncol(effect)]), 3:ncol(effect)] <- 0


  LRT.threshold <- NULL
  if(LRT.thre == TRUE){
    LRT.threshold <- LRTthre(marker, type = type, ng = ng, ns = length(Y), gv = 25,
                             cM = cM, d.eff = d.eff, simu = simu, speed = speed,
                             alpha = alpha, console = console)
  } else if (is.numeric(LRT.thre)){
    LRT.threshold <- LRT.thre
  }

  detectQTL <- function(effect, LRT.threshold, QTLdist = 10){
    LRT <- effect[, c(1, 2, ncol(effect)-1)]
    LRT[LRT[, 3] < LRT.threshold, 3] <- 0
    det0 <- c()
    for(i in unique(LRT[, 1])){
      LRTcr <- LRT[LRT[, 1] == i, 2:3]
      detcr <- rep(0, nrow(LRTcr))
      if(sum(LRTcr[, 2]) > 0){
        for(j in 1:nrow(LRTcr)){
          LRTdet <- LRTcr[LRTcr[,1 ] >= max(c(min(LRTcr[, 1]), LRTcr[j, 1]-QTLdist)) &
                            LRTcr[, 1] <= min(c(max(LRTcr[, 1]), LRTcr[j, 1]+QTLdist)), 2]
          LRT0 <- LRTcr[j, 2]
          if(LRT0 == max(LRTdet) & LRT0 > 0){detcr[j] <- 1
          } else {detcr[j] <- 0}
        }
      }
      det0 <- c(det0, detcr)
    }
    detect.QTL <- effect[det0 == 1,]
    detect.QTL
  }

  detect.QTL <- NULL
  if(detect == TRUE & is.numeric(LRT.threshold)){
    detect.QTL <- detectQTL(effect, LRT.threshold, QTLdist)
  }

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  nc0 <- length(cr0)
  if(plot.chr | (nc0 == 1 & plot.all)){
    for(k in cr0){
      graphics::par(mar = c(5, 5, 4, 2))
      graphics::par(fig = c(0, 1, 0.35, 1))
      plot(effect[effect$chr == k,]$cM, effect[effect$chr == k,]$LRT, main = "", ylab = "LRT statistic",
           xlab = "", type = "l", ylim = c(0, max(c(effect$LRT, LRT.threshold), na.rm = TRUE)), bty = "l")
      if(!is.null(LRT.threshold)){
        graphics::abline(h = LRT.threshold, col = "red", lwd = 2)
      }
      graphics::par(mar = c(2, 5, 4, 2))
      graphics::par(fig = c(0, 1, 0, 0.45), new = TRUE)
      plot(effect[effect$chr == k,]$cM, effect[effect$chr == k,]$a1, ylab = "effect", main = paste("Chromosome", k),
           xlab = "", type = "l", ylim = c(min(c(effect$a1, 0), na.rm = TRUE)*1.4, max(c(effect$a1, 0), na.rm = TRUE)*1.2), xaxt = "n", bty = "n")
      graphics::abline(h=0,col="black",lwd=2)
    }
  }
  if(plot.all & nc0 > 1){
    ncr <- table(effect$chr)
    graphics::par(mar = c(5, 5, 4, 2))
    graphics::par(fig = c(0, 1, 0.3, 1))

    lcr <- c(0, cumsum(ncr))
    x0 <- c()
    cut0 <- (max(lcr)*speed/length(lcr)/5)*(cr0-1)
    xn <- 0
    for(i in 1:length(ncr)){
      x0 <- c(x0, seq(speed, ncr[i]*speed, speed)+lcr[i]*speed+cut0[i])
      xn <- c(xn, length(x0))
    }
    xn1 <- (lcr[-1]-ncr/2)*speed+cut0
    yli <- max(c(effect$LRT, LRT.threshold), na.rm = TRUE)
    plot(x0, effect$LRT, main = "", ylab = "LRT statistic", xlab = "", type = "n",
         ylim = c(-yli/10, yli), xaxt = "n", bty = "l",axes = FALSE)
    for(k in 1:(length(xn)-1)){
      graphics::points(x0[(xn[k]+1):xn[k+1]], effect$LRT[(xn[k]+1):xn[k+1]], main = "", ylab = "LRT statistic",
                       xlab = "", type = "l", ylim = c(-yli/10, yli), xaxt = "n", bty = "l")
    }
    if(!is.null(LRT.threshold)){graphics::abline(h = LRT.threshold, col = "red", lwd = 2)}
    graphics::axis(2, seq(0, yli*1.2, 5))
    graphics::axis(1, cr0, at = xn1, cex = 1.2, tick = FALSE)
    lse <- 4000/max(x0)
    graphics::segments(x0, rep(-yli/10, length(x0)), x0, rep(-yli/5, length(x0)), lwd = lse)

    graphics::par(mar = c(2, 5, 4, 2))
    graphics::par(fig = c(0, 1, 0, 0.5), new = TRUE)
    plot(x0, effect$a1, main = "", ylab = "effect", xlab = "", type = "n",
         ylim = c(min(c(effect$a1, 0), na.rm = TRUE)*1.2, max(c(effect$a1, 0), na.rm = TRUE)*1.2), xaxt = "n", bty = "n")
    graphics::abline(h = 0, col = "blue", lwd = 2)
    for(k in 1:(length(xn)-1)){
      graphics::points(x0[(xn[k]+1):xn[k+1]], effect$a1[(xn[k]+1):xn[k+1]], main = "", ylab = "effect",
                       xlab = "", type = "l", ylim = c(min(c(effect$a1, 0), na.rm = TRUE)*1.2, max(c(effect$a1, 0))*1.2),
                       xaxt = "n", bty = "n")
    }
  }

  result <- list(effect = effect, LRT.threshold = LRT.threshold, detect.QTL = detect.QTL, model = model)
}
