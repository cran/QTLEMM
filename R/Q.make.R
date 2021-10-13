#' Generate Q Matrix
#'
#' Generate the conditional probability matrix by the information of QTL
#' and marker and the genotype data.
#'
#' @param QTL matrix. A q*2 matrix contains the QTL information, where
#' the row dimension q is the number of QTLs in the chromosomes. The
#' first column labels the chromosomes where the QTLs are located, and
#' the second column labels the positions of QTLs (in morgan (M) or
#' centimorgan (cM)). Note that chromosome and position must be divided
#' in order.
#' @param marker matrix. A k*2 matrix contains the marker information,
#' where the row dimension k is the number of markers in the chromosomes.
#' The first column labels the chromosomes where the markers are located,
#' and the second column labels the positions of QTLs (in morgan (M) or
#' centimorgan (cM)). Note that chromosomes and positions must be divided
#' in order.
#' @param geno matrix. A n*k matrix contains the k markers of the n
#' individuals. The marker genotypes of P1 homozygote (MM),
#' heterozygote (Mm) and P2 homozygote (mm) are coded as 2, 1 and 0,
#' respectively, and NA for missing value.
#' @param interval logical. When the QTL with the same position of a
#' marker, whether the marker will be skipped and not be regarded as a
#' flanking marker or not. interval=TRUE presents the marker will be
#' skipped.
#' @param type character. The population type of the dataset. Include
#' backcross (type="BC"), advanced intercross population (type="AI"), and
#' recombinant inbred population (type="RI").
#' @param ng integer. The generation number of the population type. For
#' example, the BC1 population is type="BC" with ng=1; the AI F3
#' population is type="AI" with ng=3.
#' @param cM logical. Specify the unit of marker position. cM=TRUE for
#' centi-Morgan. Or cM=FALSE for Morgan.
#'
#' @return
#' The output contains k conditional probability matrices for the k
#' flanking marker pairs (the k Q-matrices) and a conditional
#' probability matrix of each QTL for all individuals (the cp-matrix,
#' if the genotype data of testing population is input).
#'
#' @note
#' If geno=NULL, the function can be run too and the output will
#' contain k Q-matrices but no cp-matrix.
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
#'
#' @examples
#' # load the example data
#' load(system.file("extdata", "exampledata.RDATA", package = "QTLEMM"))
#'
#' # run and result
#' result <- Q.make(QTL, marker, geno)
#' head(result$cp.matrix)
Q.make <- function(QTL, marker, geno = NULL, interval = FALSE, type = "RI", ng = 2, cM = TRUE){

  if(is.null(QTL) | is.null(marker)){
    stop("Input data is missing, please cheak and fix.", call. = FALSE)
  }

  markertest <- c(ncol(marker) != 2, NA %in% marker, marker[,1] != sort(marker[,1]))

  if(!is.null(geno)){
    genotest <- table(geno)
    datatry <- try(geno*geno, silent=TRUE)
    if(class(datatry)[1] == "try-error" | FALSE %in% (names(genotest) %in% c(0, 1, 2))  | !is.matrix(geno)){
      stop("Genotype data error, please cheak your genotype data.", call. = FALSE)
    }
    markertest <- c(markertest, nrow(marker) != ncol(geno))
  }

  datatry <- try(marker*marker, silent=TRUE)
  if(class(datatry)[1] == "try-error" | TRUE %in% markertest){
    stop("Marker data error, or the number of marker does not match the genetype data.", call. = FALSE)
  }

  datatry <- try(QTL*QTL, silent=TRUE)
  if(class(datatry)[1] == "try-error" | ncol(QTL) != 2 | NA %in% QTL | max(QTL[, 1]) > max(marker[, 1])){
    stop("QTL data error, please cheak your QTL data.", call. = FALSE)
  }

  for(i in 1:nrow(QTL)){
    ch0 <- QTL[i, 1]
    q0 <- QTL[i, 2]
    if(!ch0 %in% marker[, 1]){
      stop("The specified QTL is not in the position that can estimate by the marker data.", call. = FALSE)
    }
    if(q0 > max(marker[marker[, 1] == ch0, 2]) | q0 < min(marker[marker[, 1] == ch0, 2])){
      stop("The specified QTL is not in the position that can estimate by the marker data.", call. = FALSE)
    }
  }

  if(!interval[1] %in% c(0,1) | length(interval) > 1){interval <- FALSE}
  if(!cM[1] %in% c(0,1) | length(cM) > 1){cM <- TRUE}

  if(!type[1] %in% c("AI","RI","BC") | length(type) > 1){
    stop("Parameter type error, please input AI, RI, or BC.", call. = FALSE)
  }

  if(!is.numeric(ng) | length(ng) > 1 | min(ng) < 1){
    stop("Parameter ng error, please input a positive integer.", call. = FALSE)
  }
  ng <- round(ng)

  cr <- QTL[, 1]
  QTL <- QTL[, 2]

  if(cM){
    QTL <- QTL/100
    marker[, 2] <- marker[, 2]/100
  }

  nd <- 3
  if(type == "BC"){
    checkBC <- function(x){
      return(length(table(x)))
    }
    check0 <- apply(geno, 2, checkBC)
    if(3 %in% check0){
      warning("The geno data may not be a back cross geno data. The results may be biased.")
    }
    geno[geno == 0] <- 2

    nd <- 2
    type.fun <- function(d1, d2, ng){
      bcp2M <- function(d, nG){
        r <- (1-exp(-2*d/100))/2
        pAB <- pAb <- paB <- pab <- c()

        pAB[1] <- (1-r)/2
        pAb[1] <- r/2
        paB[1] <- r/2
        pab[1] <- (1-r)/2
        if(nG > 1){
          for(i in 2:nG){
            pAB[i] <- pAB[i-1]+1/2*pAb[i-1]+1/2*paB[i-1]+(1-r)/2*pab[i-1]
            pAb[i] <- 1/2*pAb[i-1]+r/2*pab[i-1]
            paB[i] <- 1/2*paB[i-1]+r/2*pab[i-1]
            pab[i] <- (1-r)/2*pab[i-1]
          }
        }

        c(pAB[nG], pAb[nG], paB[nG], pab[nG])
      }

      bcp3M <- function(d1, d2, nG){
        r1 <- (1-exp(-2*d1/100))/2
        r2 <- (1-exp(-2*d2/100))/2
        pABC <- pABc <- pAbC <- paBC <- pabC <- paBc <- pAbc <- pabc <- c()

        pABC[1] <- (1-r1)*(1-r2)/2
        pABc[1] <- (1-r1)*r2/2
        pAbC[1] <- r1*r2/2
        paBC[1] <- r1*(1-r2)/2
        pabC[1] <- (1-r1)*r2/2
        paBc[1] <- r1*r2/2
        pAbc[1] <- r1*(1-r2)/2
        pabc[1] <- (1-r1)*(1-r2)/2

        if(nG > 1){
          for(i in 2:nG){
            pABC[i] <- pABC[i-1]+1/2*pABc[i-1]+1/2*pAbC[i-1]+1/2*paBC[i-1]+(1-r1)/2*pabC[i-1]+
              ((1-r1)*(1-r2)+r1*r2)/2*paBc[i-1]+(1-r2)/2*pAbc[i-1]+(1-r1)*(1-r2)/2*pabc[i-1]
            pABc[i] <- 1/2*pABc[i-1]+(r1*(1-r2)+(1-r1)*r2)/2*paBc[i-1]+r2/2*pAbc[i-1]+(1-r1)*r2/2*pabc[i-1]
            pAbC[i] <- 1/2*pAbC[i-1]+r1/2*pabC[i-1]+r2/2*pAbc[i-1]+r1*r2/2*pabc[i-1]
            paBC[i] <- 1/2*paBC[i-1]+r1/2*pabC[i-1]+((1-r1)*r2+r1*(1-r2))/2*paBc[i-1]+r1*(1-r2)/2*pabc[i-1]
            pabC[i] <- (1-r1)/2*pabC[i-1]+(1-r1)*r2/2*pabc[i-1]
            paBc[i] <- ((1-r1)*(1-r2)+r1*r2)/2*paBc[i-1]+r1*r2/2*pabc[i-1]
            pAbc[i] <- (1-r2)/2*pAbc[i-1]+r1*(1-r2)/2*pabc[i-1]
            pabc[i] <- (1-r1)*(1-r2)/2*pabc[i-1]
          }
        }
        c(pABC[nG], pABc[nG], pAbC[nG], paBC[nG], pabC[nG], paBc[nG], pAbc[nG], pabc[nG])
      }

      MN <- bcp2M(d1+d2, ng)
      A <- bcp3M(d1, d2, ng)
      Q0 <- c(A[1]/MN[1], A[3]/MN[1], A[2]/MN[2], A[7]/MN[2],
              A[4]/MN[3], A[5]/MN[3], A[6]/MN[4], A[8]/MN[4])
      Q1 <- matrix(Q0, 4, 2, byrow = TRUE)
      colnames(Q1) <- c("QQ", "Qq")
      rownames(Q1) <- c(22, 21, 12, 11)
      return(Q1)
    }
  } else if (type == "RI"){
    type.fun <- function(d1, d2, ng){
      RI2 <- function(d, n.Generation){
        r <- (1 - exp(-2 * d))/2
        TransMatrix <- matrix(1:25, nrow = 5)
        TransMatrix[1, ] <- c(4, 0, 2, (1-r)^2, r^2)/4
        TransMatrix[2, ] <- c(0, 4, 2, r^2, (1-r)^2)/4
        TransMatrix[3, ] <- c(0, 0, 1, r*(1-r), r*(1-r))/2
        TransMatrix[4, ] <- c(0, 0, 0, (1-r)^2, r^2)/2
        TransMatrix[5, ] <- c(0, 0, 0, r^2, ( -r)^2)/2
        Freq <- c(0, 0, 0, 1, 0)
        n.Iteration <- n.Generation - 1
        ZygoteFreq <- matrix(rep(0, (9 * n.Generation - 9)), nrow = 9)
        for (i in 1:n.Iteration) {
          Freq <- TransMatrix %*% Freq
          ZygoteFreq[, i] <- c(Freq[c(1, 3, 2, 3)], sum(Freq[4:5]),
                               Freq[c(3, 2, 3, 1)])
        }
        ZygoteFreq[, n.Iteration]
      }

      RI3 <- function(d1, d2, n.Generation){
        r1 <- (1 - exp(-2 * d1))/2
        r2 <- (1 - exp(-2 * d2))/2
        r <- (1 - exp(-2 * (d1 + d2)))/2
        f <- r1 * (1 - r2) + r2 * (1 - r1)
        TransMatrix <- matrix(rep(0, 20^2), nrow = 20)
        TransMatrix[1, ] <- c(4, 1, 0, 1, (1-r2)^2, r2^2, 0, 0, 0, 0, 1, (1-r1)^2, r1^2, 0, (1-f)^2,
                              f^2, ((1-r1) * (1-r2))^2, (r1*r2)^2, ((1-r1)*r2)^2, (r1*(1-r2))^2)/4
        TransMatrix[2, ] <- c(0, 1, 0, 0, r2*(1-r2), r2*(1-r2), 0, 0, 0, 0, 0, r1*(1-r1), r1*(1-r1),
                              0, 0, 0, rep(r1*(1-r1)*r2*(1-r2), 4))/2
        TransMatrix[3, ] <- c(0, 1, 4, 0, r2^2, (1-r2)^2, 1, 0, 0, 0, 0, r1^2, (1-r1)^2, 1, (1-f)^2,
                              f^2, (r1*r2)^2, ((1-r1)*(1-r2))^2, (r1*(1-r2))^2, ((1-r1)*r2)^2)/4
        TransMatrix[4, ] <- c(0, 0, 0, 1, r2*(1-r2), r2*(1-r2), rep(0, 8), f*(1-f), f*(1-f),
                              (1-r1)^2*r2*(1-r2), r1^2*r2*(1-r2), (1-r1)^2*r2*(1-r2), r1^2*r2*(1-r2))/2
        TransMatrix[5, ] <- c(0, 0, 0, 0, (1-r2)^2, r2^2, rep(0, 10), r1*(1-r1)*(1-r2)^2,
                              r1*(1-r1)*r2^2, r1*(1-r1)*r2^2, r1*(1-r1)*(1-r2)^2)/2
        TransMatrix[6, ] <- c(0, 0, 0, 0, r2^2, (1-r2)^2, rep(0, 10), r1*(1-r1)*r2^2,
                              r1*(1-r1)*(1-r2)^2, r1*(1-r1)*(1-r2)^2, r1*(1-r1)*r2^2)/2
        TransMatrix[7, ] <- c(0, 0, 0, 0, r2*(1-r2), r2*(1-r2), 1, rep(0, 7), f*(1-f), f*(1-f),
                              r1^2*r2*(1-r2), (1-r1)^2*r2*(1-r2), r1^2*r2*(1-r2), (1-r1)^2*r2*(1-r2))/2
        TransMatrix[8, ] <- c(0, 0, 0, 1, r2^2, (1-r2)^2, 0, 4, 1, 0, 0, (1-r1)^2, r1^2, 1, f^2,
                              (1-f)^2, ((1-r1)*r2)^2, (r1*(1-r2))^2, ((1-r1)*(1-r2))^2, (r1*r2)^2)/4
        TransMatrix[9, ] <- c(0, 0, 0, 0, r2*(1-r2), r2*(1-r2), 0, 0, 1, 0, 0, r1*(1-r1),
                              r1*(1-r1), 0, 0, 0, rep(r1*(1-r1)*r2*(1-r2), 4))/2
        TransMatrix[10, ] <- c(0, 0, 0, 0, (1-r2)^2, r2^2, 1, 0, 1, 4, 1, r1^2, (1-r1)^2, 0, f^2,
                               (1-f)^2, (r1*(1-r2))^2, ((1-r1)*r2)^2, (r1*r2)^2, ((1-r1)*(1-r2))^2)/4
        TransMatrix[11, ] <- c(rep(0, 10), 1, r1 * (1-r1), r1*(1-r1), 0, f*(1-f), f*(1-f),
                               r1*(1-r1)*(1-r2)^2, r1*(1-r1)*r2^2, r1*(1-r1)*r2^2, r1*(1-r1)*(1-r2)^2)/2
        TransMatrix[12, ] <- c(rep(0, 11), (1-r1)^2, r1^2, 0, 0, 0, (1-r1)^2*r2*(1-r2),
                               r1^2*r2*(1-r2), (1-r1)^2*r2*(1-r2), r1^2*r2*(1-r2))/2
        TransMatrix[13, ] <- c(rep(0, 11), r1^2, (1-r1)^2, 0, 0, 0, r1^2*r2*(1-r2), (1-r1)^2*r2*(1-r2),
                               r1^2*r2*(1-r2), (1-r1)^2*r2*(1-r2))/2
        TransMatrix[14, ] <- c(rep(0, 11), r1*(1-r1), r1*(1-r1), 1, f*(1-f), f*(1-f), r1*(1-r1)*r2^2,
                               r1*(1-r1)*(1-r2)^2, r1*(1-r1)*(1-r2)^2, r1*(1-r1)*r2^2)/2
        TransMatrix[15, ] <- c(rep(0, 14), (1-f)^2, f^2, rep(r1*(1-r1)*r2*(1-r2), 4))/2
        TransMatrix[16, ] <- c(rep(0, 14), f^2, (1-f)^2, rep(r1*(1-r1)*r2*(1-r2), 4))/2
        TransMatrix[17, ] <- c(rep(0, 16), ((1-r1)*(1-r2))^2, (r1*r2)^2, ((1-r1)*r2)^2, (r1*(1-r2))^2)/2
        TransMatrix[18, ] <- c(rep(0, 16), (r1*r2)^2, ((1-r1)*(1-r2))^2, (r1*(1-r2))^2, ((1-r1)*r2)^2)/2
        TransMatrix[19, ] <- c(rep(0, 16), ((1-r1)*r2)^2, (r1*(1-r2))^2, ((1-r1)*(1-r2))^2, (r1*r2)^2)/2
        TransMatrix[20, ] <- c(rep(0, 16), (r1*(1-r2))^2, ((1-r1)*r2)^2, (r1*r2)^2, ((1-r1)*(1-r2))^2)/2

        Freq <- c(rep(0, 16), 1, 0, 0, 0)
        n.Iteration <- n.Generation - 1
        ZygoteFreq <- matrix(rep(0, (20*n.Generation - 20)), nrow = 20)
        for (i in 1:n.Iteration) {
          Freq <- TransMatrix%*%Freq
          ZygoteFreq[, i] <- Freq
        }
        ZygoteFreq[, n.Iteration]
      }

      MN <- RI2(d1+d2, ng)
      A <- RI3(d1, d2, ng)
      Q0 <- c(A[1:3]/MN[1], c(A[4], A[5]+A[6], A[7])/MN[2], A[8:10]/MN[3],
              c(A[11], A[12]+A[13], A[14])/MN[4],
              c(A[15]+A[16], sum(A[17:20]), A[15]+A[16])/MN[5],
              c(A[14], A[12]+A[13], A[11])/MN[6],
              A[10:8]/MN[7], c(A[7], A[5]+A[6], A[4])/MN[8], A[3:1]/MN[9])
      Q1 <- matrix(Q0, 9, 3, byrow = TRUE)
      colnames(Q1) <- c("QQ", "Qq", "qq")
      rownames(Q1) <- c(22, 21, 20, 12, 11, 10, "02", "01", "00")
      return(Q1)
    }
  } else if (type=="AI"){
    type.fun <- function(d1, d2, ng){
      AI2 <- function(d, n.Generation){
        r <- (1 - exp(-2*d))/2
        ZygoteFreq <- matrix(rep(0, (9*n.Generation-9)), nrow = 9)
        GameteFreq <- c((1-r)/2, r/2)
        C <- GameteFreq[1]^2
        D <- GameteFreq[2]^2
        E <- 2*GameteFreq[1]*GameteFreq[2]
        F0 <- 2*GameteFreq[1]^2
        G <- 2*GameteFreq[2]^2
        ZygoteFreq[, 1] <- c(C, E, D, E, F0+G, E, D, E, C)
        if (n.Generation > 2) {
          for (i in 3:(n.Generation)) {
            GameteFreq <- (1-r)*GameteFreq+r/4
            C <- GameteFreq[1]^2
            D <- GameteFreq[2]^2
            E <- 2*GameteFreq[1]*GameteFreq[2]
            F0 <- 2*GameteFreq[1]^2
            G <- 2*GameteFreq[2]^2
            ZygoteFreq[, (i-1)] <- c(C, E, D, E, F0+G, E, D, E, C)
          }
        }
        ZygoteFreq[, (n.Generation - 1)]
      }

      AI3=function(d1, d2, n.Generation){
        d <- d1+d2
        r1 <- (1-exp(-2*d1))/2
        r2 <- (1-exp(-2*d2))/2
        r <- (1-exp(-2*d))/2
        ZygoteFreq <- matrix(rep(0, (20*n.Generation-20)), nrow = 20)
        N <- (1-r1)*(1-r2)/2
        R <- (1-r1)*r2/2
        B <- r1*r2/2
        L <- r1*(1-r2)/2
        COEF <- rep(2, 20)
        COEF[c(1, 3, 8, 10)] <- rep(1, 4)
        GS <- c(N, N, B, N, N, R, B, R, R, L, N, N, B, R, N, R,
                N, B, R, L)
        GE <- c(N, B, B, R, L, B, L, R, L, L, L, R, L, B, B, L,
                N, B, R, L)
        ZygoteFreq[, 1] <- COEF*GS*GE
        if (n.Generation > 2) {
          for (i in 3:n.Generation) {
            F23 <- c(N+L, B+R, B+R, N+L)
            F13 <- rep(c(N+B, L+R), 2)
            F12 <- rep(c(N+R, L+B), each = 2)
            NewGameteFreq <- (1-r1)*(1-r2)*c(N, R, B, L)+0.5*r1*(1-r2)*F23+
              0.5*r1*r2*F13+0.5*(1-r1)*r2*F12
            N <- NewGameteFreq[1]
            R <- NewGameteFreq[2]
            B <- NewGameteFreq[3]
            L <- NewGameteFreq[4]
            GS <- c(N, N, B, N, N, R, B, R, R, L, N, N, B, R,
                    N, R, N, B, R, L)
            GE <- c(N, B, B, R, L, B, L, R, L, L, L, R, L, B,
                    B, L, N, B, R, L)
            ZygoteFreq[, (i-1)] <- COEF*GS*GE
          }
        }
        ZygoteFreq[, (n.Generation-1)]
      }

      MN <- AI2(d1+d2, ng)
      A <- AI3(d1, d2, ng)
      Q0 <- c(A[1:3]/MN[1], c(A[4], A[5]+A[6], A[7])/MN[2], A[8:10]/MN[3],
              c(A[11], A[12]+A[13], A[14])/MN[4],
              c(A[15]+A[16], sum(A[17:20]), A[15]+A[16])/MN[5],
              c(A[14], A[12]+A[13], A[11])/MN[6],
              A[10:8]/MN[7], c(A[7], A[5]+A[6], A[4])/MN[8], A[3:1]/MN[9])
      Q1 <- matrix(Q0, 9, 3, byrow = TRUE)
      colnames(Q1) <- c("QQ", "Qq", "qq")
      rownames(Q1) <- c(22, 21, 20, 12, 11, 10, "02", "01", "00")
      return(Q1)
    }
  }

  marker <- cbind(marker, 1:nrow(marker))

  cr0 <- marker[marker[, 1]==cr[1], ]

  marker1 <- cr0[max(which(cr0[, 2] < QTL[1])), ]
  marker2 <- cr0[max(which(cr0[, 2] < QTL[1]))+1, ]
  if(interval == T){
    marker2 <- cr0[min(which(cr0[, 2]>QTL[1])), ]
  }
  marker0 <- as.numeric(c(marker1[3], marker2[3]))

  d1 <- as.numeric(QTL[1]-marker1[2])
  d2 <- as.numeric(marker2[2]-QTL[1])
  Q1 <- type.fun(d1, d2, ng)

  Q1 <- list(Q1)
  m <- 1

  n <- length(QTL)
  if(n > 1){
    for(m in 2:n){
      cr0 <- marker[marker[, 1] == cr[m], ]

      marker1 <- cr0[max(which(cr0[, 2] < QTL[m])), ]
      marker2 <- cr0[max(which(cr0[, 2] < QTL[m]))+1, ]
      if(interval == T){
        marker2 <- cr0[min(which(cr0[, 2] > QTL[m])), ]
      }
      marker0 <- as.numeric(c(marker0, marker1[3], marker2[3]))

      d1 <- as.numeric(QTL[m]-marker1[2])
      d2 <- as.numeric(marker2[2]-QTL[m])
      Q2 <- type.fun(d1, d2, ng)
      Q1[[m]] <- Q2
    }
  }


  if(!is.null(geno)){
    red.genotype <- geno[, marker0]
    cp.matrix <- matrix(0, nrow(geno), nd^n)
    D2 <- matrix(0, nd^n, n)
    for(i in 1:n){
      D2[, i] <- rep(rep(2:(3-nd), each = nd^(n-i)), nd^(i-1))
    }
    qname <- apply(D2, 1, paste, collapse = "")
    for(j in 1:nrow(geno)){
      geno.j <- red.genotype[j, ]
      pq <- c()
      for(k in 1:length(QTL)){
        g0 <- geno.j[(k*2-1):(k*2)]
        if(is.na(g0[1]) & is.na(g0[2])){
          a <- Q1[[k]]
          a <- matrix(unlist(a), nrow(a), ncol(a))
          pq0 <- apply(a, 2, mean)
        } else if (!is.na(g0[1]) & is.na(g0[2])){
          a <- row.names(Q1[[k]])
          a <- as.numeric(unlist(strsplit(a,split = "", fixed = TRUE)) )
          a <- t(matrix(a, 2))
          a <- Q1[[k]][a[, 1] == as.numeric(g0[1]), ]
          a <- matrix(unlist(a), nrow(a), ncol(a))
          pq0 <- apply(a, 2, mean)
        } else if (is.na(g0[1]) & !is.na(g0[2])){
          a <- row.names(Q1[[k]])
          a <- as.numeric(unlist(strsplit(a,split = "", fixed = TRUE)) )
          a <- t(matrix(a, 2))
          a <- Q1[[k]][a[, 2] == as.numeric(g0[2]), ]
          a <- matrix(unlist(a), nrow(a), ncol(a))
          pq0 <- apply(a, 2, mean)
        } else {
          g0 <- paste(g0[1], g0[2], sep = "")
          pq0 <- Q1[[k]]
          pq0 <- pq0[row.names(pq0) == g0]
        }
        pq <- rbind(pq, pq0)
      }

      for(i in 1:nd^n){
        a <- qname[i]
        a <- -as.numeric(unlist(strsplit(a,split = "", fixed = TRUE)))+3
        pq1 <- pq[1:length(QTL), a]
        if(is.matrix(pq1)){cp.matrix[j, i] <- prod(as.numeric(diag(pq1)))
        } else {cp.matrix[j,i] <- as.numeric(pq1)}
      }
    }
    colnames(cp.matrix) <- qname
    Q1[[m+1]] <- cp.matrix
  }

  names(Q1)[1:m] <- c(paste("Q.matrix", 1:m, sep = ""))
  if(length(Q1) == m+1){names(Q1)[m+1] <- "cp.matrix"}
  return(Q1)
}
