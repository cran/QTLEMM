#' EM Algorithm for QTL MIM
#'
#' Expectation-maximization algorithm for QTL multiple interval mapping.
#'
#' @param D.matrix matrix. The design matrix of QTL effects is a g*p
#' matrix, where g is the number of possible QTL genotypes, and p is the
#' number of effects considered in the MIM model. This design matrix can
#' be conveniently generated using the function D.make().
#' @param cp.matrix matrix. The conditional probability matrix is an n*g
#' matrix, where n is the number of individuals, and g is the number of
#' possible genotypes of QTLs. This conditional probability matrix can
#' be easily generated using the function Q.make().
#' @param y vector. A vector with n elements that contain the phenotype
#' values of individuals.
#' @param E.vector0 vector. The initial value for QTL effects. The
#' number of elements corresponds to the column dimension of the design
#' matrix. If E.vector0=NULL, the initial value for all effects will be
#' set to 0.
#' @param X matrix. The design matrix of the fixed factors except for
#' QTL effects. It is an n*k matrix, where n is the number of
#' individuals, and k is the number of fixed factors. If X=NULL,
#' the matrix will be an n*1 matrix where all elements are 1.
#' @param beta0 vector. The initial value for effects of the fixed
#' factors. The number of elements corresponds to the column dimension
#' of the fixed factor design matrix.  If beta0=NULL, the initial value
#' will be set to the average of y.
#' @param variance0 numeric. The initial value for variance. If
#' variance0=NULL, the initial value will be set to the variance of
#' phenotype values.
#' @param crit numeric. The convergence criterion of EM algorithm.
#' The E and M steps will iterate until a convergence criterion is met.
#' It must be a value between 0 and 1.
#' @param stop numeric. The stopping criterion of EM algorithm. The E and
#' M steps will halt when the iteration number reaches the stopping
#' criterion, treating the algorithm as having failed to converge.
#' @param conv logical. If set to False, it will disregard the failure to
#' converge and output the last result obtained during the EM algorithm
#' before reaching the stopping criterion.
#' @param console logical. Determines whether the process of the algorithm
#' will be displayed in the R console or not.
#'
#' @return
#' \item{E.vector}{The QTL effects are calculated by the EM algorithm.}
#' \item{beta}{The effects of the fixed factors are calculated by the EM
#' algorithm.}
#' \item{variance}{The error variance is calculated by the EM algorithm.}
#' \item{PI.matrix}{The posterior probabilities matrix after the
#' process of the EM algorithm.}
#' \item{log.likelihood}{The log-likelihood value of this model.}
#' \item{LRT}{The LRT statistic of this model.}
#' \item{R2}{The coefficient of determination of this model. This
#' can be used as an estimate of heritability.}
#' \item{y.hat}{The fitted values of trait values are calculated by
#' the estimated values from the EM algorithm.}
#' \item{iteration.number}{The iteration number of the EM algorithm.}
#'
#' @export
#'
#' @references
#'
#' KAO, C.-H. and Z.-B. ZENG 1997 General formulas for obtaining the maximum
#' likelihood estimates and the asymptotic variance-covariance matrix in QTL
#' mapping when using the EM algorithm. Biometrics 53, 653-665. <doi: 10.2307/2533965.>
#'
#' KAO, C.-H., Z.-B. ZENG and R. D. TEASDALE 1999 Multiple interval mapping
#' for Quantitative Trait Loci. Genetics 152: 1203-1216. <doi: 10.1093/genetics/152.3.1203>
#'
#' @seealso
#' \code{\link[QTLEMM]{D.make}}
#' \code{\link[QTLEMM]{Q.make}}
#' \code{\link[QTLEMM]{EM.MIM2}}
#' \code{\link[QTLEMM]{EM.MIMv}}
#'
#' @examples
#' # load the example data
#' load(system.file("extdata", "exampledata.RDATA", package = "QTLEMM"))
#'
#' # run and result
#' D.matrix <- D.make(3, type = "RI", aa = c(1, 3, 2, 3), dd = c(1, 2, 1, 3), ad = c(1, 2, 2, 3))
#' cp.matrix <- Q.make(QTL, marker, geno, type = "RI", ng = 2)$cp.matrix
#' result <- EM.MIM(D.matrix, cp.matrix, y)
#' result$E.vector
EM.MIM <- function(D.matrix, cp.matrix, y, E.vector0 = NULL, X = NULL, beta0 = NULL,
                   variance0 = NULL, crit = 10^-5, stop = 1000, conv = TRUE, console = TRUE){

  if(is.null(D.matrix) | is.null(cp.matrix) | is.null(y)){
    stop("Input data is missing, please cheak and fix", call. = FALSE)
  }

  datatry <- try(t(y)%*%cp.matrix%*%D.matrix, silent=TRUE)
  if(class(datatry)[1] == "try-error" | NA %in% D.matrix | NA %in% c(cp.matrix)){
    stop("Input data error, please check your input data.", call. = FALSE)
  }

  datatry <- try(y%*%t(y), silent=TRUE)
  if(class(datatry)[1] == "try-error" | length(y) < 2){
    stop("Input data error, please check your input data.", call. = FALSE)
  }

  Y <- c(y)
  ind <- length(Y)
  g <- nrow(D.matrix)
  eff <- ncol(D.matrix)

  Y[is.na(Y)] <- mean(Y,na.rm = TRUE)

  E.vector <- E.vector0
  beta <- beta0
  variance <- variance0

  if(is.null(E.vector)){E.vector <- rep(0, eff)}
  if(is.null(X)){
    X <- matrix(1, ind, 1)
  } else if (is.vector(X)){
    X <- matrix(X, length(X), 1)
  }
  if(is.null(beta)){
    beta <- matrix(rep(mean(Y), ncol(X)), ncol(X), 1)
  } else if (is.numeric(beta)){
    beta <- matrix(rep(beta, ncol(X)), ncol(X), 1)
  }
  if(is.null(variance)){variance <- stats::var(c(Y))}
  if(!console[1] %in% c(0,1) | length(console) > 1){console <- TRUE}
  if(!conv[1] %in% c(0,1) | length(conv) > 1){conv <- TRUE}

  datatry <- try(D.matrix%*%E.vector, silent=TRUE)
  if(class(datatry)[1] == "try-error" | NA %in% E.vector){
    stop("Parameter E.vector0 error, please check and fix.", call. = FALSE)
  }

  datatry <- try(Y%*%X%*%beta, silent=TRUE)
  if(class(datatry)[1] == "try-error" | NA %in% X | NA %in% beta){
    stop("Parameter X or bata0 error, please check and fix.", call. = FALSE)
  }

  if(!is.numeric(variance) | length(variance) > 1 | min(variance) < 0){
    stop("Parameter variance0 error, please input a positive number.", call. = FALSE)
  }
  sigma <- sqrt(variance)

  if(!is.numeric(crit) | length(crit) > 1 | min(crit) <= 0 | max(crit) >= 1){
    stop("Parameter crit error, please input a positive number between 0 and 1.", call. = FALSE)
  }

  if(!is.numeric(stop) | length(stop) > 1 | min(crit) <= 0){
    stop = 1000
    warning("Parameter stop error, adjust to 1000.")
  }

  Delta <- 1
  number <- 0

  if(length(colnames(D.matrix)) == ncol(D.matrix)){
    effectname <- colnames(D.matrix)
  }

  cat(paste("number", "var", effectname, "\n", sep = "\t")[console])

  Yt <- as.matrix(Y)
  indvec <- matrix(1, 1, ind)
  gvec <- matrix(1, 1, g)

  while (max(abs(Delta)) > crit & number < stop) {

    Et <- as.matrix(E.vector)
    bt<- as.matrix(beta)

    muji.matrix <- t(D.matrix%*%E.vector%*%indvec)+X%*%beta%*%gvec

    P0.matrix <- cp.matrix*stats::dnorm(y, muji.matrix, c(sigma))
    P0s <- apply(P0.matrix, 1, sum)
    P0.matrix[P0s == 0, ] <- rep(1, g)
    PIt <- P0.matrix*(1/P0s%*%gvec)

    PD <- t(Yt-X%*%bt)%*%PIt%*%D.matrix
    PDD <- indvec%*%PIt%*%(D.matrix^2)
    r.vector <- c(PD/PDD)

    M.matrix <- matrix(0, eff, eff)
    V.matrix <- matrix(0, eff, eff)
    for(i in 1:eff){
      Di <- D.matrix[, i]
      Dij <- D.matrix
      for(j in 1:eff){
        Dij[, j] <- Di*D.matrix[, j]
      }
      iPI <- indvec%*%PIt
      M01 <- c(iPI)%*%Dij
      V.matrix[i, ] <- c(M01)
      M.matrix[, i] <- c(M01)/c(PDD)
    }
    diag(M.matrix) <- 0

    E.t <- r.vector-M.matrix%*%E.vector
    beta.t <- solve(t(X)%*%X)%*%t(X)%*%(c(Yt)-PIt%*%D.matrix%*%E.t)

    YXb <- Y-X%*%c(beta.t)
    EVE <- t(E.t)%*%V.matrix%*%E.t
    sigma.t <- sqrt((t(YXb)%*%(YXb)-t(YXb)%*%PIt%*%D.matrix%*%E.t*2+EVE)/ind)

    Delta <- c(E.t-E.vector)
    if(NaN %in% Delta){
      break()
    }
    number <- number+1
    Ep <- round(c(E.t), 3)
    sp <- round(c(sigma.t)^2, 3)
    cat(paste(number, sp, Ep, "\n", sep = "\t")[console])

    E.vector <- c(E.t)
    beta <- c(beta.t)
    sigma <- c(sigma.t)

  }

  muji.matrix <- t(D.matrix%*%E.vector%*%indvec)+X%*%beta%*%gvec
  P0.matrix <- cp.matrix*stats::dnorm(y, muji.matrix, c(sigma))
  P0s <- apply(P0.matrix, 1, sum)
  P0.matrix[P0s == 0, ] <- rep(1, g)
  PI.matrix <- P0.matrix*(1/P0s%*%gvec)

  names(E.vector) <- effectname
  colnames(PI.matrix) <- colnames(cp.matrix)

  variance <- sigma^2

  L0 <- rep(0, ind)
  L1 <- rep(0, ind)
  Xb <- X%*%beta
  for(m in 1:g){
    L0 <- L0+(cp.matrix[, m]*stats::dnorm(Y, mean(Xb), sigma))
    L1 <- L1+(cp.matrix[, m]*stats::dnorm(Y, mean(Xb)+D.matrix[m,]%*%E.vector, sigma))
  }
  like0 <- sum(log(L0))
  like1 <- sum(log(L1))
  LRT <- 2*(like1-like0)

  y.hat <- PI.matrix%*%D.matrix%*%E.vector+Xb
  r2 <- c(stats::var(y.hat)/stats::var(y))

  if(number == stop){
    if(conv){
      E.vector <- rep(0, length(E.vector))
      beta <- 0
      variance <- 0
      PI.matrix <- matrix(0, nrow(PI.matrix), ncol(PI.matrix))
      like1 <- -Inf
      LRT <- 0
      r2 <- 0
    }
    warning("EM algorithm fails to converge, please check the input data or adjust
            the convergence criterion and stopping criterion.")
  }

  result <- list(E.vector = E.vector, beta = as.numeric(beta), variance = as.numeric(variance),
                 PI.matrix = PI.matrix, log.likelihood = like1, LRT = LRT, R2 = r2,
                 y.hat = y.hat, iteration.number = number)

  return(result)
}
