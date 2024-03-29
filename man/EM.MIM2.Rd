% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM.MIM2.R
\name{EM.MIM2}
\alias{EM.MIM2}
\title{EM Algorithm for QTL MIM with Selective Genotyping}
\usage{
EM.MIM2(
  QTL,
  marker,
  geno,
  D.matrix,
  cp.matrix = NULL,
  y,
  yu = NULL,
  sele.g = "n",
  tL = NULL,
  tR = NULL,
  type = "RI",
  ng = 2,
  cM = TRUE,
  E.vector0 = NULL,
  X = NULL,
  beta0 = NULL,
  variance0 = NULL,
  crit = 10^-5,
  stop = 1000,
  conv = TRUE,
  console = TRUE
)
}
\arguments{
\item{QTL}{matrix. A q*2 matrix contains the QTL information, where the
row dimension 'q' represents the number of QTLs in the chromosomes. The
first column labels the chromosomes where the QTLs are located, and the
second column labels the positions of QTLs (in morgan (M) or centimorgan
(cM)).}

\item{marker}{matrix. A k*2 matrix contains the marker information,
where the row dimension 'k' represents the number of markers in the
chromosomes. The first column labels the chromosomes where the markers
are located, and the second column labels the positions of markers (in
morgan (M) or centimorgan (cM)). It's important to note that chromosomes
and positions must be sorted in order.}

\item{geno}{matrix. A n*k matrix contains the genotypes of k markers
for n individuals. The marker genotypes of P1 homozygote (MM),
heterozygote (Mm), and P2 homozygote (mm) are coded as 2, 1, and 0,
respectively, with NA indicating missing values.}

\item{D.matrix}{matrix. The design matrix of QTL effects is a g*p
matrix, where g is the number of possible QTL genotypes, and p is the
number of effects considered in the MIM model. This design matrix can
be conveniently generated using the function D.make().}

\item{cp.matrix}{matrix. The conditional probability matrix is an
n*g matrix, where n is the number of genotyped individuals, and g is
the number of possible genotypes of QTLs. If cp.matrix=NULL, the
function will calculate the conditional probability matrix for selective
genotyping.}

\item{y}{vector. A vector that contains the phenotype values of
individuals with genotypes.}

\item{yu}{vector. A vector that contains the phenotype values of
individuals without genotypes.}

\item{sele.g}{character. Determines the type of data being analyzed:
If sele.g="n", it considers the data as complete genotyping data. If
sele.g="f", it treats the data as selective genotyping data and utilizes
the proposed corrected frequency model (Lee 2014) for analysis; If
sele.g="t", it considers the data as selective genotyping data and uses
the truncated model (Lee 2014) for analysis; If sele.g="p", it treats
the data as selective genotyping data and uses the population
frequency-based model (Lee 2014) for analysis. Note that the 'yu'
argument must be provided when sele.g="f" or "p".}

\item{tL}{numeric. The lower truncation point of phenotype value when
sele.g="t". When sele.g="t" and tL=NULL, the 'yu' argument must be
provided. In this case, the function will consider the minimum of 'yu'
as the lower truncation point.}

\item{tR}{numeric. The upper truncation point of phenotype value when
sele.g="t". When sele.g="t" and tR=NULL, the 'yu' argument must be
provided. In this case, the function will consider the maximum of 'yu'
as the upper truncation point.}

\item{type}{character. The population type of the dataset. Includes
backcross (type="BC"), advanced intercross population (type="AI"), and
recombinant inbred population (type="RI"). The default value is "RI".}

\item{ng}{integer. The generation number of the population type. For
instance, in a BC1 population where type="BC", ng=1; in an AI F3
population where type="AI", ng=3.}

\item{cM}{logical. Specify the unit of marker position. If cM=TRUE, it
denotes centimorgan; if cM=FALSE, it denotes morgan.}

\item{E.vector0}{vector. The initial value for QTL effects. The
number of elements corresponds to the column dimension of the design
matrix. If E.vector0=NULL, the initial value for all effects will be
set to 0.}

\item{X}{matrix. The design matrix of the fixed factors except for
QTL effects. It is an n*k matrix, where n is the number of
individuals, and k is the number of fixed factors. If X=NULL,
the matrix will be an n*1 matrix where all elements are 1.}

\item{beta0}{vector. The initial value for effects of the fixed
factors. The number of elements corresponds to the column dimension
of the fixed factor design matrix.  If beta0=NULL, the initial value
will be set to the average of y.}

\item{variance0}{numeric. The initial value for variance. If
variance0=NULL, the initial value will be set to the variance of
phenotype values.}

\item{crit}{numeric. The convergence criterion of EM algorithm.
The E and M steps will iterate until a convergence criterion is met.
It must be a value between 0 and 1.}

\item{stop}{numeric. The stopping criterion of EM algorithm. The E and
M steps will halt when the iteration number reaches the stopping
criterion, treating the algorithm as having failed to converge.}

\item{conv}{logical. If set to False, it will disregard the failure to
converge and output the last result obtained during the EM algorithm
before reaching the stopping criterion.}

\item{console}{logical. Determines whether the process of the algorithm
will be displayed in the R console or not.}
}
\value{
\item{QTL}{The QTL imformation of this analysis.}
\item{E.vector}{The QTL effects are calculated by the EM algorithm.}
\item{beta}{The effects of the fixed factors are calculated by the EM
algorithm.}
\item{variance}{The variance is calculated by the EM algorithm.}
\item{PI.matrix}{The posterior probabilities matrix after the
process of the EM algorithm.}
\item{log.likelihood}{The log-likelihood value of this model.}
\item{LRT}{The LRT statistic of this model.}
\item{R2}{The coefficient of determination of this model. This
can be used as an estimate of heritability.}
\item{y.hat}{The fitted values of trait values with genotyping are
calculated by the estimated values from the EM algorithm.}
\item{yu.hat}{The fitted values of trait values without genotyping
are calculated by the estimated values from the EM algorithm.}
\item{iteration.number}{The iteration number of the EM algorithm.}
\item{model}{The model of this analysis, contains complete a
genotyping model, a proposed model, a truncated model, and a
frequency-based model.}
}
\description{
Expectation-maximization algorithm for QTL multiple interval mapping.
It can handle genotype data which is selective genotyping.
}
\examples{
# load the example data
load(system.file("extdata", "exampledata.RDATA", package = "QTLEMM"))

# make the seletive genotyping data
ys <- y[y > quantile(y)[4] | y < quantile(y)[2]]
yu <- y[y >= quantile(y)[2] & y <= quantile(y)[4]]
geno.s <- geno[y > quantile(y)[4] | y < quantile(y)[2],]

# run and result
D.matrix <- D.make(3, type = "RI", aa = c(1, 3, 2, 3), dd = c(1, 2, 1, 3), ad = c(1, 2, 2, 3))
result <- EM.MIM2(QTL, marker, geno.s, D.matrix, y = ys, yu = yu, sele.g = "f")
result$E.vector
}
\references{
KAO, C.-H. and Z.-B. ZENG 1997 General formulas for obtaining the maximum
likelihood estimates and the asymptotic variance-covariance matrix in QTL
mapping when using the EM algorithm. Biometrics 53, 653-665. <doi: 10.2307/2533965.>

KAO, C.-H., Z.-B. ZENG and R. D. TEASDALE 1999 Multiple interval mapping
for Quantitative Trait Loci. Genetics 152: 1203-1216. <doi: 10.1093/genetics/152.3.1203>

H.-I LEE, H.-A. HO and C.-H. KAO 2014 A new simple method for improving
QTL mapping under selective genotyping. Genetics 198: 1685-1698. <doi: 10.1534/genetics.114.168385.>
}
\seealso{
\code{\link[QTLEMM]{D.make}}
\code{\link[QTLEMM]{Q.make}}
\code{\link[QTLEMM]{EM.MIM}}
}
