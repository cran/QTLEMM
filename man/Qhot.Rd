% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Qhot.R
\name{Qhot}
\alias{Qhot}
\title{QTL Hotspot}
\usage{
Qhot(DataQTL, DataCrop, ScanStep = 1, NH = 100, NP = 1000, save.pdf = TRUE)
}
\arguments{
\item{DataQTL}{data.frame. A data frame with 5 columns for QTL information.
The columns represent the serial number of QTLs, the trait names, the
chromosome numbers, the left flanking marker positions(in cM) of QTLs, and
the right flanking marker positions(in cM) of QTLs.}

\item{DataCrop}{data.frame. A data frame with 3 columns for chromosome
information consists of the chromosome names, the center positions(in cM)
and the lengths of chromosomes.}

\item{ScanStep}{numeric. A value for the length(cM) of every bin.}

\item{NH}{integer. A value for the number of spurious hotspots in the
proposed method.}

\item{NP}{integer. A value for permutation times to calculate the
threshold.}

\item{save.pdf}{logical. When set to TRUE, the PDF file of plots will be
saved in the working directory instead of being displayed in the console.}
}
\value{
\item{EQF}{The expected QTL frequency(EQF) in every bin per chromosome.}
\item{P.threshold}{The EQF thresholds for proposed method.}
\item{Q.threshold}{The EQF thresholds for the Q method.}
\item{nHot}{The numbers of detected hotspots per chromosome for the proposed
method and Q method.}

Graphical outputs for visualizing the summarized results includes the
expected QTL frequency of scan steps, and the composition of QTLs for
different traits in the detected hotspots.
}
\description{
This function generates both numerical and graphical summaries of the QTL
hotspot detection in the genomes, including information about the flanking
markers of QTLs.
}
\note{
This program may generate a large amount of graphical output. To manage this,
it's recommended to save the results in a PDF file using the "save.pdf"
argument.
}
\examples{
# load the example data
load(system.file("extdata", "QHOTexample.RDATA", package = "QTLEMM"))

# run and result
result <- Qhot(QTL.example, crop.example, 5, 20, 100, save.pdf = FALSE)
}
\references{
Wu, P.-Y., M.-.H. Yang, and C.-H. KAO 2021 A Statistical Framework for QTL
Hotspot Detection. G3: Genes, Genomes, Genetics: jkab056. <doi: 10.1093/g3journal/jkab056>
}
