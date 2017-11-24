\name{RIPSeeker-package}
\alias{RIPSeeker-package}
\alias{RIPSeeker}
\docType{package}
\title{
RIPSeeker: a statistical package for identifying protein-associated transcripts from RIP-seq experiments
}
\description{
RIPSeeker infers and discriminates RIP peaks from RIP-seq alignments using two-state HMM with negative binomial emission probability. While RIPSeeker is specifically tailored for RIP-seq data analysis, it also provides a suite of bioinformatics tools integrated within this self-contained software package comprehensively addressing issues ranging from post-alignments processing to visualization and annotation. In addition, a rule-based approach is provided as an additional function named \code{\link{rulebaseRIPSeek}} for user to obtain RPKM/FPKM (and fold-change) for the gene/transcripts expressions in RIP (and control) based on automatically retrieved online Ensembl annotation given single or paired-end alignments.
}
\details{
\tabular{ll}{
Package: \tab RIPSeeker\cr
Type: \tab Package\cr
Version: \tab 1.4.0\cr
Date: \tab 2012-11-06\cr
License: \tab GPL-2\cr
}
The front-end main function \code{\link{ripSeek}} suffices for most applications. The function takes as the only required argument the path to alignment files (BAM/BED/SAM) and outputs predicted RIP regions. Optionally, user may indicate via 'cNAME' which file(s) in the first file argument list is/are control to enable empirical false discover rate (eFDR) computation. If the arguments 'biomaRt_dataset' and/or 'goAnno' are set, ripSeek will return the annotated RIP predictions and the enriched GO terms corresponding to the genomic context of the RIP predictions. User can also specify the thresholds for statistical significance scores via \code{logOddCutoff}, \code{pvalCutoff}, \code{pvalAdjCutoff}, \code{eFDRCutoff}.
}
\author{
Yue Li <yueli@cs.toronto.edu>
}
\references{
Li, Y., Zhao, D. Y., Greenblatt, J. F., & Zhang, Z. (2013). RIPSeeker: a statistical package for identifying protein-associated transcripts from RIP-seq experiments. Nucleic Acids Research. doi:10.1093/nar/gkt142

Zhao, J., Ohsumi, T. K., Kung, J. T., Ogawa, Y., Grau, D. J., Sarma, K., Song, J. J., et al. (2010). Genome-wide Identification of Polycomb-Associated RNAs by RIP-seq. Molecular Cell, 40(6), 939D953. doi:10.1016/j.molcel.2010.12.011
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%%~~ the R documentation directory ~~
\keyword{ 
	package 
	high-throughput sequencing analysis
	RIP-seq
	Hidden Markov model
}
\seealso{
\code{\link{ripSeek}, \link{rulebaseRIPSeek}}
}
\examples{

library(RIPSeeker)

ls("package:RIPSeeker")

}
