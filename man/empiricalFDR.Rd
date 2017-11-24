\name{empiricalFDR}
\alias{empiricalFDR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compute empirical false discovery rate
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
At a p-value, find the number of regions in RIP library (denoted as "trueCount") and the number of regions in control library (denoted as "falseCount"). The empirical false discovery rate (eFDR) is estimated as the ratio of the falseCount over the trueCount.
}
\usage{
empiricalFDR(pval, pvalRIP, pvalCTL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pval}{
%%     ~~Describe \code{pval} here~~
	A scalar p-value.
}
  \item{pvalRIP}{
%%     ~~Describe \code{pvalRIP} here~~
	A column vector of p-values for the peaks identifed from RIP v.s. control comparison.
}
  \item{pvalCTL}{
%%     ~~Describe \code{pvalCTL} here~~
	A column vector of p-values for the peaks identifed from control v.s. RIP comparison.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Only when the control is available, is an empirical false discovery rate (eFDR) estimated based on the idea of "sample swap" inspired by MACS (a ChIP-seq algorithm from Zhange \emph{el al.} (2008). At each p-value, RIPSeeker finds the number of significnat RIP-regions over control (CTL) based on \code{pvalRIP} and the number of significant control regions over RIP based on \code{pvalCTL}. The eFDR is defined as the ratio of the number of "RIP" (false positive) regions identified from CTL-RIP comparison over the number of RIP regions from the RIP-CTL comparison. The maximum value for eFDR is 1 and minimum value for eFDR is max(p-value, 0). The former takes care of the case where the numerator is bigger than the denominator, and the latter for zero numerator.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A scalar probabibility value that represents the eFDR.
}
\references{
%% ~put references to the literature/web site here ~
Yong Zhang, Tao Liu, Clifford A Meyer, J\'er\^ome Eeckhoute, David S Johnson, Bradley E Bernstein, Chad Nusbaum, Richard M Myers, Myles Brown, Wei Li, and X Shirley Liu. Model-based analysis of ChIP-Seq (MACS). Genome Biology, 9(9):R137, 2008.
}
\author{
%%  ~~who you are~~
Yue Li
}
\note{
%%  ~~further notes~~
This is an internal funciton used in \code{\link{seekRIP}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{logScoreWithControl}, \link{seekRIP}, \link{computeLogOdd}, \link{scoreMergedBins}}
}
\examples{
pvalRIP <- runif(100)

pvalCTL <- runif(100)

eFDR <- empiricalFDR(pvalRIP[1], pvalRIP, pvalCTL)

pvalRIP[1]

eFDR

# more significant pval
pvalRIP[1] <- 1e-4

eFDR <- empiricalFDR(pvalRIP[1], pvalRIP, pvalCTL)

pvalRIP[1]

eFDR
}