\name{seekRIP}
\alias{seekRIP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Identify significant peaks
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Based on the posteriors derived from HMM by \code{\link{mainSeek}}, find the significant RIP regions derived from merging the adjacent RIP bins. The significance test makes use of the log odd ratio of the posteriors for RIP over the background states.
}
\usage{
seekRIP(nbhGRRIP, nbhGRCTL = NULL, padjMethod = "BH", 
	logOddCutoff = -Inf, pvalCutoff = 1, pvalAdjCutoff = 1, 
	eFDRCutoff = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nbhGRRIP}{
%%     ~~Describe \code{nbhGRRIP} here~~
	\code{\link[GenomicRanges]{GRanges}} object for the RIP library created from \code{\link{mainSeek}} containing the posterior probabilities of the hidden states for each observed read count.
}
  \item{nbhGRCTL}{
%%     ~~Describe \code{nbhGRCTL} here~~
	An optional argument as a \code{\link[GenomicRanges]{GRanges}} object for the control library created from \code{\link{mainSeek}} containing the posterior probabilities of the hidden states for each observed read count.
}
  \item{padjMethod}{
%%     ~~Describe \code{padjMethod} here~~
	Method used to adjust multiple testing performed in \code{\link[stats]{p.adjust}} (Default: "BH").
}
  \item{logOddCutoff}{
%%     ~~Describe \code{logOddCutoff} here~~
	Threshold for the  log odd ratio of posterior for the RIP over the background states (See \code{\link{seekRIP}}). Only peaks with logOdd score \emph{greater} than the \code{logOddCutoff} will be reported. Default: -Inf (i.e. no cutoff).
}
  \item{pvalCutoff}{
%%     ~~Describe \code{pvalCutoff} here~~
	Threshold for the adjusted p-value for the logOdd score. Only peaks with adjusted p-value \emph{less} than the \code{pvalAdjCutoff} will be reported. Default: 1 (i.e. no cutoff).
}
  \item{pvalAdjCutoff}{
%%     ~~Describe \code{pvalAdjCutoff} here~~
	Threshold for the adjusted p-value for the logOdd score. Only peaks with adjusted p-value \emph{less} than the \code{pvalAdjCutoff} will be reported. Default: 1 (i.e. no cutoff).
}
  \item{eFDRCutoff}{
%%     ~~Describe \code{eFDRCutoff} here~~
	Threshold for the empirical false discovery rate (eFDR). Only peaks with eFDR \emph{less} than the \code{pvalAdjCutoff} will be reported. Default: 1 (i.e. no cutoff).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
	The RIPScore is compupted in \code{\link{computeLogOdd}} as the log odd ratio of the posterior for the RIP state (\eqn{z_{i} = 2}) over the posterior for the background state (\eqn{z_{i} = 1}) in RIP library. When control is available, the RIPScore is updated by the difference between the RIPScores between RIP and control. The adjacent bins with hidden states predicted by \code{\link{nbh_vit}} as the enriched state (corresponding to the NB with larger mean) are merged. The RIPSscores are averaged over the merged bins. To assess the statistical significance of the RIPScore for each region, we assume that the RIPScore follows a \emph{Gaussian} (Normal) distribution with  mean and standard deviation  estimated using the RIPScores over all of the bins. The rationale is based on the assumption that most of the RIPScores correspond to the background state and together contribute to a stable estimate of the test statistics (TS) and p-value computed using the R built-in function \code{\link[stats]{pnorm}}.
	
	The p-value is adjusted by \code{\link[stats]{p.adjust}} with BH method by default. The same procedure is applied optionally to the control library. Only when the control is available, is an empirical false discovery rate (eFDR) estimated based on the idea of "sample swap" inspired by MACS (a ChIP-seq algorithm from Zhange \emph{el al.} (2008). At each p-value, RIPSeeker finds the number of significnat RIP-regions over control (CTL) based on \code{pvalRIP} and the number of significant control regions over RIP based on \code{pvalCTL}. The eFDR is defined as the ratio of the number of "RIP" (false positive) regions identified from CTL-RIP comparison over the number of RIP regions from the RIP-CTL comparison. The maximum value for eFDR is 1 and minimum value for eFDR is max(p-value, 0). The former takes care of the case where the numerator is bigger than the denominator, and the latter for zero numerator.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\code{\link[GenomicRanges]{GRanges}} object containing the merged bins with values slot saved for RIPScore (lodOdd), p-value (pval), adjusted p-value (pvalAdj) for RIP and optionally for control.}
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
Internal function used by \code{\link{ripSeek}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{logScoreWithControl}, \link{logScoreWithoutControl}, \link{empiricalFDR}, \link{computeLogOdd}, \link{scoreMergedBins}, \link{ripSeek}}
}
\examples{
if(interactive()) {
# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

# Parameters setting
binSize <- 1e5        # use a large fixed bin size for demo only
multicore <- FALSE    # use multicore
strandType <- "-"			# set strand type to minus strand

################ run main function for HMM inference on all chromosomes ################
mainSeekOutputRIP <- mainSeek(bamFiles=grep(pattern="SRR039214", 
    bamFiles, value=TRUE, invert=TRUE),
		binSize=binSize, strandType=strandType, 		
		reverseComplement=TRUE, genomeBuild="mm9",
		uniqueHit = TRUE, assignMultihits = TRUE, 
		rerunWithDisambiguatedMultihits = FALSE,				
		multicore=multicore, silentMain=FALSE, verbose=TRUE)
		
mainSeekOutputCTL <- mainSeek(bamFiles=grep(pattern="SRR039214", 
    bamFiles, value=TRUE, invert=FALSE),
		binSize=binSize, strandType=strandType, 		
		reverseComplement=TRUE, genomeBuild="mm9",
		uniqueHit = TRUE, assignMultihits = TRUE, 
		rerunWithDisambiguatedMultihits = FALSE,				
		multicore=multicore, silentMain=FALSE, verbose=TRUE)

# with control
ripGR.wicontrol <- seekRIP(mainSeekOutputRIP$nbhGRList$chrX, mainSeekOutputCTL$nbhGRList)

# without control
ripGR.wocontrol <- seekRIP(mainSeekOutputRIP$nbhGRList$chrX)
}
}