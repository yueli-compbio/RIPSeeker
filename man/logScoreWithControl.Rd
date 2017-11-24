\name{logScoreWithControl}
\alias{logScoreWithControl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compute RIPScore based on RIP and control posteriors and test for significance
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Compute the RIPScore using both RIP and control posteriors for each bins, merge and summarize the scores for the merged bins, and finally compute the p-value and adjusted p-value for the summary RIPScore.
}
\usage{
logScoreWithControl(nbhGRRIP, nbhGRCTL, padjMethod = "BH", getControlStats = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nbhGRRIP}{
%%     ~~Describe \code{nbhGRRIP} here~~
	\code{\link[GenomicRanges]{GRanges}} object for the RIP library created from \code{\link{mainSeek}} containing the posterior probabilities of the hidden states for each observed read count.
}
  \item{nbhGRCTL}{
%%     ~~Describe \code{nbhGRCTL} here~~
	An optional arugment as a \code{\link[GenomicRanges]{GRanges}} object for the control library created from \code{\link{mainSeek}} containing the posterior probabilities of the hidden states for each observed read count.
}
  \item{padjMethod}{
%%     ~~Describe \code{padjMethod} here~~
	Method used to adjust multiple testing performed in \code{\link[stats]{p.adjust}} (Default: "BH").
}
  \item{getControlStats}{
%%     ~~Describe \code{getControlStats} here~~
	Binary indicator to whether return statistics including computed for the control library alone. If TRUE,  then all control specific score columns will be reported with a prefix "CTL".
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The RIPScore is compupted in \code{\link{computeLogOdd}} as the log odd ratio of the posterior for the RIP state (\eqn{z_{i} = 2}) over the posterior for the background state (\eqn{z_{i} = 1}) in RIP library subtracted by the log odd ratio computed from the control library. The adjacent bins with hidden states predicted by \code{\link{nbh_vit}} as the enriched state (corresponding to the NB with larger mean) are merged. The RIPSscores are averaged over the merged bins. To assess the statistical significance of the RIPScore for each region, we assume that the RIPScore follows a \emph{Gaussian} (Normal) distribution with  mean and standard deviation  estimated using the RIPScores over all of the bins. The rationale is based on the assumption that most of the RIPScores correspond to the background state and together contribute to a stable estimate of the test statistics (TS) and p-value computed using the R built-in function \code{\link[stats]{pnorm}}. The p-value is adjusted by \code{\link[stats]{p.adjust}} with BH method by default. The same procedure is applied optionally to the control library.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\code{\link[GenomicRanges]{GRanges}} of merged bins with values slot saved for RIPScore (lodOdd), p-value (pval), adjusted p-value (pvalAdj) for RIP and optionally for control.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
%%  ~~who you are~~
Yue Li
}
\note{
%%  ~~further notes~~
Internal function used by \code{\link{seekRIP}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{logScoreWithoutControl}, \link{seekRIP}, \link{computeLogOdd}, \link{scoreMergedBins}}
}
\examples{
if(interactive()) { # check the example in seekRIP
# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

# Parameters setting
binSize <- 1e5							# use a large fixed bin size for demo only
multicore <- FALSE					# use multicore
strandType <- "-"						# set strand type to minus strand

################ run main function for HMM inference on all chromosomes ################
mainSeekOutputRIP <- 
    mainSeek(bamFiles=grep(pattern="SRR039214", 
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


################ Compute log score and test for significance WITH control ################
ripGR.wicontrol <- logScoreWithControl(mainSeekOutputRIP$nbhGRList$chrX, mainSeekOutputCTL$nbhGRList$chrX)

ripGR.wicontrol
}

}