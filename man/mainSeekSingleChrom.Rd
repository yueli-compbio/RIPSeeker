\name{mainSeekSingleChrom}
\alias{mainSeekSingleChrom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Automatic bin size selection, bin count, and HMM parameters optimization on read count vector from a single chromosome (Internal function)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This an internal function used by \code{\link{mainSeek}} to accomplish three major tasks \emph{on a single chromosome}: automatically select bin size, compute read counts within the bins, and obtain optimal HMM paramters.
}
\usage{
mainSeekSingleChrom(alignGR, K = 2, binSize = NULL, minReadCount = 10, 
	backupNumBins = 10, minBinSize = 200, maxBinSize = 1200, 
	increment = 5, pathToSavePlotsOfBinSizesVersusCosts, 
	verbose = TRUE, allowSecondAttempt = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignGR}{
%%     ~~Describe \code{alignGR} here~~
	GRanges containing the alignments on a single chromosome .
}
  \item{K}{
%%     ~~Describe \code{K} here~~
	Number of hidden states (Default: 2). By default, state 1 specifies the background and state 2 the RIP regions. The two states are recognized by the means for the two distributions (See \code{\link{nbh_em}}).
}
  \item{binSize}{
%%     ~~Describe \code{binSize} here~~
	Size to use for binning the read counts across each chromosome. If NULL, optimal bin size within a range (default: minBinSize=200, maxBinSize=1200) will be automatically selected (See \code{\link{selectBinSize}}).
}
  \item{minReadCount}{
%%     ~~Describe \code{minReadCount} here~~
	Minimum aligned read counts needed for HMM to converge (Default: 10). Note that HMM may not converge some times when majority of the read counts are zero even if some read count > 10. When that happens, a back-up function \code{\link{addDummyProb}} comes in to create a placeholder for the corresponding chromosome in GRangeList to maintain the data structure to preserve all information (successfully) obtained from other chromosomes.
}
  \item{backupNumBins}{
%%     ~~Describe \code{backupNumBins} here~~
	If read count is less than \code{minReadCount}, then use \code{backupNumBins} (Default: 10) to bin the chromosome.
}
  \item{minBinSize}{
%%     ~~Describe \code{minBinSize} here~~
	Minimum bin size to start with the bin selection (See \code{\link{selectBinSize}}). Default to 200, common minimum band size selected in RIP or RNA-seq library construction.
}
  \item{maxBinSize}{
%%     ~~Describe \code{maxBinSize} here~~
	Maximum bin size to stop with the bin selection (See \code{\link{selectBinSize}}). Default: 1200.
}
  \item{increment}{
%%     ~~Describe \code{increment} here~~
	Step-wise increment in bin size selection (See \code{\link{selectBinSize}}). Default: 5.
}
  \item{pathToSavePlotsOfBinSizesVersusCosts}{
%%     ~~Describe \code{pathToSavePlotsOfBinSizesVersusCosts} here~~
	Directory used to save the diagnostic plots for bin size selection.
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
	Binary indicator for disable (FALSE) or enable (TRUE) HMM training message from function \code{\link{nbh}} to output to the console.
}
  \item{allowSecondAttempt}{
%%     ~~Describe \code{allowSecondAttempt} here~~
	In case HMM fails to converge due to malformed paramters in EM iteraction, re-iterating the HMM process each time with a different suboptimal bin size in attempt to succeed in some trial. If all yeild nothing, fall back up to \code{\link{addDummyProb}} to return the place holder for the chromosome.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
	Argumnets passed to \code{\link{nbh}}.
}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~	
% }
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
	\item{nbhGR}{GRanges object containing the optimized HMM parameters (and the Viterbi hidden state sequence) accompanied with the read count vector following the (automatic) binning scheme.}
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
Unless a highly customized workflow is needed, \code{\link{ripSeek}} is the high-level front-end main function that should be used in most cases.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
	\code{\link{ripSeek}, \link{mainSeek}, \link{nbh_em}}
}
\examples{
# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

# Parameters setting
binSize <- 1e5							  # use a large fixed bin size for demo only
minBinSize <- NULL						# min bin size in automatic bin size selection
maxBinSize <- NULL						# max bin size in automatic bin size selection
multicore <- FALSE						# use multicore
strandType <- "-"							# set strand type to minus strand

# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

alignGal <- getAlignGal(bamFiles[1], reverseComplement=TRUE, genomeBuild="mm9")

alignGR <- as(alignGal, "GRanges")

alignGRList <- GRangesList(as.list(split(alignGR, seqnames(alignGR))))

################ run main function for HMM inference on a single chromosome ################
nbhGR <- mainSeekSingleChrom(alignGR=alignGRList$chrX, K = 2, binSize=binSize, 
			
			minBinSize = minBinSize, maxBinSize = maxBinSize)
			
nbhGR			

}