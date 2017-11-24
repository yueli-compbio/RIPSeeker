\name{viewRIP}
\alias{viewRIP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Visualize peaks from UCSC genome browser.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Upload alignments, peaks, statistical scores to UCSC genome browser for comparative visualization of the results and data available in the UCSC database.
}
\usage{
viewRIP(seekedRIP, alignGR, alignGRCTL, 
	binGR = seekedRIP, scoreType = "eFDR", 
	cutoffLine = 0.001, displayALLChr = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seekedRIP}{
%%     ~~Describe \code{seekedRIP} here~~
	\code{\link[GenomicRanges]{GRangesList}} obtained from \code{\link{ripSeek}}. Each list item represents the RIP peaks on a chromosome accompanied with statistical scores including (read) count, logOddScore, pval, pvalAdj, eFDR for the RIP and control (if available). Please refer to \code{\link{seekRIP}} for more details.
}
  \item{alignGR}{
%%     ~~Describe \code{alignGR} here~~
	\code{\link[GenomicRanges]{GRanges}} of read alignments for the RIP.
}
  \item{alignGRCTL}{
%%     ~~Describe \code{alignGRCTL} here~~
	\code{\link[GenomicRanges]{GRanges}} of read alignments for the control.
}
  \item{binGR}{
%%     ~~Describe \code{binGR} here~~
	\code{\link[GenomicRanges]{GRanges}} containing read count column corresponding to the peaks. By default, \code{alignGR} is used as \code{binGR} to display the read count in RIP condition. 
}
  \item{scoreType}{
%%     ~~Describe \code{scoreType} here~~
	Type of statistical score to display as another track in the browser (Default: eFDR). eFDR/pval/pvalAdj is displayed at -log10 scale.
}
  \item{cutoffLine}{
%%     ~~Describe \code{cutoffLine} here~~
	Draw a cutoffline in the browser to indicate the significance level above which the peaks are considered significant.
}
  \item{displayALLChr}{
%%     ~~Describe \code{displayALLChr} here~~
	Binary indicator when TRUE upload and display the information for only one chromosome rather than upload all chromosomes (Default: TRUE).
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
	Extra arguments are ignored.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
	The function is a wrapper function of \code{\link[rtracklayer]{browserSession}}, \code{\link[rtracklayer]{track}}, and \code{\link[rtracklayer]{browserView}}.
}
% \value{
% %%  ~Describe the value returned
% %%  If it is a LIST, use
% %%  \item{comp1 }{Description of 'comp1'}
% %%  \item{comp2 }{Description of 'comp2'}
% %% ...
% }
\references{
%% ~put references to the literature/web site here ~
Michael Lawrence, Vince Carey and Robert Gentleman (). rtracklayer: R interface to genome browsers and their annotation tracks. R package version 1.16.3.
}
\author{
%%  ~~who you are~~
Yue Li
}
\note{
%%  ~~further notes~~
If input contain multiple chromosomes, then multiple browser window will be open to display each chromosome. A more user-friendly way is to upload all of the information to UCSC and open a single browser for visualization, which may become one of the new features in future release.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{ripSeek}, \link[rtracklayer]{browserSession}, \link[rtracklayer]{track}, \link[rtracklayer]{browserView}}
}
\examples{
if(interactive()) { # need internet connection
# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

cNAME <- "SRR039214" 						# specify control name


# Parameters setting
binSize <- NULL							# automatically determine bin size
minBinSize <- 10000						# min bin size in automatic bin size selection
maxBinSize <- 12000						# max bin size in automatic bin size selection
multicore <- TRUE						# use multicore
strandType <- "-"						# set strand type to minus strand


################ run main function ripSeek to predict RIP ################
seekOut <- ripSeek(bamPath=bamFiles, cNAME=cNAME, 
		binSize=binSize, minBinSize = minBinSize, 
		maxBinSize = maxBinSize, strandType=strandType, 
		silentMain=TRUE, verbose=FALSE, 
		reverseComplement=TRUE, genomeBuild="mm9",
		uniqueHit = TRUE, assignMultihits = TRUE, 
		rerunWithDisambiguatedMultihits = TRUE, multicore=multicore)


################ visualization ################

viewRIP(seekOut$RIPGRList$chrX, seekOut$mainSeekOutputRIP$alignGalFiltered, 
	seekOut$mainSeekOutputCTL$alignGalFiltered, scoreType="eFDR")
}
}