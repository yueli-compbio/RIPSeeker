\name{plotCoverage}
\alias{plotCoverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plot read coverage for a GRanges object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
An internal function used by \code{\link{plotStrandedCoverage}} to plot read counts within each fixed bin across the entire chromosome. 
}
\usage{
plotCoverage(x, plotLegend = FALSE, legend.cex = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
	GRanges object with values slot saved for read counts within the corresponding ranges.
}
  \item{plotLegend}{
%%     ~~Describe \code{plotLegend} here~~
	Binary indcator. If TRUE, legend will be plotted on the top left the plot. Legend is expected to be the chromsome name and length, which must be available in the GRange object argument.  
}
  \item{legend.cex}{
%%     ~~Describe \code{legend.cex} here~~
	Font size of the legend.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
	Extra arguments passed to either the \code{\link[graphics]{plot}} or the \code{\link[graphics]{legend}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
	The read counts is plotted in blue bars as positive integer across the x-axis as the sorted positions across the chromosome. The plot can be used to examine the overall alignment properties for each chromosome.
}
% \value{
% %%  ~Describe the value returned
% %%  If it is a LIST, use
% %%  \item{comp1 }{Description of 'comp1'}
% %%  \item{comp2 }{Description of 'comp2'}
% %% ...
% The coverage plot drawn to the console.
% }
\references{
%% ~put references to the literature/web site here ~
P. Aboyoun, H. Pages and M. Lawrence (). GenomicRanges: Representation and manipulation of genomic intervals. R package version 1.8.9.
}
\author{
%%  ~~who you are~~
Yue Li
}
\note{
%%  ~~further notes~~
Users are not recommanded run this function directly but rather via a much more user friendly function \code{\link{plotStrandedCoverage}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{plotStrandedCoverage}, \link[graphics]{plot}, \link[graphics]{legend}}
}
\examples{
# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

alignGal <- getAlignGal(bamFiles[1], reverseComplement=TRUE, genomeBuild="mm9")

alignGR <- as(alignGal, "GRanges")

alignGRList <- GRangesList(as.list(split(alignGR, seqnames(alignGR))))

binSize <- 1000

binGR <- binCount(alignGRList$chrX, binSize)

plotCoverage(binGR, plotLegend=TRUE)


}