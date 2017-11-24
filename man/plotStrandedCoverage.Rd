\name{plotStrandedCoverage}
\alias{plotStrandedCoverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plot strand-specific read coverage for a GRanges object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Plot read counts within fixed bin across the entire chromosome.
}
\usage{
plotStrandedCoverage(gr, binSize = 1000, plotLegend = FALSE, ylim, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gr}{
%%     ~~Describe \code{gr} here~~
	\code{\link[GenomicRanges]{GRanges}} object containing the alignments.
}
  \item{binSize}{
%%     ~~Describe \code{binSize} here~~
	Integer indicate the size of the bin used to compute and plot the read counts.
}
  \item{plotLegend}{
%%     ~~Describe \code{plotLegend} here~~
	Binary indcator. If TRUE, legend will be plotted on the top left the plot. Legend is expected to be the chromsome name and length, which must be available in the GRange object argument.
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
	A two element scale on the y-axis, indicating the maximum read counts on the + and - strand to be plotted (e.g., ylim=c(-200, 200)).
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
	Extra arguments passed to \code{\link{plotCoverage}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
	Read count on + and - strand are displayed as red and blue bars on the positive and negative y-axis, respectively. The x-axis indicates the positions across the chromosmoe. The plot can be used to examine for each chromosome the overall alignment properties such as strand specificity (expected in non-strand-specific sequencing) and aggregation of reads.
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
P. Aboyoun, H. Pages and M. Lawrence (). GenomicRanges: Representation and manipulation of genomic intervals. R package version 1.8.9.
}
\author{
%%  ~~who you are~~
Yue Li
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{plotCoverage}, \link[graphics]{plot}, \link[graphics]{legend}}
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

plotStrandedCoverage(gr=alignGRList$chrX, binSize=binSize,

		xlab="", ylab="", plotLegend=TRUE, box.lty=0, legend.cex=2 )

}