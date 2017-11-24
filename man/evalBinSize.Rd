\name{evalBinSize}
\alias{evalBinSize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Evaluate bin size using Shimazaki cost function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given a bin size and a GRanges alignment object, the funciton computes the bin count and returns the cost of the bin size based on Shimazaki cost function.
}
\usage{
evalBinSize(binSize, alignGR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{binSize}{
%%     ~~Describe \code{binSize} here~~
	An integer that indicates the bin size applied to the binning of the chromosome.
}
  \item{alignGR}{
%%     ~~Describe \code{alignGR} here~~
	GRanges object of alignments to a single chromosome.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
	The function implements the algorithm developed by Shimazaki and Shinomoto (2007), which is based on the goodness of the fit of the time histogram to estimate the rate of neural response of an animal to certain stimuli in a spike-in experiment. The algorithm involves four simple steps:
	
	1. Divide chromosome sequence into \eqn{N} bins of width \eqn{b}.
	
	2. Count number of read counts \eqn{x_i} that enter the \eqn{i}'th bin.
	
	3. Compute: \eqn{\bar{x} = \frac{1}{N}\sum_{i=1}^{N}x_{i}} and \eqn{v = \frac{1}{N}\sum_{i=1}^{N}(x_{i} - \bar{x})^{2}}.
	
	4. Compute: \eqn{C(b) = \frac{2\bar{x}-v}{b^{2}}}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{cost}{A scalar value for the cost of the bin size.}
}
\references{
%% ~put references to the literature/web site here ~
Hideaki Shimazaki and Shigeru Shinomoto. A method for selecting the bin size of a time histogram. Neural computation, 19(6):1503-1527, June 2007.
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
	\code{\link{selectBinSize}, \link{binCount}}
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

costs <- evalBinSize(binSize, alignGRList$chrX)

}