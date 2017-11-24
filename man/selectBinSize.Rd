\name{selectBinSize}
\alias{selectBinSize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Select optimal bin size based on Shimazaki formula
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function iteratively estimates the cost of increasing bin size within a defined range and finally selects the bin size with minimum cost.
}
\usage{
selectBinSize(alignGR, minBinSize, maxBinSize = 1000, 
	increment = 5, getFullResults = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignGR}{
%%     ~~Describe \code{alignGR} here~~
	GRanges object of alignments on a single chromosome
}
  \item{minBinSize}{
%%     ~~Describe \code{minBinSize} here~~
	Minimum bin size to start with (Default: 5 * read length)
}
  \item{maxBinSize}{
%%     ~~Describe \code{maxBinSize} here~~
	Maximum bin size to end with (Default: 1000).
}
  \item{increment}{
%%     ~~Describe \code{increment} here~~
	Number of bases to increment the bin size in searching for the optimal bin size within the defined range (Default: 5).	
}
  \item{getFullResults}{
%%     ~~Describe \code{getFullResults} here~~
	Binary indicator. If TRUE, the optimal bin size (with the minimum cost), minimum cost, and all of the bin sizes considered and their costs are returned. If FALSE, only the optimal bin size is returned.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Based on the preprocessed alignments for a chromosome, RIPSeeker divides the chromosome into bins of equal size \eqn{b} and compute the number of reads that \eqn{b} needs to be determined either empirically (e.g., based on the gel-selected length of the RNA fragment) or computationally. If the bin size is too small, the read counts fluctuates greatly, making it difficult to discern the underlying read count distribution. Additionally, input size to HMM increases as bin size decreases. A very small bin size results in a very long Markov chain of read counts to model, making the computation inefficient. On the other hand, if a bin size is too large, resolution becomes poor. Consequently, one cannot detect the local RIP region with subtle but intrinsic difference from the background, and the RIP regions tend to be too wide for designing specific primer for validation.

Intuitively, selecting an appropriate bin size for each chromosome is metaphorically equivalent to choosing an optimal intervals for building a histogram (Song, 2011). Here we implement the algorithm developed by Shimazaki and Shinomoto (2007), which is based on the goodness of the fit of the time histogram to estimate the rate of neural response of an animal to certain stimuli in a spike-in experiment. This approach has been successfully applied in a recently developed ChIP-seq program (Song and Smith, 2011). Algorithm 1 describes the pseudocode adapted from Shimazaki and Shinomoto (2007) that iteratively estimates the cost \eqn{C} of increasing bin size \eqn{b} within a defined range is outlined as follows.

	For b = minBinSize to maxBinSize; do

		1. Divide chromosome sequence into \eqn{N} bins of width \eqn{b}.
	
		2. Count number of read counts \eqn{x_i} that enter the \eqn{i}'th bin.
		
		3. Compute: \eqn{\bar{x} = \frac{1}{N}\sum_{i=1}^{N}x_{i}} and \eqn{v = \frac{1}{N}\sum_{i=1}^{N}(x_{i} - \bar{x})^{2}}.
		
		4. Compute: \eqn{C(b) = \frac{2\bar{x}-v}{b^{2}}}
	
	End For
		
		Choose \eqn{b} that minimize \eqn{C(b)}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
When \code{getFullResults} is TRUE, return a list containing:

	\item{bestBinSize}{the optimal bin size (with the minimum cost)}
	\item{minCosts}{cost of the optimal bin size}
	\item{binSizes}{all of the bin sizes considered}
	\item{costs}{all of the costs}
	
When \code{getFullResults} is FALSE, only the optimal bin size (\code{bestBinSize}) is returned.
}
\references{
%% ~put references to the literature/web site here ~
	Hideaki Shimazaki and Shigeru Shinomoto. A method for selecting the bin size of a time histogram. Neural computation, 19(6):1503-1527, June 2007.
	
	Qiang Song and Andrew D. Smith. Identifying dispersed epigenomic domains from ChIP-Seq data. Bioinformatics (Oxford, England), 27(6):870-871, March 2011.
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
	\code{\link{evalBinSize}, \link{binCount}}
}
\examples{
# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

alignGal <- getAlignGal(bamFiles[1], reverseComplement=TRUE, genomeBuild="mm9")

alignGR <- as(alignGal, "GRanges")

alignGRList <- GRangesList(as.list(split(alignGR, seqnames(alignGR))))

minBinSize <- 200

maxBinSize <- 1200

gr <- alignGRList$chrX

b <- selectBinSize(gr, minBinSize, maxBinSize, increment=100, getFullResults=TRUE)

plot(b$binSizes, b$costs)

chrname <- as.character(runValue(seqnames(gr)))

chrlen <- seqlengths(gr)[chrname]

legend("topright", box.lty=0,
		
		sprintf("\%s: 1-\%d;\nTotal mapped reads: \%d;\nOptimal bin size = \%d bp",
		
		chrname, chrlen, length(gr), b$bestBinSize))	
}