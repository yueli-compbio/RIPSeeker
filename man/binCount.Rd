\name{binCount}
\alias{binCount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Count reads in nonoverlapping bins across a chromosome
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Stratify chromosome into nonoverlapping bins of the same size and count the number of reads that fall within each bin.
}
\usage{
binCount(alignGR, binSize, returnBinCountOnly = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignGR}{
%%     ~~Describe \code{alignGR} here~~
	GRanges object containing the alignments for a single chromosome.
}
  \item{binSize}{
%%     ~~Describe \code{binSize} here~~
	An integer for the bin size.
}
  \item{returnBinCountOnly}{
%%     ~~Describe \code{returnBinCountOnly} here~~
	Binary indicator. If TRUE, only the integer read count is returned; if FALSE, GRanges of bins with value slot saved as the corresponding read counts is returned.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The funciton is designed to operate to a single chromosome to facilitate parallel computing on multiple chromosomes independently. The function is used in \code{\link{evalBinSize}} to select the optimal bin size based on the read counts and in \code{\link{mainSeekSingleChrom}} to provide the read count as input for the HMM.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
If \code{returnBinCountOnly} is TRUE, then the integer read count is returned; if \code{returnBinCountOnly} is FALSE, then the GRanges of bins with value slot saved for the corresponding read counts is returned.}

\references{
%% ~put references to the literature/web site here ~
	P. Aboyoun, H. Pages and M. Lawrence. GenomicRanges: Representation
  and manipulation of genomic intervals. R package version 1.8.9.
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
	\code{\link{selectBinSize}, \link{evalBinSize}}
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


}