\name{addDummyProb}
\alias{addDummyProb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Create a dummy GRanges object as a placeholder in case nbh_em fails (Internal function)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used to generate a place holder in cases the EM fails to converge on a chromosome due to too few number of reads mapped to that chromosome. This is an internal function not expected to be directly called by the user.
}
\usage{
addDummyProb(alignGR, K = 2, randomProb = FALSE, runViterbi = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignGR}{
%%     ~~Describe \code{alignGR} here~~
	GRanges object derived from RIP-seq alignment inputs. 
}
  \item{K}{
%%     ~~Describe \code{K} here~~
	Number of hidden states (Default: 2).
}
  \item{randomProb}{
%%     ~~Describe \code{randomProb} here~~
	A binary value to indicate whether to use random probability as a place holder to present posterior probabilities. If set FALSE, (by default), equal probability is used for all states.
}
  \item{runViterbi}{
%%     ~~Describe \code{runViterbi} here~~
	A binary value to indicate whether to generate place holder for the Viterbi state sequence (Default: FALSE).
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
	Additional arguments expected to contain the bin size used for computing the bin counts in \code{\link{binCount}} function, and any other extra arguments are ignored.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
	A priviate function to fall back in case HMM fails to converge mostly due to too many zero counts in the input vector. When that occurs, a GRanges place holder object needs to be returned to keep consistent with the remaining GRanges for each chromosome. Thus, all information slot will be generated as place holder to properly create the GRangesList for the predictions on all chromosomes (each as an GRanges item in the list).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{GRanges}{A GRanges object containing the read count (in the defined bin size), alpha, beta, TRANS dummy values for the HMM}
}
% \references{
% %% ~put references to the literature/web site here ~
% }
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
	\code{\link{mainSeekSingleChrom}}
}
\examples{
# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

alignGal <- getAlignGal(bamFiles[1], reverseComplement=TRUE, genomeBuild="mm9")

alignGR <- as(alignGal, "GRanges")

alignGRList <- GRangesList(as.list(split(alignGR, seqnames(alignGR))))

x <- addDummyProb(alignGRList$chrX, binSize=10000)

x

}