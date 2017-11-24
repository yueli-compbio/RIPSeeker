\name{computeLogOdd}
\alias{computeLogOdd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compute the log odd ratio of RIP over background.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The RIPScore is computed as the log odd ratio of the posterior for the RIP state (\eqn{z_i = 2}) over the posterior for the background state (\eqn{z_i = 1})
}
\usage{
computeLogOdd(nbhGR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nbhGR}{
%%     ~~Describe \code{nbhGR} here~~
	GRanges of bins with the value slot saved for the posterior probabilities for the background and RIP state.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
To assess the statistical significance of the RIP predictions, we assign each bin a RIPScore defined as the log odd ratio of the posterior for the RIP state (\eqn{z_i = 2}) over the posterior for the background state (\eqn{z_i = 1}). When control is available, the RIPScore is updated as the difference between the RIPScores evaluated separately for RIP and control libraries. The scoring system captures the model confidence for the RIP state of each bin in the RIP library penalized by the false confidence for the RIP state of the same bin in the control library. In addition, RIPScore obviates scaling of read counts. Since sequencing depth usually differs between RIP and control libraries, scaling is necessary if the statistical score were derived from the read count differences. On the other hand, simplistic linear scaling may distort the data.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A vector of log odd scores for each bin in \code{nbhGR}.
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
\code{\link{seekRIP}, \link{scoreMergedBins}, \link{logScoreWithoutControl}, \link{logScoreWithControl}}
}
\examples{
# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

alignGal <- getAlignGal(bamFiles[1], reverseComplement=TRUE, genomeBuild="mm9")

alignGR <- as(alignGal, "GRanges")

alignGRList <- GRangesList(as.list(split(alignGR, seqnames(alignGR))))

################ run main function for HMM inference on a single chromosome ################
nbhGR <- mainSeekSingleChrom(alignGR=alignGRList$chrX, K = 2, binSize=1e5)

ripscore <- computeLogOdd(nbhGR)


}