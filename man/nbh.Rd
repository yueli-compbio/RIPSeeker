\name{nbh}
\alias{nbh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generic function of negative binomial HMM
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generic function for \code{\link{nbh.GRanges}} and \code{\link{nbh.integer}}
}
\usage{
nbh(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{object} here~~
	Object of class Integer or GRanges.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
	Extra arguments passed to either \code{\link{nbh.GRanges}} or \code{\link{nbh.integer}}.
}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
% \value{
% %%  ~Describe the value returned
% %%  If it is a LIST, use
% %%  \item{comp1 }{Description of 'comp1'}
% %%  \item{comp2 }{Description of 'comp2'}
% %% ...
% }
\references{
%% ~put references to the literature/web site here ~
Rabiner, L. R. (1989). A tutorial on hidden Markov models and selected applications in speech recognition (Vol. 77, pp. 257-286). Presented at the Proceedings of the IEEE. doi:10.1109/5.18626

Bishop, Christopher. Pattern recognition and machine learning. Number 605-631 in Information Science and Statisitcs. Springer Science, 2006.

Capp\'e, O. (2001). H2M : A set of MATLAB/OCTAVE functions for the EM estimation of mixtures and hidden Markov models. (\url{http://perso.telecom-paristech.fr/cappe/h2m/})
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
\code{\link{mainSeekSingleChrom}, \link{nbh.integer}, \link{nbh.GRanges}}
}
\examples{
# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

alignGal <- getAlignGal(bamFiles[1], reverseComplement=TRUE, genomeBuild="mm9")

alignGR <- as(alignGal, "GRanges")

alignGRList <- GRangesList(as.list(split(alignGR, seqnames(alignGR))))

binSize <- 1e5  # use a large fixed bin size for demo only

binGR <- binCount(alignGRList$chrX, binSize)

# test on GRanges object
nbhGR <- nbh(binGR, 2, runViterbi=TRUE)

# test on integer object
nbhList <- nbh(values(binGR)$count, 2, runViterbi=TRUE)








}