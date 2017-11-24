\name{galp2gal}
\alias{galp2gal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Convert GAlignmentPairs to GAlignments
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Convert GAlignmentPairs to GAlignments using CIGAR to mark flanked portion of the pairs as 'N'.
}
\usage{
galp2gal(galp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{galp}{
%%     ~~Describe \code{galp} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Each proper read pairs is combined into a single alignment record making use of the CIGAR flag `N' to indicate the number of bases between the mate pairs (i.e., the difference between the start of the right mate pair and the end of the left mate pair). In other words, the paired-end alignments are treated as gapped alignments of long fragments. The function is used within \code{\link{getAlignGal}} but can be used as a stand-alone function as well.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{gal}{GAlignments object containing for each paired alignments a single alignment record.}
}
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
\code{\link{getAlignGal}, \link{combineAlignGals}, \link[GenomicAlignments]{readGAlignments}, \link[GenomicAlignments]{readGAlignmentPairs}, \link[rtracklayer]{import}}
}
\examples{
library(Rsamtools)

extdata.dir <- system.file("extdata", package="RIPSeeker")

ex1_file <- list.files(extdata.dir, "ex1.bam", recursive=TRUE, full.names=TRUE)

galp <- readGAlignmentPairs(ex1_file, use.names=TRUE)

galp

gal <- galp2gal(galp)

gal

}
