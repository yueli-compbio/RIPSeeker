\name{combineAlignGals}
\alias{combineAlignGals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Combine alignment files into a single GAlignments object 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Import and process individual BAM/SAM/BED alignment files using \code{\link{getAlignGal}} and combine them into a single GAlignments.
}
\usage{
combineAlignGals(bamFiles, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamFiles}{
%%     ~~Describe \code{bamFiles} here~~
	A list of paths to the alignment files.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
	Arguments passed to \code{\link{getAlignGal}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
	If there is only one BAM file, then simply return the output from \code{\link{getAlignGal}}; otherwise, all processed alignments are pooled to form a single GAlignments object.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
	\item{combinedGal}{GAlignments object containning the (combined) processed alignments with the \code{values} slot saved for the "uniqueHits" binary flag defined in \code{\link{getAlignGal}} and \code{metadata} saved as a list containing argument setting for \code{reverseComplement, returnDuplicate, flagMultiHits, returnOnlyUniqueHits} defined in \code{\link{getAlignGal}}}.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
%%  ~~who you are~~
Yue Li
}
\note{
%%  ~~further notes~~
User are recommanded to pool technical replicates but keep biological replicate separate for confirmation.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{getAlignGal}, \link[GenomicAlignments]{readGAlignments}, \link[GenomicAlignments]{readGAlignmentPairs}, \link[rtracklayer]{import}}
}
\examples{
# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

# combine the alignments for technical replicates
alignGal <- 
combineAlignGals(bamFiles=grep(pattern="SRR039214", 
                  bamFiles, value=TRUE, invert=TRUE), 
                  reverseComplement=TRUE, genomeBuild="mm9")


}
