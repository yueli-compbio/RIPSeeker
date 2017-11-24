\name{addPseudoAlignment}
\alias{addPseudoAlignment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Add a psuedoalignment as a placeholder for the chromosome (Internal function)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Check whether chromosome has at least one alignment to prevent abnormal behaviour of the subsequent functions. In case no alignment is found on an entire chromosome, add a pseudo-alignment as a placeholder for that chromosome.
}
\usage{
addPseudoAlignment(alignGR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignGR}{
%%     ~~Describe \code{alignGR} here~~
GRanges object containing the alignment information.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
In case no alignment is found on an entire chromosome, add an alignment with start 1 and end 20 as a placeholder for the chromosome. This step is necessary to maintian the chromosome information.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{alignGR}{
	
	Original or augmented input GRanges object with pseudoreads, depending on whether there exists empty chromosome(s).		
}
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
	\code{\link{combineAlignGals}, \link[GenomicAlignments]{readGAlignments}, \link[GenomicAlignments]{readGAlignmentPairs}, \link[rtracklayer]{import}}
}
\examples{
# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

alignGal <- getAlignGal(bamFiles[1], reverseComplement=TRUE, genomeBuild="mm9")

alignGR <- as(alignGal, "GRanges")

alignGR

x <- addPseudoAlignment(alignGR)

x

}
