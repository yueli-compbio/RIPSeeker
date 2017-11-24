\name{combineRIP}
\alias{combineRIP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Combined predictions from (presumably) biological replicates.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A simple helper function that combines multiple prediction lists from biological replicates into a single list.
}
\usage{
combineRIP(ripPath, pattern="gff3$", 
	combineOption="intersect", 
	pvalCutoff=1, pvalAdjCutoff=1, eFDRCutoff=1, 
	logOddCutoff=-Inf, maxgap=1e3, minIntersect, genomeBuild)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{ripPath}{
%%     ~~Describe \code{x} here~~
	Path to predictions list in a select format as indicated by the file extension.
	}
	\item{pattern}{
%%     ~~Describe \code{x} here~~
	Pattern for the names of the prediction files to combine. The file names are expected to have a common extension such as "bed", "gff3", "gtf", but this is not enforced. Default: "gff3$" (i.e. the default output RIPregions.gff3 from \code{\link{ripSeek}}).
		
	}
	\item{combineOption}{
%%     ~~Describe \code{x} here~~
	Options on \emph{how} to combine the peaks including:
	
	\code{"intersect"}: is selected (default), only peaks in each biological replicate list that overlap with or are adjacent within \code{maxgap} nucleotides to at least \code{minIntersect} other replicates will be kept. If \code{minIntersect} is unspecified, then only the peaks that consistently predicted in all replicates are kept.
	
	\code{"merge"}: All overlapping peaks from the replicates will be merged into one peak.
	
	\code{"union"}: All overlapping peaks from the replicates will be merged into one peak.
	}
	\item{pvalCutoff}{
		Threshold for the p-value cutoff. Only peaks with p-value \emph{less} than the \code{logOddCutoff} will be reported. Default: 1 (i.e. no cutoff).
	}
	\item{pvalAdjCutoff}{
		Threshold for the adjusted p-value cutoff. Only peaks with adjusted p-value \emph{less} than the \code{logOddCutoff} will be reported. Default: 1 (i.e. no cutoff).
	}
	\item{eFDRCutoff}{
		Threshold for the empirical false discovery rate (eFDR). Only peaks with eFDR \emph{less} than the \code{eFDRCutoff} will be reported. Default: 1 (i.e. no cutoff).		
	}
	\item{logOddCutoff}{
		Threshold for the  log odd ratio of posterior for the RIP over the background states (See \code{\link{seekRIP}}). Only peaks with logOdd score \emph{greater} than the \code{logOddCutoff} will be reported. Default: -Inf (i.e. no cutoff).
	}
	\item{maxgap}{
		Maximum gap allowed to determine two peaks agree with each other.
	}
	\item{minIntersect}{		
		Mininum number of replicates required to have peaks either intersect or are adjacent to the peak in other replciate.
	}
	\item{genomeBuild}{
		Genome build used to obtain the chromosome information from online UCSC database to assign chromosome length to the \code{\link[GenomicRanges]{GRanges}} object created as the combined peak list.
	}	
}
% \details{
%%  ~~ If necessary, more details than the description above ~~
% }
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
	\item{gr}{\code{\link[GenomicRanges]{GRanges}} object containning chromosome locations of the combined peaks.}
}
\references{
%% ~put references to the literature/web site here ~
	P. Aboyoun, H. Pages and M. Lawrence. GenomicRanges: Representation
  and manipulation of genomic intervals. R package version 1.8.9.
  
	Michael Lawrence, Vince Carey and Robert Gentleman. rtracklayer: R
  interface to genome browsers and their annotation tracks. R package
  version 1.16.3.

}
\author{
%%  ~~who you are~~
Yue Li
}
\note{
%%  ~~further notes~~
Please run \code{\link{ripSeek}} first on all biological replicates and renamed each "RIPregions.gff3" output to correspond to different biological replicates and place all of the files into a single folder. The path of this folder can then be used as the input argument for \code{ripPath}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{combineAlignGals}, \link{ripSeek}, \link[rtracklayer]{import}, \link[rtracklayer]{import}, \link[IRanges]{reduce}, \link[IRanges]{countOverlaps}}
}
\examples{

# Retrieve system files
ripPath <- system.file("extdata/RIPregions", package="RIPSeeker")

gr1 <- combineRIP(ripPath, combineOption="intersect", genomeBuild="mm9")

gr2 <- combineRIP(ripPath, combineOption="merge", genomeBuild="mm9")

gr3 <- combineRIP(ripPath, combineOption="union", genomeBuild="mm9")

length(gr1)

length(gr2)

length(gr3)

}