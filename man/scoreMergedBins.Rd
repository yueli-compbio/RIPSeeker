\name{scoreMergedBins}
\alias{scoreMergedBins}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Average log odd scores over bins being merged into a single region
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Sum, normalize the read counts, and average the logOdd score over the bins being merged into a single enriced region.
}
\usage{
scoreMergedBins(findOverlapsHits, unmergedGRAll, mergedGRAll)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{findOverlapsHits}{
%%     ~~Describe \code{findOverlapsHits} here~~
	Output from \code{\link{findOverlaps}} as two columns indices with the first column containing the indices for unmerged GRanges and the second column the indices of the matched merged GRanges.
}
  \item{unmergedGRAll}{
%%     ~~Describe \code{unmergedGRAll} here~~
	GRanges before merging.
}
  \item{mergedGRAll}{
%%     ~~Describe \code{mergedGRAll} here~~
	GRanges after merging.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The consecutive RIP-bins predicted by the Viterbi function (See \code{\link{nbh_vit}}) are merged into a single RIP region. An aggregate RIPScore as the averaged RIPScores over the associated merged bins is assigned to each merged RIP region. In the RIPSeeker workflow, the averaged RIPScore then becomes the representative score for the region and subject to significance test carried out in \code{\link{seekRIP}}.
}


\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A merged GRanges each with scores including summed read count, averaged log odd scores, and FPK (fragment per kilobase of region length). The latter score represent a normalized read count.
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
This function is expected to be called only from \code{\link{logScoreWithoutControl}} and \code{\link{logScoreWithControl}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{seekRIP}, \link{computeLogOdd} \link{logScoreWithControl}, \link{logScoreWithoutControl}}
}
\examples{
if(interactive()) { # see example in seekRIP
# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

# Parameters setting
binSize <- 1e5					# use a large fixed bin size for demo only
multicore <- FALSE			# use multicore
strandType <- "-"				# set strand type to minus strand

################ run main function for HMM inference on all chromosomes ################
mainSeekOutputRIP <- mainSeek(bamFiles=
    grep(pattern="SRR039214", bamFiles, value=TRUE, invert=TRUE),
		binSize=binSize, strandType=strandType, 		
		reverseComplement=TRUE, genomeBuild="mm9",
		uniqueHit = TRUE, assignMultihits = TRUE, 
		rerunWithDisambiguatedMultihits = TRUE,				
		multicore=multicore, silentMain=FALSE, verbose=TRUE)
		

nbhGRRIP <- mainSeekOutputRIP$nbhGRList$chrX

logOddScore <- computeLogOdd(nbhGRRIP)

values(nbhGRRIP) <- cbind(as.data.frame(values(nbhGRRIP)), logOddScore)
	
enrichIdx <- which(values(nbhGRRIP)$viterbi_state == 2)

unmergedRIP <- nbhGRRIP[enrichIdx]	
	
mergedRIP <- reduce(unmergedRIP, min.gapwidth = median(width(unmergedRIP) ))
	
overlapIdx <- findOverlaps(mergedRIP, unmergedRIP)

# a list with query hits as names and subject hits as items
findOverlapsHits <- split(overlapIdx, queryHits(overlapIdx))

# get the score for the first merged region
x <- scoreMergedBins(findOverlapsHits[[1]], unmergedRIP, mergedRIP)

# get scores for all of the merged regions
mergedRIPList <- lapply(split(overlapIdx, queryHits(overlapIdx)),
			
			scoreMergedBins, unmergedRIP, mergedRIP)
	
names(mergedRIPList) <- NULL
	
mergedRIP <- do.call(c, mergedRIPList)

# logOddScore is the averaged logOddScore across merged bins
mergedRIP
}
}