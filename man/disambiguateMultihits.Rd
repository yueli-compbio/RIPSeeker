\name{disambiguateMultihits}
\alias{disambiguateMultihits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Assign each multihit to a unique region based on the posterior for the read-enriched hidden state
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Among multiple alignments of the same read (i.e. multihit), select the alignment corresponding to the bin with the maximum posterior for the enriched state.
}
\usage{
disambiguateMultihits(alignGal, nbhGRList, postprobCutoff = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignGal}{
%%     ~~Describe \code{alignGal} here~~
	GAlignments object with an additional column in the \code{values} slot that indicates whether the read corresponding to the current alignment is a unique hit (i.e., read mapped uniquely to a single loci) or multihit (i.e., read mapped to multiple loci).
}
  \item{nbhGRList}{
%%     ~~Describe \code{nbhGRList} here~~
	GRangesList each item containig the HMM training results on a single chromosome. Importantly, the posterior probabilities for the background and enriched states need to be present the \code{metadata} slot and used to disambiguate multihits, which is done by \code{\link{mainSeekSingleChrom}}.
}
  \item{postprobCutoff}{
%%     ~~Describe \code{postprobCutoff} here~~
	Posterior cutoff for returning only the reads with maximum posterior that is greater than the threshold (Default: 0; i.e., no cutoff).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Each multihit (i.e., read aligned to multiple loci) flagged in the \code{\link{getAlignGal}} function are assigned to a unique locus corresponding to the \eqn{j^th} bin with the highest posterior or responsibility from the RIP state. Intuitively, the RIP state corresponds to the read-enriched loci. Disambiguating multihits in this way will potentially improve the power of detecting more RIP regions but may also introduce certain bias towards the idea of "rich gets richer". After this step, RIPSeeker will rerun the functions from \code{\link{selectBinSize}} to \code{\link{nbh}} to improve the HMM model estimation with augmented read count data. Optionally, user can choose not to reiterate the training process to go straight to the next step to detect RIP regions (See \code{\link{seekRIP}}).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
GAlignments with each read mapped uniquely to a single locus.
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
\code{\link{getAlignGal}, \link{ripSeek}, \link{mainSeek}, \link{mainSeekSingleChrom}}
}
\examples{
# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

# Parameters setting
binSize <- 1e5							  # use a large fixed bin size for demo only
minBinSize <- NULL						# turn off min bin size in automatic bin size selection
maxBinSize <- NULL						# turn off max bin size in automatic bin size selection
multicore <- FALSE						# use multicore
strandType <- "-"						  # set strand type to minus strand

# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

alignGal <- combineAlignGals(bamFiles=grep(pattern="SRR039214",             
            bamFiles, value=TRUE, invert=TRUE), reverseComplement=TRUE, genomeBuild="mm9")

alignGR <- as(alignGal, "GRanges")

alignGR <- addPseudoAlignment(alignGR)

alignGRList <- GRangesList(as.list(split(alignGR, seqnames(alignGR))))

################ run mainSeekSingleChrom function for HMM inference on a single chromosome ################
nbhGRList <- lapply(alignGRList, mainSeekSingleChrom, K = 2, binSize=binSize, 
			
			minBinSize = minBinSize, maxBinSize = maxBinSize, runViterbi=FALSE)

nbhGRList <- GRangesList(nbhGRList)

alignGalFiltered <- disambiguateMultihits(alignGal, nbhGRList)

}
