\name{mainSeek}
\alias{mainSeek}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Train HMM paramters on each chromosome independently from the alignments.
}
\description{
A back-end function used by the front-end function \code{\link{ripSeek}} to train HMM paramters on all of the chromosomes indepdently. This function in turn calls another function \code{mainSeekSingleChrom} to compute HMM paramters on each chromosome separately or in parallel (if multicore is TRUE).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mainSeek(bamFiles, reverseComplement = FALSE, 
	genomeBuild = "mm9", uniqueHit = TRUE, 
	assignMultihits = TRUE, strandType = NULL, 
  paired=FALSE, rerunWithDisambiguatedMultihits = TRUE,
 	silentMain = FALSE, multicore = TRUE, 
 	returnAllResults = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamFiles}{
%%     ~~Describe \code{bamFiles} here~~
	A list of paths to individual BAM files. BED and SAM files are also accepted.
}
  \item{reverseComplement}{
%%     ~~Describe \code{reverseComplement} here~~
	Whether the reads came from the original or the opposite strand of the RNA being sequenced. If former, then \code{reverseComplement} should be FALSE; otherwise TRUE, in which case the strand signs will be switched from + to -, - to +, and * is unchanged.
}
  \item{genomeBuild}{
%%     ~~Describe \code{genomeBuild} here~~
	When the input alignment format is BED, genomeBuild is only required in \code{\link{getAlignGal}} to determine the chromosome lengths for the GAlignments obejct using function \code{\link{SeqinfoForUCSCGenome}}. BAM and SAM header have chromosome information, and thus genomeBuild is not needed.
}
  \item{uniqueHit}{
%%     ~~Describe \code{uniqueHit} here~~
	Binary indicator. If uniqueHit=TRUE, only reads mapped to single unique loci are used to train the HMM. Otherwise, all of the reads including multihits will be used for the HMM. A multihit is a read mapped to more than one loci. The flags for uniqueHits and multihits are the metadata values of GAlignments object constructed in \code{\link{getAlignGal}}.
}
  \item{assignMultihits}{
%%     ~~Describe \code{assignMultihits} here~~
	Binary indicator used by ripSeek to tell the function whether disambiguate multihits by assigning them to unique loci with the maximum posterior probability obtained from running HMM (See \code{\link{nbh_em}})
}
  \item{strandType}{
%%     ~~Describe \code{strandType} here~~
	A character variable indicate which strand the RIPSeeker needs to operate on. The options are NULL, '+', '-', '*'. If NULL or '*', then all of the reads will be used (preferable for non-strand specific sequencing). If '+' or '-', only reads from '+' or '-' strand will be used, respectively. Note that the sign is assumed to be THE SAME AS the strand sign of the processed alignment object and will be the opposite sign if reverseComplement is TRUE (See \code{reverseComplement} above).
}
  \item{paired}{
%%     ~~Describe \code{paired} here~~
  Binary to indicate whether the library is paired-end (TRUE) or single-end (FALSE by default) (see \code{\link{getAlignGal}}).
}
  \item{rerunWithDisambiguatedMultihits}{
%%     ~~Describe \code{rerunWithDisambiguatedMultihits} here~~
	After multihits have been asigned to unique loci, \code{rerunWithDisambiguatedMultihits} (Default: TRUE) indicates whether to re-run the HMM on the augmented read alignmnet data. If FALSE, the HMM step will not be re-run, and the workflow will proceed to RIP detection (See \code{\link{seekRIP}}) using the nondisambiguated alignments, which can either be the alignments containing only the uniqueHits (if uniqueHit=TRUE) or the alignments containing both the uniqueHits and multiHits (if uniqueHit=FALSE).
}
  \item{silentMain}{
%%     ~~Describe \code{silentMain} here~~
	Binary indicator to indicate whether to disable the verbose from the \code{\link{mainSeekSingleChrom}} function. If FALSE (by default), the EM training process will be output to the console for user to keep track of the training progress.
}
  \item{multicore}{
%%     ~~Describe \code{multicore} here~~
	Binary indicator to indicate whether to use \code{\link[parallel]{mclapply}}function to compute HMM on chromosomes in parallel. The multicore function will speed up the computation by a factor proportional to the total number of CPU cores on the machine but may impose larger memory overhead than the singe-threading approach.
}
  \item{returnAllResults}{
%%     ~~Describe \code{returnAllResults} here~~
	Binary indicator to indicate whether to return all (HMM trained parameters, original, and disambiguated GAlignments) or just the HMM results.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
	Arguments passed to \code{\link{mainSeekSingleChrom}}.
}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
	A list containing:
	\item{nbhGRList}{GRangesList each item containig the HMM training results on a single chromosome.}
	\item{alignGal}{Original alignment data in GAlignments object}
	\item{alignGalFiltered}{Disambiguated alignmnet data with multihits assigned to unique loci.}
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
	\code{\link{ripSeek}, \link{mainSeekSingleChrom}, \link[parallel]{mclapply}}
}
\examples{
# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

# Parameters setting
binSize <- 1e5							  # use a large fixed bin size for demo only
minBinSize <- NULL						# min bin size in automatic bin size selection
maxBinSize <- NULL						# max bin size in automatic bin size selection
multicore <- FALSE						# use multicore
strandType <- "-"							# set strand type to minus strand

################ run mainSeekSingleChrom function for HMM inference on all chromosomes ################
mainSeekOut <- mainSeek(bamFiles=grep(pattern="SRR039214", 
    bamFiles, value=TRUE, invert=TRUE),
		binSize=binSize, minBinSize = minBinSize, 
		maxBinSize = maxBinSize, strandType=strandType, 		
		reverseComplement=TRUE, genomeBuild="mm9",
		uniqueHit = TRUE, assignMultihits = TRUE, 
		rerunWithDisambiguatedMultihits = FALSE,				
		multicore=multicore, silentMain=FALSE, verbose=TRUE)

}
