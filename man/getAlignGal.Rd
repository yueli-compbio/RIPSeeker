\name{getAlignGal}
\alias{getAlignGal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Import and processs in BAM/SAM/BED format
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Import and process single-end or paired-end alignments in a BAM/SAM/BED file to retain valid alignments defined by the arguments below. Multihits (same read mapped to multiple loci) are flagged for the subsequent disambiguation with function \code{\link{disambiguateMultihits}}). The final output is a GAlignments object.
}
\usage{
getAlignGal(alignFilePath, format, genomeBuild, 
	deleteGeneratedBAM = FALSE, reverseComplement = FALSE, 
	returnDuplicate = FALSE, flagMultiHits = TRUE, 
	returnOnlyUniqueHits = FALSE, paired = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignFilePath}{
%%     ~~Describe \code{alignFilePath} here~~
	Path to the alignment file.
}
  \item{format}{
%%     ~~Describe \code{format} here~~
	The alignmnet format can be determiend automatically from the file extension or specified by the user. The supported formats are BAM, SAM, and BED.
}
  \item{genomeBuild}{
%%     ~~Describe \code{genomeBuild} here~~
	Genome build used to obtain the chromosome information from online UCSC database in order to construct GAlignments object. Since the BAM/SAM header provides the chromosome information, the argument needs to be set only in the absence of the header information for some BAM/SAM files or when BED file is used. Examples for the common \code{genomeBuild} are "mm9" for mouse or "hg19" for human reference genomes. Note that an appropriate genome build that has been used in the alignment is important for desirable outcome. For instance, user should use "mm10" if the alignments are based on "mm10" rather than "mm9" genome build.
}
  \item{deleteGeneratedBAM}{
%%     ~~Describe \code{deleteGeneratedBAM} here~~
	Binary indicator to indicate whether the converted BAM from the original SAM input file needs to be deleted from the local disk (Default: FALSE).
}
  \item{reverseComplement}{
%%     ~~Describe \code{reverseComplement} here~~
	Binary indicator to indicate whether the reads were sequenced from the opposite strand of the original RNA molecule. \code{reverseComplement} only applies to strand-specific sequencing in which case only the strand generated during \emph{second strand synthesis} is sequenced. Thus, if \code{reverseComplement=TRUE}, the strand signs of the alignments are switched (i.e. + to -, - to +, and * unchanged); otherwise (\code{reverseComplement=FALSE}) retian the original the strand signs.
}
  \item{returnDuplicate}{
%%     ~~Describe \code{returnDuplicate} here~~
	Indicator (TRUE, FALSE, NA) to instruct whether the duplicate alignmnets need to be returned (Default: FALSE). Duplicate reads are a set of reads that align to exactly the same genomic coordinate. Because transcripts are usually hundreds or thousands of base pairs long and thus much longer than the read (25-100 nt), the chance that the same 25-100 nt portion of the transcript being sequenced multiple times is very small and may very likely be due to PCR artifact. This argument is acutally passed to 'isDuplicate' in \code{\link[Rsamtools]{scanBamFlag}}.
}
  \item{flagMultiHits}{
%%     ~~Describe \code{flagMultiHits} here~~
	Binary indicator for whether to add additional binary column named "uniqueHits" to indicate whether the corresponding aligned reads are unique hit (uniqueHits==TRUE)  or multihit (uniqueHits==FALSE). Multihits represent multiple alignments of the same read due to gene duplications or repetitive elements of the genome. The multhits typically constitute a substantial proportion of the total mapped reads. Rather than being removed, these multihits are flagged (\code{flagMultiHits=TRUE} by default) and in the later step assigned to a unique region by (\code{\link{disambiguateMultihits}}).
}
  \item{returnOnlyUniqueHits}{
%%     ~~Describe \code{returnOnlyUniqueHits} here~~
	Binary indicator to return only the unique hits and discard all of the multihits (Default: FALSE).
}
  \item{paired}{
%%     ~~Describe \code{paired} here~~
	Binary indicator to indicate whether the alignments are paired-end (Default: FALSE). For paired-end alignments, properly paired reads are combined into a single alignment record making use of the CIGAR flag `N' to indicate the number of bases between the mate pairs (i.e., the length of the insert fragment). In other words, the paired-end alignments are treated as gapped alignments of long fragments (See \code{\link{galp2gal}}).
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
	Extra arguments are ignored.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
	The BAM file is imported using \code{\link[GenomicAlignments]{readGAlignments}} for single-end or \code{\link[GenomicAlignments]{readGAlignmentPairs}} for paired-end alignments. The SAM file is converted to BAM first and then imported as above. The BED file is first imported by \code{\link[rtracklayer]{import}} as GRanges object and subsequently converted to \code{GAlignments} via the constructor function \code{\link[GenomicAlignments]{GAlignments}}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
	\item{alignGal}{GAlignments object containning the processed alignments with the \code{values} slot saved for the "uniqueHits" binary flag (See \code{flagMultiHits} above) and \code{metadata} saved as a list containing argument setting for \code{reverseComplement, returnDuplicate, flagMultiHits, returnOnlyUniqueHits}	
}
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


}
