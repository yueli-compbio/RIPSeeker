\name{computeRPKM}
\alias{computeRPKM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compute RPKM based on gene annotations
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given a list of single-end or paired-end read alignment files in BAM/SAM/BED format, compute the read counts and normalized read counts as expression of annotated transcript in the unit of "reads per kilobase of exon per million mapped reads" (RPKM).
}
\usage{
computeRPKM(bamFiles, RIPSeekerRead = TRUE, paired = FALSE, 
	countMode = "IntersectionNotEmpty", featureGRanges, 
	idType = "ensembl_transcript_id", featureType = "exon", 
	ignore.strand = FALSE, txDbName = "biomart", 
	moreGeneInfo = FALSE, saveData, justRPKM = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamFiles}{
%%     ~~Describe \code{bamFiles} here~~
	A list of one or more BAM/SAM/BED alignment files.
}
  \item{RIPSeekerRead}{
%%     ~~Describe \code{RIPSeekerRead} here~~
	Binary flag. If TRUE, then import and process the alignment files using the built-in function \code{\link{combineAlignGals}} from RIPSeeker package; if FALSE, then import the files by directly calling the required functions. The flag makes using the function outside of \code{\link{RIPSeeker}} package become possible.
}
  \item{paired}{
%%     ~~Describe \code{paired} here~~
	Binary to indicate whether the alignments files are paired-end. The alignments file must be either paired-end or single-end but not both.
}
  \item{countMode}{
%%     ~~Describe \code{countMode} here~~
	An argument used to set the \code{mode} argument in the underlying function \code{\link{summarizeOverlaps}} employed to compute the read counts for each feature. The possible mode includes "Union", "IntersectionStrict", and "IntersectionNotEmpty". All three modes avoid double counting the reads by either discarding reads that completely fall into multiple features or counting the read only once for the feature that uniquely and completely includes it. Please refer to \code{\link{summarizeOverlaps}} for details. 
}
  \item{featureGRanges}{
%%     ~~Describe \code{featureGRanges} here~~
	GRanges of features as an optional argument for function to compute RPKM/FPKM just for those features without retrieving online annotations.
}
  \item{idType}{
%%     ~~Describe \code{idType} here~~
	A character string that specifies the type of the annotations, which can "ensembl_transcript_id", "ensembl_gene_id", "ucsc", etc. Refer to \code{\link[biomaRt]{listFilters}} for more information.
}
  \item{featureType}{
%%     ~~Describe \code{featureType} here~~
	Features that will be groupped by genes/transcripts in a GRangesList. The available options are "exon" (Default), "intron", "fiveUTR", "threeUTR", and "CDS" corresponding to the functions \code{\link[GenomicFeatures]{exonsBy}}, \code{\link[GenomicFeatures]{cdsBy}}, \code{\link[GenomicFeatures]{intronsByTranscript}}, \code{\link[GenomicFeatures]{fiveUTRsByTranscript}}, \code{\link[GenomicFeatures]{threeUTRsByTranscript}}, and \code{\link[GenomicFeatures]{cdsBy}}, respectively.
}
  \item{ignore.strand}{
%%     ~~Describe \code{ignore.strand} here~~
	Whether to ignore strand when counting the reads (Default: FALSE).
}
  \item{txDbName}{
%%     ~~Describe \code{txDbName} here~~
	Name of the transcript database to use to retreive the annotation. The available options are "biomart" (Default) or "UCSC" corresponding to the functions \code{\link[GenomicFeatures]{makeTxDbFromBiomart}} and \code{\link[GenomicFeatures]{makeTxDbFromUCSC}}, respectively.
}
  \item{moreGeneInfo}{
%%     ~~Describe \code{moreGeneInfo} here~~
	Binary indicator to indicate whether to download more information for each genes/transcripts rather than having only the gene/transcript IDs (Default: FALSE).
}
  \item{saveData}{
%%     ~~Describe \code{saveData} here~~
	Path of output file.
}
  \item{justRPKM}{
%%     ~~Describe \code{justRPKM} here~~
	Binary for whether to return only the \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}}.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
	Extra arguments passed to functions \code{\link[GenomicFeatures]{makeTxDbFromBiomart}}, \code{\link[GenomicFeatures]{makeTxDbFromUCSC}}, \code{\link[biomaRt]{useMart}}, \code{\link{combineAlignGals}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
	The function is a wrapper function making use of several external functions from several well maintained and freely available Bioconductor packages including GenomicFeatures, GenomicRanges, biomaRt and Rsamtools packages. The paired-end alignments are converted into single-end using function \code{\link{galp2gal}} and then subject to read count computation by \code{\link[GenomicAlignments]{summarizeOverlaps}}, which does not yet directly support paired-end alignments.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{rpkmSEobject}{A \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}} object with \code{assays} slot saved for counts, \code{rowRanges} holds the features, \code{metadata} for RPKM/FPKM (normalized) gene expression.}
\item{rpkmDF}{Data frame with or without the detailed gene information columns depending on whether \code{moreGeneInfo} is TRUE or FALSE. \code{rpkmDF} is only returned within in a list when \code{justRPKM} is FALSE.}
\item{featureGRanges}{The features in GRanges object that are used to compute the gene expression. \code{featureGRanges} is only returned within in a list when \code{justRPKM} is FALSE.}
}
\references{
%% ~put references to the literature/web site here ~
M. Carlson, H. Pages, P. Aboyoun, S. Falcon, M. Morgan, D. Sarkar and M. Lawrence. GenomicFeatures: Tools for making and manipulating transcript centric annotations. R package version 1.8.2.

P. Aboyoun, H. Pages and M. Lawrence (). GenomicRanges: Representation and manipulation of genomic intervals. R package version 1.8.9.

Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. Steffen Durinck, Paul T. Spellman, Ewan Birney and Wolfgang Huber, Nature Protocols 4, 1184-1191 (2009).

BioMart and Bioconductor: a powerful link between biological databases and microarray data analysis. Steffen Durinck, Yves Moreau, Arek Kasprzyk, Sean Davis, Bart De Moor, Alvis Brazma and Wolfgang Huber, Bioinformatics 21, 3439-3440 (2005).

Martin Morgan and Herv\'e Pag\`es (). Rsamtools: Binary alignment (BAM), variant call (BCF), or tabix file import. R package version 1.8.5. http://bioconductor.org/packages/release/bioc/html/Rsamtools.html
}
\author{
%%  ~~who you are~~
Yue Li
}
\note{
%%  ~~further notes~~
Also works for RNA-seq alignments. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link[GenomicFeatures]{makeTxDbFromBiomart}, \link[GenomicFeatures]{makeTxDbFromUCSC}, \link[biomaRt]{useMart}, \link[GenomicFeatures]{exonsBy}, \link[GenomicFeatures]{cdsBy}, \link[GenomicFeatures]{intronsByTranscript}, \link[GenomicFeatures]{fiveUTRsByTranscript}, \link[GenomicFeatures]{threeUTRsByTranscript}, \link[GenomicFeatures]{cdsBy}, \link{combineAlignGals}, \link[GenomicAlignments]{summarizeOverlaps}, \link[Rsamtools]{ScanBamParam}, \link[GenomicAlignments]{readGAlignmentPairs}, \link[GenomicAlignments]{readGAlignments} }
}
\examples{
if(interactive()) { # need internet connection
	
# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

# use biomart
txDbName <- "biomart"
biomart <- "ENSEMBL_MART_ENSEMBL"		# use archive to get ensembl 65
dataset <- "mmusculus_gene_ensembl"		
host <- "dec2011.archive.ensembl.org" 	# use ensembl 65 for annotation

resultlist <- computeRPKM(bamFiles=grep(pattern="SRR039214", 
        bamFiles, value=TRUE, invert=TRUE), #featureGRanges=featureGRanges, 
				dataset=dataset, moreGeneInfo=TRUE, justRPKM=FALSE,
				idType="ensembl_transcript_id", txDbName=txDbName, 
				biomart=biomart, host=host, by="tx")
		
}

}
