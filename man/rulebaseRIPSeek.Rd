\name{rulebaseRIPSeek}
\alias{rulebaseRIPSeek}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compute the RPKM and foldchange between two conditions for the annotated genes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function takes alignments in two conditions (with replicates) as input and computes the gene expression in the two conditions in the unit of RPKM (reads per kilobase of exon per million mapped reads) or FPKM for paired-end alignments (where "F" stands for the fragment the mate-pairs are derived from), and then the foldchange ratio of the RPKM of each gene in RIP or treatment condition in general over control condition. The control files (i.e. the denominator in the foldchange ratio) is specified by user in the "cNAME" argument.
}
\usage{
rulebaseRIPSeek(bamFiles, cNAME, featureGRanges, rpkmCutoff = 0.4, 
	fcCutoff = 3, moreRIPGeneInfo = TRUE, idType = "ensembl_transcript_id", 
	myMin = .Machine$double.xmin, saveRData, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamFiles}{
%%     ~~Describe \code{bamFiles} here~~
	A list of one or more BAM/SAM/BED alignment files.
}
  \item{cNAME}{
%%     ~~Describe \code{cNAME} here~~
	An identifer pattern found in the control alignment files. Once specified, these files will be used as control as the denomenator of the foldchange ratio and the remaining files as RIP, the numarator of the foldchange ratio.
}
  \item{featureGRanges}{
%%     ~~Describe \code{featureGRanges} here~~
	GRanges of features as an optional argument for function to compute RPKM/FPKM just for those features without retrieving online annotations.
}
  \item{rpkmCutoff}{
%%     ~~Describe \code{rpkmCutoff} here~~
	Cutoff for RPKM in RIP above which the genes will be reported if the \code{fcCutoff} is also satisfied (Default: 0.4).
}
  \item{fcCutoff}{
%%     ~~Describe \code{fcCutoff} here~~
	Cutoff for foldchange in RIP relative to the control above which the genes will be reported if the \code{rpkmCutoff} is also satisfied (Default: 3).
}
  \item{moreRIPGeneInfo}{
%%     ~~Describe \code{moreRIPGeneInfo} here~~
	Binary indicator to indicate whether to download more information for each genes/transcripts rather than having only the gene/transcript IDs (Default: TRUE).
}
  \item{idType}{
%%     ~~Describe \code{idType} here~~
	A character string that specifies the type of the annotations, which can "ensembl_transcript_id" (Default), "ensembl_gene_id", "ucsc", etc. Refer to \code{\link[biomaRt]{listFilters}} for more information.
}
  \item{myMin}{
%%     ~~Describe \code{myMin} here~~
	Add a small value to both the numerator and denomenator as "pseudocount" to prevent the case where the denomenator is zero and the ratio becomes infinity regardless the value of the numerator (Default: \code{.Machine$double.xmin}).
}
  \item{saveRData}{
%%     ~~Describe \code{saveRData} here~~
	Path of output RData and tab-delim results.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
	Extra arguments passed to \code{\link{computeRPKM}} and/or\code{\link[biomaRt]{useMart}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
	The function uses \code{\link{computeRPKM}} to download annotation and compute RPKM/FPKM of the annotated genes in the list of files. The alignments file are separated into control as identified by the "cNAME" and the RIP (or any treatment) that do not have the \code{cNAME} in their file names. The alignments in either group are pooled together. If \code{moreRIPGeneInfo} is specified, the function witll query the Ensembl database. The chromosome ID in the numerical format used in Ensembl is prefixed with "chr" and the strand 1 and -1 converted to + and - for convenience.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list containing the following items:
\item{nRPKM}{RPKM of genes in RIP or treatment condition ('n' stands for numerator in the foldchange ratio).}
\item{dRPKM}{RPKM of genes in control condition ('d' stands for denomenator in the foldchange ratio)}
\item{rpkmDF}{Data frame containing read count, RPKM for the RIP (or treatment) and control, foldchange, and optional gene information including "chromosome_name", "start_position", "end_position", "strand", "external_gene_id", "ensembl_transcript_id", "ensembl_gene_id", "ucsc", "description"}
\item{rpkmCutoff}{Cutoff used for RPKM as book keeping value.}
\item{fcCutoff}{Cutoff used for foldchange as book keeping value.}
\item{featureGRanges}{\code{\link[GenomicRanges]{GRanges}} object of the features for which the RPKM and foldchange are computed.}
}
\references{
%% ~put references to the literature/web site here ~
Zhao, J., Ohsumi, T. K., Kung, J. T., Ogawa, Y., Grau, D. J., Sarma, K., Song, J. J., et al. (2010). Genome-wide Identification of Polycomb-Associated RNAs by RIP-seq. Molecular Cell, 40(6), 939D953. doi:10.1016/j.molcel.2010.12.011


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
if(interactive()) {
	
# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

cNAME <- "SRR039214" 						# specify control name

# output file directory
outDir <- paste(getwd(), "ripSeek_example")

# use biomart
txDbName <- "biomart"
biomart <- "ENSEMBL_MART_ENSEMBL"		# use archive to get ensembl 65
dataset <- "mmusculus_gene_ensembl"		
host <- "dec2011.archive.ensembl.org" 	# use ensembl 65 for annotation

resultlist <- rulebaseRIPSeek(bamFiles, "SRR039214", dataset=dataset, 
		txDbName=txDbName, biomart=biomart, host=host, by="tx")

}


}
