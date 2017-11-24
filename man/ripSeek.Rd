\name{ripSeek}
\alias{ripSeek}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
HMM-based de novo RIP predictions using alignment data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is the main interface to most essential functions of RIPSeeker package.
}
\usage{
ripSeek(bamPath, cNAME, binSize = NULL, strandType = NULL, 
	paired=FALSE, biomaRt_dataset, goAnno, exportFormat = "gff3", 
	annotateFormat = "txt", annotateType = "TSS", outDir, 
	padjMethod = "BH", logOddCutoff = 0, pvalCutoff = 1, 
	pvalAdjCutoff = 1, eFDRCutoff = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamPath}{
%%     ~~Describe \code{bamPath} here~~
	Either a path to all of the bam files or a list of paths to individual BAM files. BED and SAM files are also accepted.
}
  \item{cNAME}{
%%     ~~Describe \code{cNAME} here~~
	An identifer pattern found in the control alignment files. Once specified, these files will be used as control and the remaining files as RIP for discriminative analysis (see \code{\link{seekRIP}}).
}
  \item{binSize}{
%%     ~~Describe \code{binSize} here~~
	Size to use for binning the read counts across each chromosome. If NULL, optimal bin size within a range (default: minBinSize=200, maxBinSize=1200) will be automatically selected (see \code{\link{selectBinSize}}).
}
  \item{strandType}{
%%     ~~Describe \code{strandType} here~~
	Type of strand can be +, -, or * as in GAlignments, GAlignmentPairs, or GRanges (see \code{\link{GenomicRanges}}).
}
  \item{paired}{
%%     ~~Describe \code{strandType} here~~
  Binary to indicate whether the library is paired-end (TRUE) or single-end (FALSE by default) (see \code{\link{getAlignGal}}).
}
  \item{biomaRt_dataset}{
%%     ~~Describe \code{biomaRt_dataset} here~~
	The dataset name used in biomaRt for retrieving genomic information for a given species name (see \code{\link{annotateRIP}}).
}
  \item{goAnno}{
%%     ~~Describe \code{goAnno} here~~
	GO dataset name used for GO enrichment analysis (See \code{\link{annotateRIP}}).
}
  \item{exportFormat}{
%%     ~~Describe \code{exportFormat} here~~
	Format to export the RIP predictions. The commonly used ones are GFF and BED, which can be directly imported as a track to a genomic viewer such as Integrative Genomic Viewer, SAVANT or USCSC browser.
}
  \item{annotateFormat}{
%%     ~~Describe \code{annotateFormat} here~~
	Format to export the annotated RIP predictions. The default "txt" is a tab-delimited format, recommanded for viewing in Excel.
}
  \item{annotateType}{
%%     ~~Describe \code{annotateType} here~~
	Type of genomic information in association with the RIP predictions that can be retrieved from Ensembl database (Default: TSS; See \code{\link{annotateRIP}}).
}
  \item{outDir}{
%%     ~~Describe \code{outDir} here~~
	Output directory to save the results. The output data include ...
}
  \item{padjMethod}{
%%     ~~Describe \code{padjMethod} here~~
	Method to adjust multiple testing (Benjamini-Hocherge method by default).
}
  \item{logOddCutoff}{
%%     ~~Describe \code{logOddCutoff} here~~
	Threshold for the  log odd ratio of posterior for the RIP over the background states (See \code{\link{seekRIP}}). Only peaks with logOdd score \emph{greater} than the \code{logOddCutoff} will be reported. Default: 1.
}
  \item{pvalCutoff}{
%%     ~~Describe \code{pvalCutoff} here~~
	Threshold for the p-value for the logOdd score. Only peaks with p-value \emph{less} than the \code{pvalCutoff} will be reported. Default: 1 (i.e. no cutoff).
}
  \item{pvalAdjCutoff}{
%%     ~~Describe \code{pvalAdjCutoff} here~~
	Threshold for the adjusted p-value for the logOdd score. Only peaks with adjusted p-value \emph{less} than the \code{pvalAdjCutoff} will be reported. Default: 1 (i.e. no cutoff).
}
  \item{eFDRCutoff}{
%%     ~~Describe \code{eFDRCutoff} here~~
	Threshold for the empirical false discovery rate (eFDR). Only peaks with eFDR \emph{less} than the \code{eFDRCutoff} will be reported. Default: 1 (i.e. no cutoff).
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
	Arguments passed to \code{\link{mainSeek}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This is the main front-end function of RIPSeeker and in many cases the only function that users need to get RIP predictions and all relevant information.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list is returned with the following items:
% \item{mainSeekOutputRIP}{A list comprising two item: (1) GRangesList of the HMM trained parameters for each chromosome on RIP; (2) GAlignments object of processed the RIP alignment inputs from \code{\link{combineAlignGals}}.}
% \item{mainSeekOutputCTL}{Same as \code{mainSeekOutputRIP} but for the control library.}

\item{mainSeekOutputRIP}{A (inner) list comprising three items:

	\code{nbhGRList}: \code{\link[GenomicRanges]{GRangesList}} of the HMM trained parameters for each chromosome on RIP.
  
	\code{alignGal, alignGalFiltered}: \code{\link[GenomicAlignments]{GAlignments}} objects of the RIP alignment outputs from \code{\link{combineAlignGals}} and \code{\link{disambiguateMultihits}}, respectively. The former may contain multiple alignments due to the same reads whereas the latter contains a one-to-one mapping from read to alignment after disambiguating the multihits.
}

\item{mainSeekOutputCTL}{Same as \code{mainSeekOutputRIP} but for the control library (if available).}

\item{RIPGRList}{The results as \code{\link[GenomicRanges]{GRangesList}} generated from the RIP peak detection. Each list item represents the RIP peaks on a chromosome accompanied with statistical scores including (read) count, logOddScore, pval, pvalAdj, eFDR for the RIP and control (if available). Please refer to \code{\link{seekRIP}} for more details.}

\item{annotatedRIPGR}{If \code{annotatedRIPGR} is TRUE, the additional genomic information will be retreived according to the genomic coordinates of the peaks in RIPGRList. The results are saved in this separate GRanges object as the final results that user will find the most useful.}
}
\references{
%% ~put references to the literature/web site here ~
Zhao, J., Ohsumi, T. K., Kung, J. T., Ogawa, Y., Grau, D. J., Sarma, K., Song, J. J., et al. (2010). Genome-wide Identification of Polycomb-Associated RNAs by RIP-seq. Molecular Cell, 40(6), 939D953. doi:10.1016/j.molcel.2010.12.011

The RIPSeeker manuscript has been submitted to NAR for review.
}
\author{
%%  ~~who you are~~
Yue Li
}
\note{
%  ~~further notes~~
You may only want to know the expression/abundance of known transcripts/genes or the foldchange between two conditions. In that case, use \code{\link{rulebaseRIPSeek}} and \code{\link{computeRPKM}}, respectively. Both singl-end and paired-end alignments are supported in these functions.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{rulebaseRIPSeek}}
}
\examples{
if(interactive()) { # need internet connection

# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

cNAME <- "SRR039214" 						# specify control name

# output file directory
outDir <- paste(getwd(), "ripSeek_example", sep="/")

# Parameters setting
binSize <- NULL							# automatically determine bin size
minBinSize <- 10000						# min bin size in automatic bin size selection
maxBinSize <- 12000						# max bin size in automatic bin size selection
multicore <- TRUE						# use multicore
strandType <- "-"						# set strand type to minus strand

biomart <- "ENSEMBL_MART_ENSEMBL"		# use archive to get ensembl 65
dataset <- "mmusculus_gene_ensembl"		# mouse dataset id name	
host <- "dec2011.archive.ensembl.org" 	# use ensembl 65 for annotation

goAnno <- "org.Mm.eg.db"


################ run main function ripSeek to predict RIP ################
seekOut <- ripSeek(bamPath=bamFiles, cNAME=cNAME, 
		binSize=binSize, minBinSize = minBinSize, 
		maxBinSize = maxBinSize, strandType=strandType, 
		outDir=outDir, silentMain=FALSE,
		verbose=TRUE, reverseComplement=TRUE, genomeBuild="mm9",
		biomart=biomart, host=host,
		biomaRt_dataset = dataset,
		goAnno = goAnno,
		uniqueHit = TRUE, assignMultihits = TRUE, 
		rerunWithDisambiguatedMultihits = TRUE, multicore=multicore)


################ visualization ################

viewRIP(seekOut$RIPGRList$chrX, seekOut$mainSeekOutputRIP$alignGalFiltered, 
	seekOut$mainSeekOutputCTL$alignGalFiltered, scoreType="eFDR")

}

}
