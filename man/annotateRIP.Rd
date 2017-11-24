\name{annotateRIP}
\alias{annotateRIP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Annotate RIP peaks with genomic information and perform GO enrichement
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given the genomic coordinates of each predicted RIP regions, query the Ensembl database whether each region is nearby or overlaps any known (noncoding) genes.
}
\usage{
annotateRIP(sigGRanges, biomaRt_dataset, featureType = "TSS", 
	goAnno, strandSpecific = FALSE, exportFormat = "txt", 
	hasGOdb = !missing(goAnno), goPval = 0.1, outDir, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sigGRanges}{
%%     ~~Describe \code{sigGRanges} here~~
	\code{\link[GenomicRanges]{GRanges}} object indicating the chromosomal coordinates of each RIP peaks.
}
  \item{biomaRt_dataset}{
%%     ~~Describe \code{biomaRt_dataset} here~~
	Ensembl dataset available from biomaRt (See \code{\link[biomaRt]{listDatasets}}). For instance, the human and mouse annotations are \code{hsapiens_gene_ensembl} and \code{mmusculus_gene_ensembl}, respectively.
}
  \item{featureType}{
%%     ~~Describe \code{featureType} here~~
	TSS, miRNA, Exon, 5'UTR, 3'UTR, transcript or Exon plus UTR defined in \code{\link[ChIPpeakAnno]{getAnnotation}}.
}
  \item{goAnno}{
%%     ~~Describe \code{goAnno} here~~
	Optional argugment that specifies a GO dataset used for GO enrichement analysis performed by \code{\link[ChIPpeakAnno]{getEnrichedGO}}. For instance, the human and mouse GO datasets are \code{org.Hs.eg.db} and \code{org.Mm.eg.db}.
}
  \item{strandSpecific}{
%%     ~~Describe \code{strandSpecific} here~~
	Indicate whether the annotations should be strand-specific (Default: FALSE)
}
  \item{exportFormat}{
%%     ~~Describe \code{exportFormat} here~~
	Format to export using \code{\link{exportGRanges}} (Default: "txt", i.e. tab-delim file).
}
  \item{hasGOdb}{
%%     ~~Describe \code{hasGOdb} here~~
	A binary flag that indicates whether GO enrichement is performed in order to export the results. \code{hasGOdb} can be FALSE either because \code{goAnno} is not specifiy or because the GO database does not exist.
}
  \item{goPval}{
%%     ~~Describe \code{goPval} here~~
	P-value cutoff to determine the significance of enriched GO terms by \code{\link[ChIPpeakAnno]{getEnrichedGO}}.
	
}
  \item{outDir}{
%%     ~~Describe \code{outDir} here~~
	Output directory.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
	Extra arguments passed to \code{\link[biomaRt]{useMart}} to specify the database and to passed \code{\link[ChIPpeakAnno]{getEnrichedGO}} to specify the GO enrichment procedure.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
To access the up-to-date Ensembl database, RIPSeeker employs \code{\link[biomaRt]{useMart}} and \code{\link[ChIPpeakAnno]{getAnnotation}} from biomaRt and ChIPpeakAnno Bioconductor packages to dynamically establish internet connection to the database and retrieve the up-to-date annotations. Then, \code{\link[ChIPpeakAnno]{annotatePeakInBatch}} from ChIPpeakAnno is used to efficiently annotate all of the predicted regions based on the Ensembl annotation. A predicted region may overlap multiple genes, all of which will be reported as separate records. Moreover, \code{\link[ChIPpeakAnno]{getEnrichedGO}} from ChIPpeakAnno is applied to the annotated predictions to discover enriched Gene Ontology (GO) terms involving the protein-associated transcriptome.

In order to use old annotation (e.g., mm9 v.s. mm10), user also needs to specify the host and biomart arguments accepted within \code{\link[biomaRt]{useMart}}. To access to mouse annotation from Ensembl version 65, for instance, user needs to call annotateRIP(..., dataset="mmusculus_gene_ensembl", biomart="ENSEMBL_MART_ENSEMBL", host="dec2011.archive.ensembl.org", ...), which will run useMart(dataset="mmusculus_gene_ensembl", biomart="ENSEMBL_MART_ENSEMBL", host="dec2011.archive.ensembl.org", ...) to get the mm9 annotation from Ensembl (v65).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

\item{sigGRangesAnnotated}{\code{sigGRanges} augmented with genomic information including "ensembl_gene_id", "external_gene_id", and "description"}

\item{enrichedGO}{Output from \link[ChIPpeakAnno]{getEnrichedGO}. All three main GO categories ("Biological Process", "Molecular Function", "Cellular Component") are combined together and returned. The argument is only returned when \code{hasGOdb} is TRUE.}

If \code{outDir} is specified, then the above \code{sigGRangesAnnotated} is saved as RIPregions_annotated.txt and RIPregions_annotated.RData, and \code{enrichedGO} as RIPregions_enrichedGO.txt in the \code{outDir} directory.

}
\references{
%% ~put references to the literature/web site here ~
Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. Steffen Durinck, Paul T. Spellman, Ewan Birney and Wolfgang Huber, Nature Protocols 4, 1184-1191 (2009).

BioMart and Bioconductor: a powerful link between biological databases and microarray data analysis. Steffen Durinck, Yves Moreau, Arek Kasprzyk, Sean Davis, Bart De Moor, Alvis Brazma and Wolfgang Huber, Bioinformatics 21, 3439-3440 (2005).
  
Lihua Julie Zhu, Herve Pages, Claude Gazin, Nathan Lawson, Jianhong Ou, Simon Lin, David Lapointe and Michael Green (2012). ChIPpeakAnno: Batch annotation of the peaks identified from either ChIP-seq, ChIP-chip experiments or any experiments resulted in large number of chromosome ranges.. R package version 2.4.0.
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
\code{\link[biomaRt]{useMart}, \link[ChIPpeakAnno]{getAnnotation}, \link[ChIPpeakAnno]{getEnrichedGO}}
}
\examples{
if(interactive()) { # need internet connection
# Retrieve system files
extdata.dir <- system.file("extdata", package="RIPSeeker") 

bamFiles <- list.files(extdata.dir, ".bam$", recursive=TRUE, full.names=TRUE)

bamFiles <- grep("PRC2", bamFiles, value=TRUE)

# Parameters setting
binSize <- NULL							# automatically determine bin size
minBinSize <- 10000						# min bin size in automatic bin size selection
maxBinSize <- 12000						# max bin size in automatic bin size selection
multicore <- FALSE						# use multicore
strandType <- "-"							# set strand type to minus strand

biomart <- "ENSEMBL_MART_ENSEMBL"		# use archive to get ensembl 65
dataset <- "mmusculus_gene_ensembl"		# mouse dataset id name	
host <- "dec2011.archive.ensembl.org" 	# use ensembl 65 for annotation

goAnno <- "org.Mm.eg.db"


################ run main function for HMM inference on all chromosomes ################
mainSeekOutputRIP <- mainSeek(
    bamFiles=grep(pattern="SRR039214", bamFiles, value=TRUE, invert=TRUE),
		binSize=binSize, minBinSize = minBinSize, 
		maxBinSize = maxBinSize, strandType=strandType, 		
		reverseComplement=TRUE, genomeBuild="mm9",
		uniqueHit = TRUE, assignMultihits = TRUE, 
		rerunWithDisambiguatedMultihits = TRUE,				
		multicore=multicore, silentMain=FALSE, verbose=TRUE)
		
# use defined binSize from RIP
RIPBinSize <- lapply(mainSeekOutputRIP$nbhGRList, function(x) median(width(x)))


mainSeekOutputCTL <- mainSeek(
    bamFiles=grep(pattern="SRR039214", bamFiles, value=TRUE, invert=FALSE),
		binSize=RIPBinSize, strandType=strandType, 		
		reverseComplement=TRUE, genomeBuild="mm9",
		uniqueHit = TRUE, assignMultihits = TRUE, 
		rerunWithDisambiguatedMultihits = TRUE,				
		multicore=multicore, silentMain=FALSE, verbose=TRUE)

################ significance test on Viterbi predicted peaks ################
ripGR <- seekRIP(mainSeekOutputRIP$nbhGRList$chrX, mainSeekOutputCTL$nbhGRList)


################ Annotate peaks ################

annotatedRIPGR <- annotateRIP(sigGRanges = ripGR,
				biomaRt_dataset = dataset, goAnno = goAnno, 
				strandSpecific = !is.null(strandType),
				host=host, biomart=biomart)

head(annotatedRIPGR$sigGRangesAnnotated)	
}
}