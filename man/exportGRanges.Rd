\name{exportGRanges}
\alias{exportGRanges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Export GRanges object in a specified format
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A wrapper function of \code{\link[rtracklayer]{export}} with additional support for exporting tab-delimted format with no re-arrangement of the original GRanges output.
}
\usage{
exportGRanges(gRanges, outfile, exportFormat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gRanges}{
%%     ~~Describe \code{gRanges} here~~
	\code{\link[GenomicRanges]{GRanges}} object to export.
}
  \item{outfile}{
%%     ~~Describe \code{outfile} here~~
	File path for output.
}
  \item{exportFormat}{
%%     ~~Describe \code{exportFormat} here~~
	Desirable format including "txt" and other formats specified in \code{\link[rtracklayer]{export}}.
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
Output the text to the file stream defined in \code{outfile}.
}
\references{
%% ~put references to the literature/web site here ~
Michael Lawrence, Vince Carey and Robert Gentleman. rtracklayer: R interface to genome browsers and their annotation tracks. R package version 1.16.3.
}
\author{
%%  ~~who you are~~
Yue Li
}
\note{
%%  ~~further notes~~
The function is used in ripSeek to export desired format and can be used as general purpose function.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link[rtracklayer]{export}}
}
\examples{
if(interactive()) { # need permission to write to the current dir
gr <-
    GRanges(seqnames =
            Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
            ranges =
            IRanges(1:10, width = 10:1, names = head(letters,10)),
            strand =
            Rle(strand(c("-", "+", "*", "+", "-")),
                c(1, 2, 2, 3, 2)),
            score = 1:10,
            GC = seq(1, 0, length=10))

	outfile <- paste(getwd(), "/gr.txt", sep="/")

	exportGRanges(gr, outfile=outfile, exportFormat="txt")
}

}