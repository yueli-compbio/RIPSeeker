\name{nbh.GRanges}
\alias{nbh.GRanges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Optimize HMM parameters based on the read counts on a chromosome
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Inheritance function from \code{\link{nbh}} that receives an object of GRanges class with additional column of read counts (for each strand) and call \code{\link{nbh.integer}} to derive the most probable sequence of hidden states
}
\usage{
\method{nbh}{GRanges}(x, K, ...)
}



%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{binGR} here~~
	GRanges with 'values' slot used for bin counts in 1D vector of integers.
}
  \item{K}{
%%     ~~Describe \code{K} here~~
	Number of hidden states.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
	Extra arguments passed to nbh.integer for the actual HMM computation.
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
	\item{binGR}{GRanges of bin counts with metadata slot saved for the optimized HMM parameters including alpha, beta for the K negative binomial mixture components and TRANS (the transition probabilities)}
}
\references{
%% ~put references to the literature/web site here ~
Rabiner, L. R. (1989). A tutorial on hidden Markov models and selected applications in speech recognition (Vol. 77, pp. 257-286). Presented at the Proceedings of the IEEE. doi:10.1109/5.18626

Bishop, Christopher. Pattern recognition and machine learning. Number 605-631 in Information Science and Statisitcs. Springer Science, 2006.

Capp\'e, O. (2001). H2M : A set of MATLAB/OCTAVE functions for the EM estimation of mixtures and hidden Markov models. (\url{http://perso.telecom-paristech.fr/cappe/h2m/})
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
\code{\link{mainSeekSingleChrom}, \link{nbh.integer}}
}
\examples{
if(interactive()) ?nbh # see nbh for example of nbh running on GRanges object
}