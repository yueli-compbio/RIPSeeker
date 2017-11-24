\name{nbh.integer}
\alias{nbh.integer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
HMM posterior decoding and NB parameter optimization
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Inherithance function from \code{\link{nbh}} that receives a vector of integers and compute optimal HMM parameters via EM algorithm.
}
\usage{
\method{nbh}{integer}(x, K, NBM_NIT_MAX = 250, 
	NBM_TOL = 0.01, NBH_NIT_MAX = 250, 
	NBH_TOL = 0.001, runViterbi = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{count} here~~
	A vector of integers, conceptaully representing the read counts within bins of chromosome.
}
  \item{K}{
%%     ~~Describe \code{K} here~~
	Number of hidden states.
}
  \item{NBM_NIT_MAX}{
%%     ~~Describe \code{NBM_NIT_MAX} here~~
	Maximum number of EM iterations (Default: 250) for the negative binomial mixture model (NBM) intialization step (See \code{\link{nbm_em}}).
}
  \item{NBM_TOL}{
%%     ~~Describe \code{NBM_TOL} here~~
	Threshold as fraction of increase in likelihood (given the current NBM parameters) comparing with the likelihood from the last iteration. EM for the NBM stops when the improvement is below the threshold (Default: 0.01).
}
  \item{NBH_NIT_MAX}{
%%     ~~Describe \code{NBH_NIT_MAX} here~~
	Maximum number of EM iterations (Default: 250) for the negative binomial hidden Markov model (NBH).
}
  \item{NBH_TOL}{
%%     ~~Describe \code{NBH_TOL} here~~
	Threshold as fraction of increase in likelihood (given the current NBH parameters) comparing with the likelihood from the last iteration. EM for the NBH stops when the improvement is below the threshold (Default: 0.001).
}
  \item{runViterbi}{
%%     ~~Describe \code{runViterbi} here~~
	Binary indicator. If TRUE, Viterbi algorithm will be applied to derive the maximum likelihood hidden state sequence using the optimized HMM paramters obtained from the EM (See \code{\link{nbh_em}}).
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
	Extra arguments are ignored.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
	The function consists of three major steps: (1) negarive binomail mixture model used to initialized HMM parameters; (2) optimization of HMM paramters using EM algorithm; (3) Viterbi maximum-liklihood estimation of hidden state sequence. Step (1) involves optimization of NBM parameters assuming the data points are independently sampled from a mixture of K NB distributions (See \code{\link{nbh_init}}). Given the optimized paramters for K-NBM, step (2) drops the independence assumption by introducing the transition probibility between hidden variables, which is initlaized as the mixing proportions of NBM (See \code{\link{nbh_init}}). Given the optimized HMM paramters, step (3) derives the maximum liklihood hidden state sequence using Viterbi algorithm. Step (3) is run only when runViterbi is TRUE.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list containing:
	\item{initAlpha}{Initialized alpha of NBM from \code{\link{nbh_init}}.}
	\item{initBeta}{Initialized beta of NBM from \code{\link{nbh_init}}.}
	\item{initTRANS}{Initialized mixing proportion of NBM from \code{\link{nbh_init}}.}
	\item{postprob}{Posteriors of the K hidden states for each observed count derived from \code{\link{nbh_em}} (e.g., posteriors of background and enriched state in a two-state HMM).}
	\item{alpha}{Optimized alpha of the NB mixture components in the HMM using \code{\link{nbh_em}}.}
	\item{TRANS}{Optimized transition probability of the HMM using \code{\link{nbh_em}}.}
	\item{viterbi_state}{Sequence of discrete values representing the hidden states derived from the maxmium likelihood estimation using Viterbi algorithm (See \code{\link{nbh_vit}}).}
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
\code{\link{mainSeekSingleChrom}, \link{nbh}, \link{nbh.GRanges}}
}
\examples{
if(interactive()) ?nbh # see nbh for example of nbh running on integer object
}