\name{nbh_gen}
\alias{nbh_gen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Simulate data from a negative binomial HMM.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generate count data and the hidden states based on the user-supplied HMM paramters.
}
\usage{
nbh_gen(TRANS, alpha, beta, Total)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TRANS}{
%%     ~~Describe \code{TRANS} here~~
	Expected a squared matrix of probabilities (\eqn{0 \le p \le 1}) with row and column length equal to that of alpha and beta and row sum and column sum both equal to 1 (within some numerical deviation of 1e-6).
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
	Expected a vector of positive values with length equal to that of beta and the row/column of TRANS.
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
	Expected a vector of positive values with length equal to that of alpha and the row/column of TRANS.
}
  \item{Total}{
%%     ~~Describe \code{Total} here~~
	Total number of data points to generate.
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
	\item{count}{Simulation count data.}
	\item{label}{Hidden states associated with the simulated data.}
}
\references{
%% ~put references to the literature/web site here ~
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
\code{\link{nbh_em}, \link{nbm_chk}, \link{randindx}}
}
\examples{
# Simulate data using user-supplied transition prob, alpha and beta for the NB HMM parameters
TRANS_s <- matrix(c(0.9, 0.1, 0.3, 0.7), nrow=2, byrow=TRUE)
alpha_s <- c(2, 4)
beta_s  <- c(1, 0.25)
Total <- 100
x <- nbh_gen(TRANS_s, alpha_s, beta_s, Total)
}