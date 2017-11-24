\name{nbh_chk}
\alias{nbh_chk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Check the parameters of the negative binomial HMM
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function verifies the numerical range and dimension of the NBH paramters alpha, beta, and TRANS and returns the number of hidden states. It is used in \code{\link{nbh_em}} before running EM.
}
\usage{
nbh_chk(TRANS, alpha, beta)
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
\item{N}{Number of components or equivalently the length of alpha, beta, or wght.}
}
\references{
%% ~put references to the literature/web site here ~
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
\code{\link{nbh_em}, \link{nbm_chk}}
}
\examples{
# two hidden states
TRANS <- matrix(c(0.9, 0.1, 0.3, 0.7), nrow=2, byrow=TRUE)

alpha <- c(2, 4)

beta  <- c(1, 0.25)

nbh_chk(TRANS, alpha, beta)

}











