\name{nbm_chk}
\alias{nbm_chk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Check the parameters of the negative binomial mixture model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function verifies the numerical range and dimension of the NBM paramters alpha, beta, and wght and returns the number of components. It is used in \code{\link{nbm_em}} before running EM.
}
\usage{
nbm_chk(alpha, beta, wght)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
	Expected a vector of positive values with length equal to that of beta and wght.
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
	Expected a vector of positive values with length equal to that of alpha and wght.
}
  \item{wght}{
%%     ~~Describe \code{wght} here~~
	Expected a vector of probabilities (\eqn{0 \le p \le 1}) with length equal to that of alpha and beta and summed to 1 (within some numerical deviation of 1e-6).
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
Christopher Bishop. Pattern recognition and machine learning. Number 605-631 in Information Science and Statisitcs. Springer Science, 2006.

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
\code{\link{nbm_em}, \link{nbh_chk}}
}
\examples{
# two mixing components
wght <- c(0.5,0.5)

alpha <- c(1, 20)

beta <- c(1, 1)

nbm_chk(alpha, beta, wght)

}