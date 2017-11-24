\name{nbm_em}
\alias{nbm_em}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Expectation conditional maximization of likelihood for negative binomial mixture model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given an input read count vector of integers, the function optimzes the parameters for the negative binomial mixture model of K components using expectation conditional maximization.
}
\usage{
nbm_em(count, alpha, beta, wght, NBM_NIT_MAX = 250, NBM_TOL = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{count}{
%%     ~~Describe \code{count} here~~
	A vector of integers, conceptaully representing the read counts within bins of chromosome.
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
	Initial values for \eqn{\alpha_k} for all K NB.
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
	Initial values for \eqn{\beta_k} for all K NB.
}
  \item{wght}{
%%     ~~Describe \code{wght} here~~
	Initial values for \eqn{\pi_k} for all K NB.
}
  \item{NBM_NIT_MAX}{
%%     ~~Describe \code{NBM_NIT_MAX} here~~
	Maximum number of EM iterations (Default: 250).
}
  \item{NBM_TOL}{
%%     ~~Describe \code{NBM_TOL} here~~
	Threshold as fraction of increase in likelihood (given the current NBM parameters) comparing with the likelihood from the last iteration. EM for the NBM stops when the improvement is below the threshold (Default: 0.01).	
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Given a \eqn{K}-NBM, the goal is to maximize the likelihood function with respect to the parameters comprising of \eqn{\alpha_k} and \eqn{\beta_k} for the K NB components and the mixing coefficients \eqn{\pi_k}, which are the priors \eqn{p(z=k)}. Because there is no analytical solution for the maximum likelihood (ML) estimators of the above quantities, a modified EM procedures called Expectation Conditional Maximization is employed (Meng and Rubin, 1994).

In E-step, the posterior probability is evaluated using NB density functions with initialized \eqn{\alpha_k}, \eqn{\beta_k}, and \eqn{\pi_k}. In the CM step, \eqn{\pi_k} is evaluated first followed by Newton updates of \eqn{\alpha_k} and \eqn{\beta_k}. EM iteration terminates when the percetnage of increase of log likelihood drop below \code{NBM_TOL}, which is deterministic since EM is guaranteed to converge. For more details, please see the manuscript of RIPSeeker.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list containing:
\item{alpha}{alpha_k for all K components of NB.}
\item{beta}{beta_k for all K components of NB.}
\item{wght}{pi_k for all K components of NB.}
\item{logl}{Log likelihood in each EM iteration.}
\item{postprob}{Posterior probabilities for each observed data point in the last EM iteration.}
}
\references{
%% ~put references to the literature/web site here ~
Bishop, Christopher. Pattern recognition and machine learning. Number 605-631 in Information Science and Statisitcs. Springer Science, 2006.

X. L. Meng, D. B. Rubin, Maximum likelihood estimation via the
ECM algorithm: A general framework, Biometrika, 80(2):267-278 (1993).

J. A. Fessler, A. O. Hero, Space-alternating generalized expectation-maximization algorithm, IEEE Tr. on Signal Processing, 42(10):2664 -2677 (1994).

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
\code{ \link{nbh_init}, \link{nbh}, \link{nbh.GRanges}, \link{nbh_em}}
}
\examples{
# Simulate data
TRANS_s <- matrix(c(0.9, 0.1, 0.3, 0.7), nrow=2, byrow=TRUE)
alpha_s <- c(2, 4)
beta_s  <- c(1, 0.25)
Total <- 1000
x <- nbh_gen(TRANS_s, alpha_s, beta_s, Total);

N <- 2

cnt <- x$count
label <- x$label

Total <- length(cnt)

# dummy initialization
wght0 <- c(0.5,0.5)

alpha0 <- c(1, 20)

beta0 <- c(1, 1)

NIT_MAX <- 50
TOL <- 1e-100

# initialize param with nbm

nbm <- nbm_em(cnt, alpha0, beta0, wght0, NIT_MAX, TOL)

map.accuracy <- length(which(max.col(nbm$postprob) == label))/Total

print(map.accuracy)
}