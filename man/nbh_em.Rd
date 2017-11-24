\name{nbh_em}
\alias{nbh_em}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Expectation conditional maximization of negative binomial HMM parameters using forward-backward algorithm
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given an input read count vector of integers, the function optimizes the parameters for the negative binomial HMM of K hidden states using expectation conditional maximization with forward-backward algorithm to acheive the exact inference.
}
\usage{
nbh_em(count, TRANS, alpha, beta, NBH_NIT_MAX = 250, 
	NBH_TOL = 1e-05, MAXALPHA = 1e+07, MAXBETA = 1e+07)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{count}{
%%     ~~Describe \code{count} here~~
	A vector of integers, conceptaully representative of the read counts within bins of chromosome.
}
  \item{TRANS}{
%%     ~~Describe \code{TRANS} here~~
	Transition probability matrix, a squared matrix of probabilities (\eqn{0 \le p \le 1}) with row and column length equal to that of alpha and beta and row sum and column sum both equal to 1 (within some numerical deviation of 1e-6).
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
	Shape parameter of the NB as a vector of positive values with length equal to that of beta and the row/column of TRANS.
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
	Inverse scale parameter of the NB as a vector of positive values with length equal to that of beta and the row/column of TRANS.
}
  \item{NBH_NIT_MAX}{
%%     ~~Describe \code{NBH_NIT_MAX} here~~
	Maximum number of EM iterations (Default: 250) for the negative binomial hidden Markov model (NBH).
}
  \item{NBH_TOL}{
%%     ~~Describe \code{NBH_TOL} here~~
	Threshold as fraction of increase in likelihood (given the current NBH parameters) comparing with the likelihood from the last iteration. EM for the NBH stops when the improvement is below the threshold (Default: 0.001).
}
  \item{MAXALPHA}{
%%     ~~Describe \code{MAXALPHA} here~~
	The maximum value of alpha in case the update goes beyond the numerical upper limit of the system. Once alpha becomes larger than \code{MAXALPHA}, the EM itaration is prematurely terminated to prevent malfunction.
}
  \item{MAXBETA}{
%%     ~~Describe \code{MAXBETA} here~~
	The maximum value of beta in case the update goes beyond the numerical upper limit of the system. Once beta becomes larger than \code{MAXBETA}, the EM itaration is prematurely terminated to prevent malfunction.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Given a K-state HMM with NB emission (NBH), the goal is to maximize the likelihood function with respect to the parameters comprising of \eqn{\alpha_k} and \eqn{\beta_k} for the K NB components and the transition probabilities \eqn{A_jk} between any state \eqn{j} and \eqn{k}, which are the priors \eqn{p(z=k)}. Because there is no analytical solution for the maximum likelihood (ML) estimators of the above quantities, a modified EM procedures called Expectation Conditional Maximization is employed (Meng and Rubin, 1994).

In E-step, the posterior probability is evaluated by forward-backward algorithm using NB density functions with initialized alpha, beta, and TRANS. In the CM step, \eqn{A_jk} is evaluated first followed by Newton updates of \eqn{\alpha_k} and \eqn{\beta_k}. EM iteration terminates when the percetnage of increase of log likelihood drop below \code{NBH_TOL}, which is deterministic since EM is guaranteed to converge. For more details, please see the manuscript of RIPSeeker.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list containing:
	\item{alpha}{optimized alpha_k for NB at state K}
	\item{beta}{optimized beta_k for NB at state K}
	\item{TRANS}{optimized transition probability matrix}
	\item{logl}{Log likelihood in each EM iteration.}
	\item{postprob}{Posterior probabilities for each observed data point at the last EM iteration.}
	\item{dens}{the negative binomial probabilities computed at the last EM iteration}
}
\references{
%% ~put references to the literature/web site here ~
	Rabiner, L. R. (1989). A tutorial on hidden Markov models and selected applications in speech recognition (Vol. 77, pp. 257-286). Presented at the Proceedings of the IEEE. doi:10.1109/5.18626


Christopher Bishop. Pattern recognition and machine learning. Number 605-631 in Information Science and Statisitcs. Springer Science, 2006.

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
\code{ \link{nbh_init}, \link{nbh}, \link{nbh.GRanges}, \link{nbh_vit},\link{nbm_em}}
}
\examples{
# Simulate data
TRANS_s <- matrix(c(0.9, 0.1, 0.3, 0.7), nrow=2, byrow=TRUE)
alpha_s <- c(2, 4)
beta_s  <- c(1, 0.25)
Total <- 100

x <- nbh_gen(TRANS_s, alpha_s, beta_s, Total);

count <- x$count
label <- x$label

Total <- length(count)

# dummy initialization
TRANS0 <- matrix(rep(0.5,4), 2)

alpha0 <- c(1, 20)

beta0 <- c(1, 1)

NIT_MAX <- 50
TOL <- 1e-100
nbh <- nbh_em(count, TRANS0, alpha0, beta0, NIT_MAX, TOL)

map.accuracy <- length(which(max.col(nbh$postprob) == label))/Total

vit <- nbh_vit(count, nbh$TRANS, nbh$alpha, nbh$beta)

vit.accuracy <- length(which(vit$class == label))/Total

# Plots
par(mfrow=c(2,2), cex.lab=1.2, cex.main=1.2)

plot(count, col="blue", type="l", main=sprintf("A. Simulated Data (Total = \%i)",Total))

plot(as.numeric(nbh$logl), xlab="EM Iteration", ylab="Log-Likelihood", 
main="B. Log-Likelihood via EM");grid()


# Marginal postprob
plot(nbh$postprob[,2], col="blue", type="l", ylim = c(0,1),
ylab="Marginal Posteriror or True State")
points(label-1, col="red")
title(main = sprintf("C. MAP Prediciton Accuracy = \%.2f\%s", 100 * map.accuracy, "\%"))


# Viterbi states
plot(vit$class - 1, col="dark green", type="l", ylim = c(0,1),
ylab="Viterbi or True State")
points(label-1, col="red")
title(main = sprintf("D. Viterbi Prediciton Accuracy = \%.2f\%s", 100 * vit.accuracy, "\%"))

		
}