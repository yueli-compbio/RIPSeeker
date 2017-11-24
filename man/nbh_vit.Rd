\name{nbh_vit}
\alias{nbh_vit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Derive maximum likelihood hidden state sequence using Viterbi algorithm
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given read counts and HMM parameters (optimized by \code{\link{nbh_em}}), derive the sequence of hidden states that maximizes the joint likelihood of observed and latent data.
}
\usage{
nbh_vit(count, TRANS, alpha, beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{count}{
%%     ~~Describe \code{count} here~~
	A vector of integers, conceptaully representative of the read counts within bins of chromosome.
}
  \item{TRANS}{
%%     ~~Describe \code{TRANS} here~~
	Optimized transition probability matrix, a squared matrix of probabilities (\eqn{0 \le p \le 1}) with row and column length equal to that of alpha and beta and row sum and column sum both equal to 1 (within some numerical deviation of 1e-6).
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
	Optimized shape parameter of the NB as a vector of positive values with length equal to that of beta and the row/column of TRANS.
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
	Optimized inverse scale parameter of the NB as a vector of positive values with length equal to that of beta and the row/column of TRANS.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
Given a K-state HMM with NB emission (NBH), the goal is to find the latent states corresponding to the observed data that maximize the joint likelihood \eqn{ln p(X, Z) = ln p(x_1, \ldots, x_N, z_1, \ldots, z_N)}. The optimal solution is obtained via Viterbi algorithm, which essentially belongs to the more general framework of Dynamic Programming.

Briefly, starting from the second node of the Markov chain, we select state of the first node that maximizes \eqn{ln p(x_1, x_2, z_2 | z_1)} for \emph{every} state of \eqn{z_2}. Then, we move on to the next node and the next until reaching to the last node. In the end, we make choice for the state of the last node that together leads to the maximum \eqn{ln p(X, Z)}. Finally, we backtrack to find the choices of states in all of the intermeidate nodes to form the final solution.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list containing:
	\item{class}{ML sequence of latent states}
	\item{logl}{Log-likelihood corresponding to the latents states \code{class}}
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
\note{
%%  ~~further notes~~
The function is expected to run after learning the model parameters of HMM using \code{\link{nbh_em}} and (optionally) disambiguating the multihits using \code{\link{nbh_vit}}. However, nothing prevents user from running it with a random set of HMM parameters. Also, note that Viterbi algorithm finds the most probable \emph{sequence of states}, which is not the same as maximizing the posterior probabilities for all the individual latent variables. For instance, a observed data point may be classified as from state 2 in the most probable chain in spite its marginal posterior probability for state 2 is zero.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{ \link{nbh_init}, \link{nbh}, \link{nbh.GRanges}, \link{nbh_em},\link{nbm_em}}
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