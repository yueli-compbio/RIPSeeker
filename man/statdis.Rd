\name{statdis}
\alias{statdis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Returns the stationary distribution of a Markov chain.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Given a transition matrix \eqn{A}, returns the stationary distribution of a Markov chain by computing the eigen vectors of \eqn{A}.
}
\usage{
statdis(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
	Transition probability matrix, a squared matrix of probabilities (\eqn{0 \le p \le 1}) with row and column length equal to that of alpha and beta and row sum and column sum both equal to 1 (within some numerical deviation of 1e-6).
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
\item{w}{Stationary weights for the distributions of K components based on the transition probability matrix.}
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
\code{\link{nbh_em}}
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

# Plot the marginal distribution (in the stationnary regime)
# Compute negative binomial distributions for all model states
t <- 0:max(count)

tmp <- nbh_em(t, nbh$TRANS, nbh$alpha, nbh$beta, 1)

dens <- tmp$dens

w <- statdis(nbh$TRANS)

# Plot estimate of marginal probabilities
marprob <- apply(t(dens) * (t(w) \%*\% matrix(1, ncol=length(t))), 2, sum)

plot(t, marprob, pch=8, col="blue", main="Estimated marginal distribution")

# Plot empirical estimated probabilities
dhist <- matrix(0, ncol=length(t))

for(i in t){
	dhist[1+i] <- sum(count == i)/Total	
}

points(t, dhist, pch=3, col="red")

}