\name{nbh_init}
\alias{nbh_init}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Initialize negative binomial HMM parameters using negative binomial mixture model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function finds a sensible set of initial NB HMM parameters by fitting a NB mixture model of K components using the read count data.
}
\usage{
nbh_init(count, K, NBM_NIT_MAX = 250, NBM_TOL = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{count}{
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Because the EM algorithm in HMM tends to fall into local optimal with poor initialization, NB mixture model with K mixture components (K-NBM) is first applied to the data to obtain a reasonable estimate for the HMM parameters. Given the read count vector, the function applied the lower level function \code{\link{nbm_em}} (NB mixture model) to find alpha, beta, and mixing proportion of the K NB mixture components. Alpha and beta are the parameters of the NB mixture components initialized as the last K quantiles of the nonzero read counts and 1, respectively. The mixing proportions or component weights (wght) of the NB distributions are first initialized as uniform and after EM optimization are used to form a symmetrical transition probability matrix such that probability of state 1 transitioning to state 2 is equal to the probability of state 2 transitioning to state 1. Such matrix is used as the initial transition probability for the HMM model tranining (See \code{\link{nbh_em}}).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list containing:
\item{alpha}{Alpha paramter of the K NB components optimized using \code{\link{nbm_em}}}
\item{beta}{Beta paramter of the K NB components optimized using \code{\link{nbm_em}}}
\item{TRANS}{Transition probability intialized as a symmetrical matrix of mixing proportion of the K NB components optimized using \code{\link{nbm_em}}.}
}
\references{
%% ~put references to the literature/web site here ~
Rabiner, L. R. (1989). A tutorial on hidden Markov models and selected applications in speech recognition (Vol. 77, pp. 257-286). Presented at the Proceedings of the IEEE. doi:10.1109/5.18626

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
\code{ \link{nbm_em}, \link{nbh}, \link{nbh.GRanges}, \link{nbh_em}}
}
\examples{
# Simulate data
Total_train <- 1000

Total_test <- 200

TRANS_s <- matrix(c(0.9, 0.1, 0.5, 0.5), nrow=2, byrow=TRUE)
alpha_s <- c(2, 2)
beta_s  <- c(1, 0.25)

train <- nbh_gen(TRANS_s, alpha_s, beta_s, Total_train)

test <- nbh_gen(TRANS_s, alpha_s, beta_s, Total_test)

nbhInit <- nbh_init(train$count, ncol(TRANS_s))

count <- train$count
label <- train$label

# NBH initialization
nbhInit <- nbh_init(count, ncol(TRANS_s))

TRANS0 <- nbhInit$TRANS
alpha0 <- nbhInit$alpha
beta0 <- nbhInit$beta

# NBH EM
nbh <- nbh_em(count, TRANS0, alpha0, beta0)

map.accuracy <- length(which(max.col(nbh$postprob) == label))/Total_train

vit <- nbh_vit(count, nbh$TRANS, nbh$alpha, nbh$beta)

vit.accuracy <- length(which(vit$class == label))/Total_train

vit_test <- nbh_vit(test$count, nbh$TRANS, nbh$alpha, nbh$beta)

vit_test.accuracy <- length(which(vit_test$class == test$label))/Total_test

nbh_wt_KMLinit <- list(mapAccuracy_train=map.accuracy, vitAccuracy_train=vit.accuracy,
		vitLogl_train=vit$logl, vitAccuracy_test=vit_test.accuracy, 
		vitLogl_test=vit_test$logl)	
}