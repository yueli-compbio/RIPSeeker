\name{randindx}
\alias{randindx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generates random indexes with a specified probability distribution
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns an array of T indexes distributed as specified by p (which should be a normalized probability vector).
}
\usage{
randindx(p, Total, NO_CHK)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
%%     ~~Describe \code{p} here~~
	A row vector of normalized probabilities that dictate the transition probability from the current state to the next state. For example, \eqn{p = [0.2, 0.8]} indicates that the current state transitoins to state 1 at 0.2 and 2 at 0.8. The current state itself can either be the state 1 or 2.
}
  \item{Total}{
%%     ~~Describe \code{Total} here~~
	Total number of states needed to be generated using the input transition vector.
}
  \item{NO_CHK}{
%%     ~~Describe \code{NO_CHK} here~~
	Check whether the first argument is a valid row vector of normalized probabilities.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
	The function is used by \code{nbh_gen} to generate random data point based on the user-supplied transition probability matrix.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
	\item{I}{Index/Indices or state(s) sampled following the transition.}probability.
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
\code{\link{nbh_gen}}
}
\examples{
# Total contains the length of data to simulate
Total <- 100

# number of states
N <- 2

# transition probabilities between states
TRANS <- matrix(c(0.9, 0.1, 0.3, 0.7), nrow=2, byrow=TRUE)

label <- matrix(0, Total, 1)

# Simulate initial state
label[1] <- randindx(matrix(1,ncol=N)/N, 1, 1)

# Use Markov property for the following time index
for(t in 2:Total) {
	
	label[t] <- randindx(TRANS[label[t-1],], 1, 1)
}	

plot(label)

}