# nbh_gen	Simulates data from a negative binomial HMM.
# 			Use: [count,label] = nbh_gen(TRANS,alpha,beta,T|label) where
# 			count is an array of length Total which contains the simulated
# 			data and label contains the corresponding simulated state
# 			sequence. If the last argument is of length one, it is
# 			taken as the number of observations, otherwise it is
# 			considered as a specified state sequence.
# 			Requires Matlab's Statistics toolbox or GNU Octave.
nbh_gen <- function(TRANS, alpha, beta, Total)
{	
	# Inputs arguments
	if(missing(TRANS)){stop("TRANS is missing")}
	if(missing(alpha)){stop("alpha is missing")}
	if(missing(beta)){stop("beta is missing")}
	if(missing(Total)){stop("Total is missing")}
	
	N <- nbh_chk(TRANS, alpha, beta)
	
	# Make sure that the output will be a column vector
	alpha <- matrix(alpha, N, 1)
	beta <- matrix(beta, N, 1)
		
	if(length(Total)==1) {
		# Total contains the length of data to simulate
		# First simulate labels
		label <- matrix(0, Total, 1)
		
		# Simulate initial state
		label[1] <- randindx(matrix(1,ncol=N)/N, 1, 1)
		
		# Use Markov property for the following time index
		for(t in 2:Total) {
			
			label[t] <- randindx(TRANS[label[t-1],], 1, 1)
		}
	} else {
		# Total directly contains a sequence of labels
		label <- matrix(Total, length(Total), 1)		
	}
	
	rates <- rep(0, Total)
	
	# First draw the rates, then the Poisson data
	for(i in 1:Total) {
		rates[i] <- rgamma(1, shape=alpha[label[i]], rate=beta[label[i]])
	}
	count <- rpois(Total, rates)
	
	# return count and label
	list(count=count, label=label)
	
}