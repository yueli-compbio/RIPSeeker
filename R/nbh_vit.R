# nbh_vit   A posteriori sequence estimation for negative binomial HMM
# 			using dynamic programming (also known as Viterbi algorithm).
# 			Use: [class,logl] = nbh_vit(count,TRANS,alpha,beta).

nbh_vit <- function(count, TRANS, alpha, beta)
{
	
	# Inputs arguments
	if(missing(count)){stop("TRANS is missing")}
	if(missing(TRANS)){stop("TRANS is missing")}
	if(missing(alpha)){stop("alpha is missing")}
	if(missing(beta)){stop("beta is missing")}
	
	# Data length
	Total <- length(count)	
	if(any(count < 0) || any(count != round(count))){
		stop("Data does not contain positive integers.")
	}

	count <- matrix(count, nrow=Total, 1)
	# Number of mixture components
	N <- nbh_chk(TRANS, alpha, beta)
	alpha <- matrix(alpha, 1, N)
	beta <- matrix(beta, 1, N)
	
	# 0: Compute density values on a log scale
	# Compute log(count!), the second solution is usually much faster
	# except if max(count) is very large
	cm <- max(count)
	if(cm > 50000){
		dnorm <- as.matrix(lgamma(count + 1))
	} else {
		tmp <- cumsum(rbind(0, log(as.matrix(1:max(count)))))
		dnorm <- as.matrix(tmp[count+1])
	}
	
	logdens <- ( matrix(1, nrow=Total) %*% (alpha * log(beta/(1+beta)) - lgamma(alpha)) -
			count %*% log(1+beta) + lgamma(count %*% matrix(1, ncol=N) + matrix(1, nrow=Total) %*% alpha)
				- dnorm %*% matrix(1, ncol=N) )
	
	
	# Dynamic programming recursion
	# HMM transition probabilities in log
	TRANS <- log(TRANS + .Machine$double.xmin)
	
	# Partial loglikelihood array and bactracking array
	PLOGL <- matrix(0, Total, N)
	BCKTR <- matrix(0, Total-1, N)
	
	# Use uniform a priori probability for the initial state
	PLOGL[1,] <- logdens[1,] - log(N)
	
	for(t in 2:Total) {		
		# get z_t-1 that max logl for each state of z_t		
		tmp <- (t(PLOGL[t-1,,drop=F]) %*% matrix(1, ncol=N)) + TRANS		
		PLOGL[t,] <- apply(tmp, 2, max)
		
		# use transpose to find max position for each column
		BCKTR[t-1,] <- max.col(t(tmp))
				
		PLOGL[t,] <- PLOGL[t,] + logdens[t,]
	}
	
		
	# Backtracking
	class <- matrix(0,nrow=Total)
		
	class[Total] <- which.max(PLOGL[Total, ])
	
	for(t in (Total-1):1) {
		
		class[t] <- BCKTR[t, class[t+1]]		
	}
	
	list(class=class, logl=max(PLOGL[Total, ]))
}