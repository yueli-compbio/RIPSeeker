# nbh_init  Initalize parameters for nbh_em
# 			Function nbm_em (NB mixture model) is used to find alpha,
# 			beta, and wght (mixprop); wght (1xN) is repeated N times row-wise
# 			to represent the initial TRANS for the subsequent nbh_em training
# 			Use: nbh0 <- nbh_init(count, K)
#			nbh0: list(TRANS, alpha, beta)
nbh_init <- function(count, K, NBM_NIT_MAX=250, NBM_TOL=1e-3)
{
	
	if(missing(count)) {stop("data count is missing")}	
	if(missing(K)) {stop("Number of cluster K is missing")}
	
	message(sprintf("\nStarting NB mixture model (nbm_em) for K=%d clusters:", K), appendLF=TRUE)
		
	alpha <- tail(quantile(count[count>0], probs=seq(0, 1, 1/K)), K)
	
	beta <- rep(1, K)
	
	wght <- rep(0.5, K)
	
	nbm <- nbm_em(count, alpha, beta, wght, NBM_NIT_MAX=NBM_NIT_MAX, NBM_TOL=NBM_TOL)
	
	alpha <- nbm$alpha
	
	beta <- nbm$beta
	
	TRANS <- matrix(rep(nbm$wght, K), ncol=K, byrow=TRUE)
	
	
	# Order the parameters s.t. they are in increasing order of mean values
	# thus, the parameters for hypo- and hyper-state are always 
	# the first and last, respectively
	
	myorder <- order(alpha/beta)
	
	alpha <- alpha[myorder]
	
	beta <- beta[myorder]
	
	TRANS <- TRANS[myorder, myorder]		
	
	
	nbhInit <- list(TRANS=TRANS, alpha=alpha, beta=beta)
		
	
	nbhInit
}
