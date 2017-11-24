# nbh.integer: the main HMM function that runs nbh_init, nbh_em, nbh_vit for input object of class "integer"
nbh.integer <- function(x, K, NBM_NIT_MAX=250, NBM_TOL=1e-2, 
		NBH_NIT_MAX=250, NBH_TOL=1e-3, runViterbi = FALSE, ...)
{	
	count <- x
	
	message(sprintf("\n***1. Initializing negative binomial HMM (nbh) with %d states:\n", K))
	
	nbhInit <- nbh_init(count, K, NBM_NIT_MAX, NBM_TOL)
	
	message(sprintf("\n***2. Traininig nbh with forward-backward algorithm:\n"))

	nbhFB <- nbh_em(count, nbhInit$TRANS, nbhInit$alpha, nbhInit$beta, NBH_NIT_MAX, NBH_TOL)
	
	if(runViterbi) {
	
		message(sprintf("\n***3. Deriving maximum-likelihood hidden state sequence with Viterbi algorithm:\n"))
		
		nbhVit <- nbh_vit(count, nbhFB$TRANS, nbhFB$alpha, nbhFB$beta)
		
		message("done!")
		
		list(	initTRANS=nbhInit$TRANS, initAlpha=nbhInit$alpha, initBeta=nbhInit$beta,
				postprob=nbhFB$postprob, viterbi_state=nbhVit$class, 
				TRANS=nbhFB$TRANS, alpha=nbhFB$alpha, beta=nbhFB$beta	)
	} else {
		
		list(	initTRANS=nbhInit$TRANS, initAlpha=nbhInit$alpha, initBeta=nbhInit$beta,
				postprob=nbhFB$postprob, 
				TRANS=nbhFB$TRANS, alpha=nbhFB$alpha, beta=nbhFB$beta	)	
		
	}
}
