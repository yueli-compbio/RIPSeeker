# evalBinSize: evaluate bin size using Shimazaki cost function
# This function should be called only within selectBinSize.R
# Author: Yue
###############################################################################

evalBinSize <- function(binSize, alignGR)
{					
	
	# selecting chromosome with at least a non-zero count
	# esp. when anlayses are separate for individual chromosomes
	ki <- binCount(alignGR, binSize, returnBinCountOnly = TRUE)
		
	# apply Shimazaki cost formula
	k <- mean(ki)
	
	v <- sum((ki - k)^2)/length(ki)
	
	cost <- (2*k - v)/binSize^2
		
	return(cost)	
}
