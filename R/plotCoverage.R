# Function Name: 	plotCoverage
# Description: 		plot genomic coverage across the entire GRanges
# Input: 			GenomicRanges object with bin counts
# Output:			plot
#
# Author: Yue Li
###############################################################################

plotCoverage <- function(x, plotLegend=FALSE, legend.cex=1, ...)
{	
	stopifnot(!missing(x))
	
	plot(start(x), values(x)$count, type = "s", col = "blue", ...)
    
	if(plotLegend) {
		
		chrname <- as.character(runValue(seqnames(x)))
	
		chrlen <- seqlengths(x)[chrname]
	
		# hack legend to ignore extra arguments
		formals(legend) <- c(formals(legend), alist(... = ))
		
		legend("topleft", legend=sprintf("%s:1-%d", chrname, chrlen), cex=legend.cex, ...)
	}
}

plotStrandedCoverage <- function(gr, binSize=1000, plotLegend=FALSE, ylim, ...)
{	
	stopifnot(!missing(gr))
	
	gr.posv <- gr[strand(gr) == '+', ]
	
	gr.negv <- gr[strand(gr) == '-', ]
			
	xpos <- binCount(gr.posv, binSize)
		
	xneg <- binCount(gr.negv, binSize)
		
	posCount <- values(xpos)$count
	
	negCount <- values(xneg)$count
	
    if(missing(ylim)) ylim <- min(max(posCount), max(negCount)) * c(-1, 1)
    
    plotCoverage(x=xpos, ylim = ylim, plotLegend=plotLegend, ...)
    
    lines(start(xneg), -negCount, type = "s", col = "red")
              
    abline(h = 0, col = "dimgray")
}
