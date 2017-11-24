# Function Name: 	exportGRanges
# Description: 		x
# Input: 			y
# Output:			z
#
# Author: Yue Li
###############################################################################

exportGRanges <- function(gRanges, outfile, exportFormat)
{	
	stopifnot(!missing(gRanges))

	if(exportFormat == "txt") {
		
		write.table(as.data.frame(gRanges), file=outfile,
				sep="\t", quote=F, row.names=FALSE)
	} else {
		
		# use export function that export variaous formats (?export)
		export(gRanges, con=outfile, format=exportFormat)
	}
}