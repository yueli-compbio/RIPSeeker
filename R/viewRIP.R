# Function Name:   viewRIP
# Description: 		x
# Input: 			y
# Output:			z
#
# Author: Yue Li
###############################################################################

viewRIP <- function(seekedRIP, alignGR, alignGRCTL, binGR=seekedRIP, scoreType = "eFDR",
                        cutoffLine = 1e-3, displayALLChr = FALSE, ...)
{	
  
  session <- browserSession()
  
  # if not logOddScore, then must be probability score
  # so transform prob to log10 scale to ease visualization
  if(scoreType != "logOddScore") {
    
    values(seekedRIP)[,scoreType] <- -log10(values(seekedRIP)[,scoreType])
    
    # transform the cutoff line also
    if(cutoffLine < 1) cutoffLine <- -log10(cutoffLine)
  }
  
  
  scoreTrack <- seekedRIP[,scoreType]
  
  readcnttrack <- binGR[,"count"]
  
  names(values(scoreTrack)) <- "score"
  
  names(values(readcnttrack)) <- "score"
  
  
  # load tracks
  #	track(session, scoreType, yLineOnOff = TRUE, priority=3,
  #			yLineMark = cutoffLine, 
  #			color=as.integer(col2rgb("purple"))) <- RangedData(scoreTrack)
  #	
  chrnames <- as.character(runValue(seqnames(seekedRIP)))
  
  i <- 1
  
  if(!missing(alignGR)) {
    
    alignGR <- as(alignGR, "GRanges")
    
    if(!displayALLChr) {
      
      alignGR <- alignGR[as.character(seqnames(alignGR)) == chrnames]		
    }		
    
    track(session, "aligned RIP reads", priority=i) <- RangedData(alignGR)
    
  }
  
  if(!missing(alignGRCTL)) {
    
    i <- i + 1
    
    alignGRCTL <- as(alignGRCTL, "GRanges")
    
    if(!displayALLChr) {
      
      alignGRCTL <- alignGRCTL[as.character(seqnames(alignGRCTL)) == chrnames]		
    }		
    
    track(session, "aligned CTL reads", priority=i) <- RangedData(alignGRCTL)
    
  }
  
  i <- i + 1
  
  track(session, "read counts", yLineOnOff = TRUE, 
        color=as.integer(col2rgb("dark green")), priority=i) <- RangedData(readcnttrack)
  
  i <- i + 1
  
  track(session, scoreType, yLineOnOff = TRUE, priority=i,
        yLineMark = cutoffLine, 
        color=as.integer(col2rgb("purple"))) <- RangedData(scoreTrack)
  
  
  
  trackNames(session) ## list the track names
  
  # launch browser view
  # 	browserView(session, range(seekedRIP))
  displayRange <- GRangesForUCSCGenome(genome(seekedRIP)[1],
                                       as.character(seqnames(seekedRIP)[1]),
                                       IRanges(start(range(seekedRIP)),
                                               end=end(range(seekedRIP))))
  
  browserView(session, displayRange)
}
