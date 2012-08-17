plotRatios <- function(filenames_, plotLabels_=NULL,title_,ylim_, ref_, channels=c("i114", "i115", "i116", "i117"), plotLogRatios_=TRUE, pdf_=FALSE, debug_=FALSE, expectedRatios_=NULL)
{  
  # assign labels
  if(is.null(plotLabels_) | length(filenames_) != length(plotLabels_)) {
    plotLabels_ <- filenames_
  }
    
  # create the plot
  if(pdf_) {
    # pdf(sub("tsv", "pdf", filename_))
    pdf(paste(title_, ".pdf"), paper="a4r")
  }
 
  # initialize plot using the title as file name
  par(mfrow=c(1,length(filenames_))) 
  
  # make file list
  for(f in 1:length(filenames_)) {
    filename_ = filenames_[f]  
    plotLabel_ = plotLabels_[f]
    
    #channels <- c("i114", "i115", "i116", "i117")
    if(is.na( match(ref_, channels) ) ) {
      stop("Ref channel is not in channel list")
    }
    
    experiment = read.table(filename_, header=TRUE)
    
    #print(summary(itraq4plex))
    
    # remove zero entries for now
    clear <- experiment
    for(c in 1:length(channels)) {
      clear <- clear[clear[[channels[c]]] != 0 ,]
    }
    
    ratios <- c()
    logratios <- c()
    
    for(c in 1:length(channels)) {
      if(channels[c] != ref_) {
        C <- channels[c]
        
        ratioColumn <- paste("R", C, "/", ref_) 
        logRatioColumn <- paste("LR", C, "/", ref_)
        # compute ratios
        clear[[ratioColumn]] <- clear[[C]]/clear[[ref_]]
        # compute logratios
        clear[[logRatioColumn]] <- log(clear[[C]]) - log(clear[[ref_]])
        ratios <- cbind(ratios, ratioColumn)
        logratios <- cbind(logratios, logRatioColumn)
      }
    }
    
    if(plotLogRatios_) {
      plotset <- subset(clear, select=logratios) 
    }
    else {
      plotset <- subset(clear, select=ratios)
    }
    
    if(debug_) {
      print(logratios)
      print(summary(plotset))
      print(summary(clear))
    }
        
    boxplot(plotset,ylim=ylim_)
    
    title(main=title_, sub=plotLabel_)  
    
    if(!is.null(expectedRatios_)) {
      for(r in 1:length(expectedRatios_)) {
        abline(h = log(expectedRatios_[r]), lty=4)
      }
    }
  }
  
  if(pdf_) {
    dev.off()
  }
  
}