source("R/analyse-spikeins.R") # main plotting code

analyse <- function(global = FALSE)
{
  # also antisense python code uses full raw read count instead of featurecounts so overestimates
  # NB featurecounts generated count files will need the first line removed 
  # Use e.g. sed -i '' '1d' <filename> to remove first line (on a Mac, elsewhere omit the '')
  # path <- "~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/antisense-260816/"
  # replicates <- c("Col-01_S1", "Col-02_S1", "Col-03_S2", "Vir-1_S2", "Vir-2_S2", "Vir-3_S1")
  # conditions <- c("Col","Col","Col","Vir","Vir","Vir")
  # data <- lapply(replicates, function(r) read_data(paste(path,r,"-sensecounts.txt",sep=""), paste(path,r,"-anticounts.txt", sep="")))
  # 
  # # do global correction first, overwrite with individual if required: this covers unspliced genes
  # # which we don't have an individual correction for
  # # do correction with global per replicate ratios
  # by_rep_ratios <- c(0.0031,0.000854,0.0111,0.000203, 0.000269, 0.00015)
  # data <- lapply(1:length(data), function(i) insert_global_ratio(data[[i]], by_rep_ratios[i]))
  # 
  # result <- data
  # if (!global)
  # {
  #   # do correction with individual per replicate ratios
  #   # load data in form of spliced reads which match intron-exon structure on sense strand
  #   splicepath = "/Users/kmourao/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/spliced_antisense/"
  #   splicecounts <- lapply(replicates, function(r) read_spliced_data(splicepath, r))
  #   by_gene_ratios <- splicecounts
  #   result <- lapply(1:length(data), function(i) insert_indiv_ratio(data[[i]], by_gene_ratios[i]))
  # }
  # 
  ## plot all data + curve fitted to non-polyA data
  plot_histograms_bygene(result[1:3], replicates[1:3], "CorrectedCol-fwd.pdf", "CorrectedCol-rev.pdf", rcurve, overlay=TRUE)
  plot_histograms_bygene(result[4:6], replicates[4:6], "CorrectedVir-fwd.pdf", "CorrectedVir-rev.pdf", rcurve, overlay=TRUE)

  res <- result
  
  # add everything back together (yuk, must be a better way to do this)
  rescol <- data.frame(res[[1]]$Geneid, res[[1]]$adjfwd + res[[2]]$adjfwd + res[[3]]$adjfwd, 
                       res[[1]]$adjrev + res[[2]]$adjrev + res[[3]]$adjrev)
  resvir <- data.frame(res[[1]]$Geneid, res[[4]]$adjfwd + res[[5]]$adjfwd + res[[6]]$adjfwd, 
                       res[[4]]$adjrev + res[[5]]$adjrev + res[[6]]$adjrev)
  names(resvir) <- c("Geneid", "adjfwd", "adjrev")
  names(rescol) <- c("Geneid", "adjfwd", "adjrev")

  return(result)
}


##############################################################################################################
# plot_histogram_bygene: Plot a histogram of list of dataset result, optionally with fitted curve to rcurve
# result: list of datasets to plot histogram for
# replicates: list of replicates
# name: name of file to save plot to, must end in .pdf
# rcurve: dataset to plot fitted curve for (can be different to r!)
plot_histograms_bygene <- function(result, replicates, name1, name2, rcurve=NULL, overlay=FALSE)
{
  #png(file=name1,width=1000,height=800)
  pdf(file=name1,width=5,height=8)
  
  if (!overlay)
  {
    par(mfrow = c(3,2))
  }
  else
  {
    par(mfrow=c(3,1))
  }
  
  f <- lapply(1:length(result), function(i) plot_histogram_pair_bygene(result[[i]], "fwd", replicates[i], overlay))
  
  dev.off()
  
  #png(file=name2,width=1000,height=800)
  pdf(file=name2,width=5,height=8)
  
  if (!overlay)
  {
    par(mfrow = c(3,2))
  }
  else
  {
    par(mfrow=c(3,1))
  }
  
  f <- lapply(1:length(result), function(i) plot_histogram_pair_bygene(result[[i]], "rev", replicates[i], overlay))
    
  #if (!is.null(rcurve))
  #{
  #  fitnormal(rcurve, h, xlims[2])
  #}
  
  dev.off()

}

plot_histogram_pair_bygene <- function(d, dir, rep, overlay)
{
  raw = paste(dir, ".anti", sep="")
  adj = paste("adj", dir, sep="")
  xlims = c(0,6)
  ylims = c(0,1500)
  
  if (!overlay)
  {
    r <- as.numeric(d[[raw]]) 
    h <- hist(log10(r), breaks=50, col="lightblue", xlab="log10(Antisense counts)", 
              ylab="Number of genes", main=paste("Raw antisense counts by gene (",rep, " - ", dir,")",sep=""),
              xlim=xlims, ylim=ylims)
    
    r <- as.numeric(d[[adj]]) 
    h <- hist(log10(r), breaks=50, col="lightblue", xlab="log10(Antisense counts)", 
              ylab="Number of genes", main=paste("Adjusted antisense counts by gene (",rep, " - ", dir,")", sep=""),
              xlim=xlims, ylim=ylims)
  }
  else
  {
    r <- as.numeric(d[[raw]]) 
    h1 <- hist(log10(r), breaks=50, col=rgb(0,0,0,0.3), xlab="log10(Antisense counts)", 
              ylab="Number of genes", main=paste("Antisense counts by gene (",rep, " - ", dir,")",sep=""),
              xlim=xlims, ylim=ylims)
    
    r <- as.numeric(d[[adj]]) 
    h2 <- hist(log10(r), breaks=50, col=rgb(0,0,1,0.5), xlab="log10(Antisense counts)", 
              ylab="Number of genes", main=paste("Antisense counts by gene (",rep, " - ", dir,")", sep=""),
              xlim=xlims, ylim=ylims, add=T)
    
    legend("topright", c("Raw counts", "After correction"), col=c(rgb(0,0,0,0.3), rgb(0,0,1,0.5)), lwd=10)
    
    box()
  }
  return(NULL) # remember to return something for apply calls
}