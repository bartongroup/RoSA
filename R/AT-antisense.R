source("analyse-spikeins.R") # main plotting code

analyse <- function(global = FALSE)
{
  # also antisense python code uses full raw read count instead of featurecounts so overestimates
  # NB featurecounts generated count files will need the first line removed 
  # Use e.g. sed -i '' '1d' <filename> to remove first line (on a Mac, elsewhere omit the '')
  path <- "~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/antisense-260816/"
  replicates <- c("Col-01_S1", "Col-02_S1", "Col-03_S2", "Vir-1_S2", "Vir-2_S2", "Vir-3_S1")
  conditions <- c("Col","Col","Col","Vir","Vir","Vir")
  data <- lapply(replicates, function(r) read_data(paste(path,r,"-sensecounts.txt",sep=""), paste(path,r,"-anticounts.txt", sep="")))
  
  # do global correction first, overwrite with individual if required: this covers unspliced genes
  # which we don't have an individual correction for
  # do correction with global per replicate ratios
  by_rep_ratios <- c(0.005,0.0008,0.0112,0.0011, 0.0008, 0.0008)
  data <- lapply(1:length(data), function(i) insert_global_ratio(data[[i]], by_rep_ratios[i]))
  
  if (!global)
  {
    # do correction with individual per replicate ratios
    # load data in form of spliced reads which match intron-exon structure on sense strand
    splicepath = "/Users/kmourao/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/spliced_antisense/"
    splicecounts <- lapply(replicates, function(r) read_spliced_data(splicepath, r))
    by_gene_ratios <- splicecounts
    result <- lapply(1:length(data), function(i) insert_indiv_ratio(data[[i]], by_gene_ratios[i]))
  }
  
  ## plot all data + curve fitted to non-polyA data
  plot_histograms_bygene(result[1:3], replicates[1:3], "CorrectedCol-fwd.pdf", "CorrectedCol-rev.pdf", rcurve)
  plot_histograms_bygene(result[4:6], replicates[4:6], "CorrectedVir-fwd.pdf", "CorrectedVir-rev.pdf", rcurve)

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
plot_histograms_bygene <- function(result, replicates, name1, name2, rcurve=NULL)
{
  #png(file=name,width=1000,height=800)
  pdf(file=name1,width=10,height=8)
  
  par(mfrow = c(3,2))
  f <- lapply(1:length(result), function(i) plot_histogram_pair_bygene(result[[i]], "fwd", replicates[i]))
  
  dev.off()
  
  pdf(file=name2,width=10,height=8)
  
  par(mfrow = c(3,2))
  f <- lapply(1:length(result), function(i) plot_histogram_pair_bygene(result[[i]], "rev", replicates[i]))
    
  #if (!is.null(rcurve))
  #{
  #  fitnormal(rcurve, h, xlims[2])
  #}
  
  dev.off()
}

plot_histogram_pair_bygene <- function(d, dir, rep)
{
  raw = paste(dir, ".sense", sep="")
  adj = paste("adj", dir, sep="")
  xlims = c(0,7)
  ylims = c(0,2500)
  
  r <- as.numeric(d[[raw]]) 
  h <- hist(log10(r), breaks=50, col="lightblue", xlab="log10(Antisense counts)", 
            ylab="Number of genes", main=paste("Raw antisense counts by gene (",rep, " - ", dir,")",sep=""),
            xlim=xlims, ylim=ylims)
  
  xlims = c(0,6)
  ylims = c(0,1500)
  r <- as.numeric(d[[adj]]) 
  h <- hist(log10(r), breaks=50, col="lightblue", xlab="log10(Antisense counts)", 
            ylab="Number of genes", main=paste("Adjusted antisense counts by gene (",rep, " - ", dir,")", sep=""),
            xlim=xlims, ylim=ylims)
  return(NULL) # remember to return something for apply calls
}


############################################################################################################
# Functions to read in sense and antisense counts for full genes, or just spliced reads

# Read sense and antisense counts data
read_data <- function(sensename, antiname)
{
  sense <- read.csv(sensename, sep="\t")
  anti <- read.csv(antiname, sep="\t")
  names(sense)[c(7,8)] = c("fwd.sense", "rev.sense")
  names(anti)[c(7,8)] = c("fwd.anti", "rev.anti")
  data <- merge(sense,anti,by="Geneid")
  data <- data[c("Geneid","fwd.sense","rev.sense","fwd.anti","rev.anti","Length.x")]
  return(data)
}

# Read sense and antisense spliced counts data
read_spliced_data <- function(path, repname)
{
  data <- read.csv(paste(path,repname,"-antisplicedoutput.csv",sep=""))
  names(data)[1] <- "Geneid"
  return(data)
}

############################################################################################################
# Apply antisense correction rastio on a per gene basis
insert_indiv_ratio <- function(d, ratio)
{
  d <- merge(x = d, y = ratio, by = "Geneid", all.x = TRUE)
  d <- apply_correction(d, "ratio", "ratio")
  d$newratiof[!is.na(d$adjfwd) & d$adjfwd>0] <- d$adjfwd[!is.na(d$adjfwd) & d$adjfwd>0]/d$fwd.sense[!is.na(d$adjfwd) & d$adjfwd>0]
  d$newratior[!is.na(d$adjrev) & d$adjrev>0] <- d$adjrev[!is.na(d$adjrev) & d$adjrev>0]/d$rev.sense[!is.na(d$adjrev) & d$adjrev>0]
  
  # note that counts coming in from spliced counts will potentially be quite different from full gene counts!
  
  return(d)
}

# Apply antisense correction ratio globally across all genes
insert_global_ratio <- function(d, ratio)
{
  d$correction <- rep(ratio, length(d$Geneid))
  d <- apply_correction(d, "correction", "correction")
  return(d)
}

# Apply correction to both strands
apply_correction <- function(d, fwdcorrection, revcorrection)
{
  d$adjfwd[!is.na(d[fwdcorrection])] <- pmax(0, round(d$rev.anti[!is.na(d[fwdcorrection])] - 
                                       d$fwd.sense[!is.na(d[fwdcorrection])] * d[fwdcorrection][!is.na(d[fwdcorrection])]))
  d$adjrev[!is.na(d[fwdcorrection])] <- pmax(0, round(d$fwd.anti[!is.na(d[revcorrection])] - 
                              d$rev.sense[!is.na(d[revcorrection])] * d[revcorrection][!is.na(d[revcorrection])]))
  return(d)
}