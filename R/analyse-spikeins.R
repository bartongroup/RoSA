##############################################################################################################
# Code to create plots of spike-in data, and associated spliced reads data
##############################################################################################################

library("LSD")

##############################################################################################################
# plot_spikeins: Plot spike-in data sense/antisense and draw linear fit line
# sp: spikein dataset
# name: name of file to output to, must end in .pdf
# label: label for this dataset e.g. Col, Vir
# xmax: max extent of x axis
# ymax: max extent of y axis
# aslog: TRUE if plotting log/log else FALSE
plot_spikeins <- function(sp, name, label, xmax, ymax, aslog=TRUE)
{
  if (!is.null(name))
  {
    pdf(file=name,width=7,height=6)
  }
  
  if (!aslog)
  {
    plot(sp$sense, sp$anti, pch=21, cex=0.8, col='blue', bg='cyan', 
         xlab = 'Sense counts', 
         ylab = 'Antisense counts', xlim = c(0,xmax), ylim = c(0,ymax))
    abline(lm(sp$anti[sp$anti!=0] ~ 0 + sp$sense[sp$anti!=0]), col='blue')
  }
  else
  {
    plot(log10(sp$sense+1), log10(sp$anti+1), pch=21, cex=0.8, col='blue', bg='cyan',
         xlab = 'log(Sense counts+1)', cex.lab=1.5, cex.axis=1.5,
         ylab = 'log(Antisense counts+1)', xlim = c(0,log10(xmax)), ylim = c(0,log10(ymax)))
    use <- is.finite(log10(sp$sense)) & is.finite(log10(sp$anti))
    
    model = lm(log10(sp$anti)[use] ~ log10(sp$sense)[use])
    if (!is.na(model$coefficients[2])) # don't draw a line if gradient is infinite
    {
      abline(model, col='blue')
    }
  }
  title(main=paste(label, " spike-ins: antisense vs sense counts", sep=""), cex.main=1.5)
  
  if (!is.null(name))
  {
    dev.off()
  }
  
  
}

##############################################################################################################
# plot_antisensedata_and_spikeins: plot data and spikeins on same graph, sense vs antisense
# data: full dataset of read counts
# sp: spike in counts
# label: condition or replicate name for graph labels
# name: name of pdf file to output to, must end in .pdf
# aslog: plot as log/log if TRUE, else linear
plot_antisensedata_and_spikeins <- function(data, sp, label, name, aslog=TRUE,
                                            xmax = 2e+07, ymax=5e+04, legx=0, legy=0)
{
  plot_title <- paste(label, ": sense vs antisense counts", sep="")
  
  # symbols
  data_pch <- 21  # symbol for data points - circle
  sp_pch <- 21    # symbol for spike-in points
  
  # #legend
  # legend_labels <- c(paste(label, " spliced reads", sep=""), paste(label, " spike-ins", sep=""), paste(label, " spike-in ratio", sep=""))
  # legend_colours <- c("blue", "black", "black")
  # legend_bg <- c(NA, "black", NA)
  # legend_lty <- c(NA, NA, 1)
  # legend_pch <- c(data_pch, sp_pch, NA)
  
  legend_labels <- c(paste(label, " spike-ins", sep=""), paste(label, " spike-in ratio", sep=""))
  legend_colours <- c("black", "black")
  legend_bg <- c("black", NA)
  legend_lty <- c(NA, 1)
  legend_pch <- c(sp_pch, NA)
  
  if (!aslog & legx==0 & legy==0)
  {
    legx = 0.6 * xmax
    legy = 0.4 * ymax
  }
  else if (legx==0 & legy==0)
  {
    legx = 0.6 * log10(xmax)
    legy = 0.4 * log10(ymax)
  }
  
  pdf(file=name,width=7,height=6)
  if (!aslog)
  {
    plot(data$sensecounts, data$anticounts, pch=data_pch, cex=0.5, col="black", 
         xlab='Sense counts', ylab='Antisense counts', xlim = c(0,xmax), ylim = c(0,ymax))
    par(new=TRUE)
    plot(sp$sense, sp$anti, pch=sp_pch, cex=0.8, col='blue', bg='cyan', axes = FALSE, xlab = '', ylab = '', 
         xlim = c(0,xmax), ylim = c(0,ymax))
    abline(lm(sp$anti ~ 0 + sp$sense), col='blue')
    legend(legx, legy, legend_labels, col = legend_colours, pt.bg = legend_bg, lty = legend_lty, pch = legend_pch)
  }
  else
  {
#    plot(log10(data$sensecounts+1), log10(data$anticounts+1), pch=data_pch, cex=0.3, col="black", 
#         xlab='log(Sense counts+1)', ylab='log(Antisense counts+1)', xlim = c(0,log10(xmax)), ylim = c(0,log10(ymax)))
    heatscatter(log10(data$sensecounts/data$Length.x+1),log10(data$anticounts/data$Length.x+1),colpal="bl2gr2rd",
                cor=FALSE,cex.main=1, xlim=c(0,log10(xmax/1000)), ylim=c(0,1.5),
                xlab='log(Sense counts+1)', ylab='log(Antisense counts+1)', main='')
    
    par(new=TRUE)
    plot(log10(sp$sense/sp$length+1), log10(sp$anti/sp$length+1), pch=sp_pch, cex=1, col='white', bg='black', axes = FALSE, xlab = '', ylab = '', 
         xlim = c(0,log10(xmax/1000)), ylim = c(0,1.5))
    use <- is.finite(log10(sp$anti)) & is.finite(log10(sp$sense))
    abline(lm(log10(sp$anti[use]/sp$length[use]) ~ log10(sp$sense[use]/sp$length[use])), col='black')
    legend(legx, legy, legend_labels, col = legend_colours, pt.bg = legend_bg, lty = legend_lty, pch = legend_pch)
    
    
  }
  title(main=plot_title)
  dev.off()
  
  # plot an MA plot too
  # M = log2(antisense) - log2(sense)
  # A = 1/2 * (log2(antisense) + log2(sense))
  
  Msplice = log(data$anticounts+1) - log(data$sensecounts+1)
  Asplice = 0.5 * (log(data$anticounts+1) + log(data$sensecounts)+1)
  #Asplice = log(data$Length.x)
  Mspikes = log(sp$anti + 1) - log(sp$sense+1)
  Aspikes = 0.5 * (log(sp$anti+1) + log(sp$sense+1))
  #Aspikes = log(sp$length)
  
  plot(Asplice, Msplice, pch=data_pch, cex=0.3, col="black",
       xlim=c(0,14), ylim=c(-12,3))
  heatscatter(Asplice,Msplice,colpal="bl2gr2rd",
              cor=FALSE,cex.main=1, xlim=c(0,14), ylim=c(-12,3))
  par(new=TRUE)
  plot(Aspikes, Mspikes, pch=sp_pch, cex=1, col='black', bg='black',
       xlim=c(0,14), ylim=c(-12,3))
  
  #plot(log10(data$Length.x), log10(data$anticounts+1))
}

##############################################################################################################
# plot_all: Plot 2 datsets sense vs antisense + spike-in ratios
# d1,d2: read counts data
# sp1,sp2: corresponding spikein counts
# label1,label2: labels to use in legend for each dataset
# name: name of pdf file to create, must end in .pdf
# xmax: max extent of x axis
# ymax: max extent of y axis
# aslog: plot as log/log if TRUE, else linear
# legx: x position of legend top left corner
# legy: y position of legend top left corner
plot_all <- function(d1,d2,sp1,sp2,label1,label2,name,xmax,ymax,aslog=TRUE,legx=0,legy=0)
{
  pdf(file=name,width=7,height=5)
  
  # symbols
  data_pch <- 21  # symbol for data points - circle
  
  # legend
  legend_labels <- c(label1, label2, paste(label1," spike-in ratio", sep=""), 
                     paste(label2," spike-in ratio", sep=""))
  legend_colours <- c("red", "black", "red", "black")
  legend_bg <- c(NA, 'black', NA, NA)
  legend_lty <- c(NA, NA, 1, 1)
  legend_pch <- c(data_pch, data_pch, NA, NA)
  
  if (!aslog & legx==0 & legy==0)
  {
    legx = 0.55 * xmax
    legy = 0.4 * ymax
  }
  else if (legx==0 & legy==0)
  {
    legx = 0.55 * log10(xmax)
    legy = 0.4 * log10(ymax)
  }
  
  if (aslog)
  {
    plot(log10(d1$sensecounts+1), log10(d1$anticounts+1), pch=data_pch, cex=0.3, col="black", bg="black", 
         xlab='log(Sense counts+1)', ylab='log(Antisense counts+1)', xlim = c(0,log10(xmax)), ylim = c(0,log10(ymax)))
    par(new=TRUE)
    use <- is.finite(log10(sp1$sense)) & is.finite(log10(sp1$anti))
    abline(lm(log10(sp1$anti)[use] ~ log10(sp1$sense)[use]), col='black')
    par(new=TRUE)
    use <- is.finite(log10(sp2$sense)) & is.finite(log10(sp2$anti))
    abline(lm(log10(sp2$anti)[use] ~ log10(sp2$sense)[use]), col='red')
    par(new=TRUE)
    plot(log10(d2$sensecounts+1), log10(d2$anticounts+1), pch=data_pch, cex=0.3, col="red", axes = FALSE, 
         xlab = '', ylab = '', xlim = c(0,log10(xmax)), ylim = c(0,log10(ymax)))
    
    legend(legx, legy, legend_labels, col = legend_colours, pt.bg = legend_bg, lty = legend_lty, pch = legend_pch)
  }
  else
  {
    plot(d1$sensecounts, d1$anticounts, pch=data_pch, cex=0.5, col="black", bg="black", 
         xlab='Sense counts', ylab='Antisense counts', xlim = c(0,xmax), ylim = c(0,ymax))
    par(new=TRUE)
    abline(lm(sp1$anti ~ 0 + sp1$sense), col='black')
    par(new=TRUE)
    abline(lm(sp2$anti ~ 0 + sp2$sense), col='red')
    par(new=TRUE)
    plot(d2$sensecounts, d2$anticounts, pch=data_pch, cex=0.5, col="red", axes = FALSE, 
         xlab = '', ylab = '', xlim = c(0,60000), ylim = c(0,3e+06))
    
    legend(legx, legy,legend_labels,  col = legend_colours, pt.bg = legend_bg, lty = legend_lty, pch = legend_pch)
  }
  title(main=paste(label1," & ", label2, ": Sense vs antisense counts", sep=""))
  
  dev.off()
}

##############################################################################################################
# plot_spikeins_by_id: Plot spike-ins by id (e.g. sorted by number of sense counts)
# sp1,sp2: spikein counts
# label1,label2: labels to use in legend for each dataset
# aslog: plot as log/log if TRUE, else linear
plot_spikeins_by_id <- function(sp1, sp2, label1, label2, aslog=TRUE)
{
  sorted_sp1 <- sp1[order(sp1$sense),]
  sorted_sp2 <- sp2[order(sp2$sense),]
  
  if (aslog)
  {
    plot_spikeins_by_id_and_strand(log10(sorted_sp1$sense),log10(sorted_sp2$sense),label1,label2,
                                   "Sense", "log(Number of assigned (sense) reads)", log10(1600000),log10(1500000))
    plot_spikeins_by_id_and_strand(log10(sorted_sp1$anti),log10(sorted_sp2$anti),label1,label2,
                                   "Antisense", "log(Number of assigned (antisense) reads)",log10(9000),log10(8000))
  }
  else
  {
    plot_spikeins_by_id_and_strand(sorted_sp1$sense,sorted_sp2$sense,label1,label2,
                                   "Sense", "Number of assigned (sense) reads", 1600000,1500000)
    plot_spikeins_by_id_and_strand(sorted_sp1$anti,sorted_sp2$anti,label1,label2,
                                   "Antisense", "Number of assigned (antisense) reads",9000,8000)
  }
}

# helper function for plot_spikeins_by_id
plot_spikeins_by_id_and_strand <- function(sp1, sp2, label1, label2, strand, ylabel, ymax, legy)
{
  plot(sp1, pch=23, cex=1, col="black", bg="red", xlim=c(0,100), ylim=c(0,ymax),
       xlab="Spike-ins (ordered by #assigned sense reads)", ylab=ylabel)
  par(new=TRUE)
  plot(sp2, col="blue", pch=20, axes=FALSE, xlim=c(0,100), ylim=c(0,ymax), xlab="", ylab="")
  title(main="Sense reads assigned to spike-ins")
  legend(0, legy, c(label1, label2), col=c("blue", "blue"), pch=c(23,20), pt.bg=c("red",NA))
}

##############################################################################################################
# plotrep: plot a replicate's spike-in data as sense v antisense, using parameters provided in d
# d: plot details
#   d$rep: data per replicate, organised as sense and anti columns
#   d$colour: colour for the corresponding replicate
#   d$symbol: symbol for the corresponding replicate
# xlimit: max x value for plot
# ylimit: max y value for plot
# xlabel: label for x-axis
# ylabel: label for y-axis
# plottitle: title for plot
plotrep <- function(d, xlimit, ylimit, xlabel, ylabel, plottitle)
{
  if (d$isfirst)
  {
    # first data is plotted on new graph
    par(new=FALSE)
    
    # first call also sets up x and y limits, axis labels and plot title
    plot(log10(d$rep$sense), log10(d$rep$anti), pch=d$symbol, cex=0.5, col=d$colour,
         xlab = xlabel, ylab = ylabel, xlim = xlimit, ylim = ylimit)
    title(main = plottitle)
    par(new=TRUE)
  }
  else
  {
    # later calls do not reset labels, axes or title
    plot(log10(d$rep$sense), log10(d$rep$anti), pch=d$symbol, cex=0.5, col=d$colour,
         xlab = '', ylab = '', axes=FALSE, xlim = xlimit, ylim = ylimit)
    par(new=TRUE)
  }
  
  # plot linear fit lines
  use <- is.finite(log10(d$rep$sense)) & is.finite(log10(d$rep$anti))
  abline(lm(log10(d$rep$sense)[use] ~ log10(d$rep$anti)[use]), col=d$colour)
  
  # keep this plot for further data
  par(new=TRUE)
}

##############################################################################################################
# plot_spikeins_by_rep: Plot spikeins by replicate
# details: 
#  details$rep = data for each condition, organised in sense and anti columns
#  details$isfirst = TRUE for first condition, FALSE for the rest
#  details$condition = vector of condition names
#  details$colour = vector of colours for each condition
#  details$symbol = vector of symbols for each condition
# pdfname: name of pdf to output plot to
# title: title of plot
# maxx: max x value for plot
# maxy: max y value for plot
# w: width of pdf in inches
# h: height of pdf in inches
# legendx: x position for top left of legend
# legendy: u position for top left of legend
plot_spikeins_by_rep <- function(details, pdfname, title, maxx, maxy, w=7, h=6, legendx=6, legendy=4)
{
  pdf(file=pdfname,width=w,height=h)
  temp <- apply(details, 1, function(d) plotrep(d,c(0,log10(maxx)), c(0,log10(maxy)), 
                                                'log(Antisense counts)', 'log(Sense counts)', title))
  legend(legendx, legendy, details$condition, col=details$colour, pch=details$symbol)
  dev.off()
}

##############################################################################################################
# load_spike_in_rep: load spikein data for 1 replicate, return as list of lanes data
# repname: name of replicate
# path: path to replicate
# numlanes: number of lanes to collate - expect lanes to be in files labelled 
#                                        <path><repname>L00<lane number>_spikein_<fwd/rev>counts.csv
load_spike_in_rep <- function(repname, path, numlanes)
{
  suffixf <- "_spikein_fwdcounts.csv"
  suffixr <- "_spikein_revcounts.csv"
  filef <- paste(path, repname, sep="")
  filer <- paste(path, repname, sep="")
  
  replist <- lapply(1:numlanes, function(l) load_counts_fwdandrev(paste(filef, "L00", l, suffixf, sep=""), 
                                                                  paste(filer, "L00", l, suffixr, sep="")))
  
  # assumes all geneids in same order, so sum all sense/anti counts into two columns by row
  sumreps <- Reduce('+', replist)
  reps <- data.frame(replist[[1]]$Geneid, sumreps$sense, sumreps$anti)
  
  # rename columns
  names(reps) <- c("Geneid", "sense", "anti")
  return(reps)
}

##############################################################################################################
# load_counts_fwdandrev: load spike-in data from fwd and rev files
# file1: forward counts file, stripped of initial lines, tab delimited
# file2: reverse counts file, ditto
load_counts_fwdandrev <- function(file1, file2)
{
  f <- read.csv(file1, sep="\t")
  r <- read.csv(file2, sep="\t")
  
  names(f)[7]<- "sense"
  names(r)[8]<- "anti"
  
  l <- merge(f,r,by="Geneid")
  
  return(l)
}

##############################################################################################################
# load_spike_in_rep_direct: Load spikein rep directly from rep counts (rather than lanes)
# expects files to be in path/repname-<fwd/rev>counts.txt
# repname: name of replicate (used in filename)
# path: path to replicate
load_spike_in_rep_direct <- function(repname, path)
{
  suffixf <- "-fwdcounts.txt"
  suffixr <- "-revcounts.txt"
  filef <- paste(path, repname, sep="")
  filer <- paste(path, repname, sep="")
  
  rep <- load_counts_fwdandrev(paste(filef, suffixf, sep=""), paste(filer, suffixr, sep=""))
  rep <- data.frame(rep$Geneid, rep$sense, rep$anti)
  
  names(rep) <- c("Geneid", "sense", "anti")
  return(rep)
}

##############################################################################################################
# calcratios: Calculate antisense:sense ratios for each dataset in a list of conditions
# paths: list of paths to a set of conditions for an experiment
# conditions: list of conditions for an experiment
# expects counts files to found in path/condition-<fwd/rev>counts.txt
# labels: names of experiments, one for each path
calcratios <- function(paths, conditions, labels)
{
  ratios <- sapply(1:length(paths), function(p) sapply(conditions[[p]], function(c) calcratio(c,paths[p])))
  results <- cbind(unlist(labels),unlist(conditions),unlist(ratios))
  row.names(results) <- NULL
  colnames(results) <- c("experiment", "condition", "ratio")
  
  results <- data.frame(results, stringsAsFactors=FALSE)
  return(results)
}

##############################################################################################################
# calcratio: calculate the antisense:sense ratio for a given condition/replicate/lane
# expects files to be in path/condition-<fwd/rev>counts.txt
# condition: name of condition to load
# path: path to condition files
calcratio <- function(condition, path)
{
  rep <- load_spike_in_rep_direct(condition, path)
  m <- lm(rep$anti[rep$anti!=0] ~ 0 + rep$sense[rep$anti!=0])
  return(m$coefficients)
}

##############################################################################################################
# plot_histogram: Plot a histogram of dataset r, optionally with fitted curve to rcurve
# r: dataset to plot histogram for
# xmax: max extent of histogram x-axis
# ymax: max extent of histogram y-axis
# name: name of file to save plot to, must end in .pdf
# rcurve: dataset to plot fitted curve for (can be different to r!)
plot_histogram <- function(r, xmax, ymax, name=NULL, rcurve=NULL)
{
  if (!is.null(name))
  {
    pdf(file=name,width=7,height=5)
  }
  h <- hist(r, breaks=75, col="lightblue", xlab="Antisense:sense ratio", 
            ylab="Number of replicates", main="Distribution of antisense:sense ratios by replicate",
            xlim=c(0,xmax), ylim=c(0,ymax))
  
  if (!is.null(rcurve))
  {
    fitnormal(rcurve, h, xmax)
  }
  
  if (!is.null(name))
  {
    dev.off()
  }
}

##############################################################################################################
# fitnormal: Fit a normal curve to a histogram
# data: the data used to generate the histogram
# hist: the histogram plot
# xmax: maximum value of x to draw curve to (without this would stop at max data point)
fitnormal <- function(data, hist, xmax)
{
  xfit<-seq(0, xmax, length=1000) 
  yfit<-dnorm(xfit, mean=mean(data), sd=sd(data)) 
  yfit <- yfit*diff(hist$mids[1:2])*length(data) 
  lines(xfit, yfit, col="black", lwd=2)
}



##############################################################################################################
# plot_spikeins: Plot spike-in data sense/antisense and draw linear fit line
# sp: spikein dataset
# name: name of file to output to, must end in .pdf
# label: label for this dataset e.g. Col, Vir
# xmax: max extent of x axis
# ymax: max extent of y axis
# aslog: TRUE if plotting log/log else FALSE
plot_spikeins_talk <- function(sp, name, label, xmax, ymax, aslog=TRUE)
{
  if (!is.null(name))
  {
    pdf(file=name,width=7,height=6)
  }
  
  thecol = 'blue'
  thebg = 'cyan'
  thepch = 21
  if ((label=='Col-01-') | (label=="Vir-1-"))
  {
    thecol='black'
    thebg='black'
    thepch = 17
    
  }
  if ((label=='Col-03-') | (label=="Vir-3-"))
  {
    thecol='red'
    thebg='red'
    thepch = 8
    
  }
  
  
  if (!aslog)
  {
    plot(sp$sense, sp$anti, pch=21, cex=0.8, col='blue', bg='cyan', 
         xlab = 'Sense counts', 
         ylab = 'Antisense counts', xlim = c(0,xmax), ylim = c(0,ymax))
    abline(lm(sp$anti[sp$anti!=0] ~ 0 + sp$sense[sp$anti!=0]), col='blue')
  }
  else
  {
    plot(log10(sp$sense+1), log10(sp$anti+1), pch=thepch, cex=0.8, col=thecol, bg=thebg,
         xlab = 'log(Sense counts+1)', cex.lab=1.5, cex.axis=1.5,
         ylab = 'log(Antisense counts+1)', xlim = c(0,log10(xmax)), ylim = c(0,log10(ymax)))
    
    use <- is.finite(log10(sp$sense)) & is.finite(log10(sp$anti))
    
    model = lm(log10(sp$anti)[use] ~ log10(sp$sense)[use])
    if (!is.na(model$coefficients[2])) # don't draw a line if gradient is infinite
    {
      abline(model, col=thecol)
    }
  }
  #title(main=paste(label, " spike-ins: antisense vs sense counts", sep=""), cex.main=1.5)
  if (label=='Col-03-') 
  {
    legend(x=0.3,y=4, legend=c("Col-01", "Col-02", "Col-03"), col=c('black', 'blue', 'red'), pch=c(17,21,8), bty='n', cex=1.25,
           title="Replicates:", y.intersp=1.5)
    title(main="Col condition: antisense vs sense counts", cex.main=1.5)
  }
  else if (label=="Vir-3-")
  {
    legend(x=0.3,y=4, legend=c("Vir-1", "Vir-2", "Vir-3"), col=c('black', 'blue', 'red'), pch=c(17,21,8), bty='n', cex=1.25,
           title="Replicates:", y.intersp=1.5)
    title(main="Vir condition: antisense vs sense counts", cex.main=1.5)
  }
  
  if (!is.null(name))
  {
    dev.off()
  }
  
  
}