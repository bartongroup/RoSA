source("R/analyse-spikeins.R")
library("easyGgplot2", lib.loc="~/Library/R/3.2/library")

#================================================================================
# Constants
GENEID = "Geneid"
SENSE  = "sense"
ANTI   = "anti"
LENGTH = "length"

#================================================================================
#' Analyse spike-in antisense counts e.g. analyse("/Users/kmourao/Documents/ENCODE/anti-analysis")
#' 
#' expects data to be organised as follows
#' each group of data in its own directory
#' fwd and rev counts for each sample to be named as <sample_id><fwd suffix> and <sample_id><rev suffix>
#' start directory is top level directory containing the group directories
#' e.g. start_dir/
#'          Lecuyer/
#'              ENCSBL456ABCfwdcounts.txt
#'              ENCSBL456ABCrevcounts.txt
#'              ...
#'          Wold/
#'              ENCSBL987ABCfwdcounts.txt
#'              ENCSBL987ABCrevcounts.txt
#'              ...
#'              
#' Expects each counts file to have a header row, with gene ids in column genecol, 
#' forward counts in column fwdcol and reverse counts in column revcol  
#' 
#' @param start_dir Top level directory for data
#' @param fwdsuffix String at end of file name to indicate forward counts, default fwdcounts.txt
#' @param revsuffix String at end of file name to indicate reverse counts, default revcounts.txt
#' @param skiplines Number of lines to skip at top of each counts file, default 0
#' @param genecol Column number of gene id column, default 1
#' @param fwdcol Column number of (fwd) sense counts column, default 7
#' @param revcol Column number of (rev) antisense counts column, default 8
#' 
#' @export
analyse <- function(start_dir,                  
                    fwdsuffix="fwdcounts.txt", 
                    revsuffix="revcounts.txt", 
                    skiplines=0,
                    genecol=1,
                    lengthcol=6,
                    fwdcol=7,
                    revcol=8)
{
  groups = list.dirs(recursive=FALSE, full.names=FALSE, path=start_dir)
  ratios = lapply(groups, function(g) analyse_by_group(file.path(start_dir,g), g,
                                                      fwdsuffix, revsuffix, skiplines, genecol, lengthcol, fwdcol, revcol))
  
  allratios = do.call("rbind",ratios)
  
  # remove everything with ratio = -1
  # those are 0 sense count cases
  allratios$ratio = as.numeric(allratios$ratio)
  allratios = allratios[allratios$ratio != -1,]
  
  # boxplot the ratios
  par(mfrow=c(1,1), new=FALSE)
  bp = boxplot(as.numeric(log10(allratios$ratio))~allratios$group,data=allratios, main="Antisense:sense ratios", 
               xlab="Group", ylab="Log10(Antisense:sense ratios)", col=c("blue","purple","lightblue","royalblue"))
  # label outliers in boxplot with replicate id. Left aligned.
  #text(bp$group, allratios[allratios$ratio %in% bp$out,]$ratio, 
  #     labels=allratios[allratios$ratio %in% bp$out,]$rep, pos=4)
  
#  sc = stripchart(as.numeric(allratios$ratio)~allratios$group, vertical=TRUE, method="jitter",
    #              col="brown3", pch=16)
  

  
  plot <- ggplot2.stripchart(data=allratios, xName='group',yName='ratio',
                     groupName='group',
                     groupColors=c("black","black","black","black"),
                     showLegend=FALSE,
                     backgroundColor="white", xtitle="", ytitle="Antisense:sense ratios",
                     mainTitle="Antisense:sense ratios by group",
                     addBoxplot=TRUE, boxplotFill=c("#E69F00", "#56B4E9","#E50EEE", "#56B4E9"))
  print(plot)
  
  logdata = allratios
  logdata$ratio = log10(logdata$ratio)
  plot <- ggplot2.stripchart(data=logdata, xName='group',yName='ratio',
                             groupName='group',
                             groupColors=c("black","black","black","black"),
                             showLegend=FALSE,
                             backgroundColor="white", xtitle="", ytitle=expression('log'[10]*'(Antisense:sense ratios)'),
                             mainTitle="Antisense:sense ratios by group",
                             addBoxplot=TRUE, boxplotFill=c("darkorange2","darkorange2","darkorange2","darkorange2"))
  print(plot)

  rcurve <- as.numeric(allratios$ratio[allratios$ratio < 0.002])
  plot_histogram(as.numeric(allratios$ratio), 0.012, 275, "ENCODEhist.pdf", rcurve) #)
  
#  curve(dnorm(x, mean(as.numeric(allratios$ratio)) / length(allratios$ratio), 
#              sd(as.numeric(allratios$ratio))), add=TRUE, col="darkblue", lwd=2)
  
  bp = boxplot(as.numeric(allratios$ratio),data=allratios, main="Antisense:sense ratios", 
               xlab="All Groups", ylab="Antisense:sense ratios", col="blue")
  # label outliers in boxplot with replicate id. Left aligned.
  text(bp$group, allratios[allratios$ratio %in% bp$out,]$ratio, 
       labels=allratios[allratios$ratio %in% bp$out,]$rep, pos=4)
  
  plot <- ggplot2.stripchart(data=log10(allratios$ratio),
                             showLegend=FALSE,
                             backgroundColor="white", xtitle="All groups", ytitle="Antisense:sense ratios", 
                             mainTitle="Antisense:sense ratios",
                             addBoxplot=TRUE, boxplotFill="#E69F00")
  print(plot)
  
}

recurse_dirs <- function (start_dir, 
             fwdsuffix="fwdcounts.txt", 
             revsuffix="revcounts.txt")
{
  # get the directories
  dirlist = list.dirs(path=start_dir, recursive=TRUE, full.names=FALSE)
  
  # put them in a dataframe, calculate the lengths of each path and sort 
  # allows us to process lowest level dirs first, so that higher
  # level dirs will have access to aggregated data from lower levels
  dirs = data.frame(dirlist,stringsAsFactors = FALSE)
  dirs$length = nchar(dirs[[1]])
  dirs = dirs[with(dirs, order(-length, dirlist)), ]
  dirs = dirs[dirs$length > 0,]
  
  ratios = list()
  for (d in dirs$dirlist)
  {
    # name the group as the bottom level directory
    group = tail(strsplit(d, "/")[[1]],1)
    path = paste(start_dir, d, sep="")
    ratios[[d]] = analyse_by_group(path, group, fwdsuffix, revsuffix)
  }
  
  return(ratios)
}

#================================================================================
#'Analyse spike-in data from a single group directory
#'  
#' @param start_dir Group directory for data
#' @param group Name of group
#' @param fwdsuffix String at end of file name to indicate forward counts, default fwdcounts.txt
#' @param revsuffix String at end of file name to indicate reverse counts, default revcounts.txt
#' @param skiplines Number of lines to skip at top of each counts file, default 0
#' @param genecol Column number of gene id column, default 1
#' @param fwdcol Column number of (fwd) sense counts column, default 7
#' @param revcol Column number of (rev) antisense counts column, default 8
#' 
analyse_by_group <- function(start_dir, 
                             group, 
                             fwdsuffix="fwdcounts.txt", 
                             revsuffix="revcounts.txt", 
                             skiplines=0,
                             genecol=1,
                             lengthcol=6,
                             fwdcol=7,
                             revcol=8)
{
  # get files containing all the fwd and rev counts
  fwdnames = list.files(path=start_dir, pattern=paste("*", fwdsuffix, sep=""), full.names=TRUE)
  revnames = list.files(path=start_dir, pattern=paste("*", revsuffix, sep=""), full.names=TRUE)
  
  # return early if there are no counts files in this directory
  if (length(fwdnames) == 0)
  {
    print(paste("No counts files found in directory", group))
    return(NULL)
  }
  
  # read in all the fwd and rev counts
  fwdlabel = paste(fwdsuffix, ".*", sep="")
  prefixes = lapply(fwdnames, function(x) { gsub(fwdlabel, "", x) })
  prefixes = lapply(prefixes, function(x) { basename(x) })
  
  revfiles = lapply(revnames, read.table, "\t", header=TRUE, skip=skiplines, stringsAsFactors=FALSE)
  fwdfiles = lapply(fwdnames, read.table, "\t", header=TRUE, skip=skiplines, stringsAsFactors=FALSE)
  
  # check here that we have columns in each file matching gene/fwd/rev col - if not indicates a preprocessing problem
  rcolcheck = unlist(lapply(revfiles, function(x) { ncol(x) >= max(genecol, fwdcol, revcol)}))
  fcolcheck = unlist(lapply(fwdfiles, function(x) { ncol(x) >= max(genecol, fwdcol, revcol)}))
  if (!all(rcolcheck) | (!all(fcolcheck)))
  {
    print("Some files do not have enough columns. These files will be ignored:")
    print(fwdnames[!(fcolcheck)])
    print(revnames[!(rcolcheck)])
    
    fwdfiles = fwdfiles[fcolcheck]
    revfiles = revfiles[rcolcheck]
    
    prefixes = prefixes[fcolcheck & rcolcheck]
  }
  
  # merge sense and antisense counts
  allcounts = mapply(merge_fwd_rev, fwdfiles, revfiles, genecol, lengthcol, fwdcol, revcol, SIMPLIFY=FALSE)
  
  # set up holder for ratio results
  ratios = data.frame("rep"=as.character(prefixes), "ratio"=0, "group"=group, stringsAsFactors=FALSE)
  
  # calc ratios and plot scatters of sense vs antisense
  #par(mfrow=c(1,1))
  #ratios$ratio = mapply(calc_ratio_talk, allcounts, prefixes, SIMPLIFY=FALSE)
  par(mfrow=c(3,2))
  ratios$ratio = mapply(calc_ratio, allcounts, prefixes, SIMPLIFY=FALSE)
  
  return(ratios)
}



#================================================================================
#' Calculate antisense:sense spikein ratios for a sample,
#' printing scatterplot with fitted line of sense vs antisense counts
#'
#'@param counts The counts for the spikeins with sense and anti columns
#'@param prefix The sample id
calc_ratio <- function(counts, prefix)
{
  if (sum(counts$sense, na.rm=TRUE) == 0)
  {
    #no sense counts (!) just output a warning
    print(paste("0 sense counts in sample:",prefix))
    ratio = -1
  }
  else
  {
    #plot_spikeins(counts, NULL, prefix, 3e+06, 7e+03, aslog=FALSE)
    m <- lm(counts$anti[counts$anti!=0] ~ 0 + counts$sense[counts$anti!=0])
    ratio = (m$coefficients)
  }
  return(ratio)
}

#' #================================================================================
#' #' Calculate antisense:sense spikein ratios for a sample, 
#' #' printing scatterplot with fitted line of sense vs antisense counts
#' #'
#' #'@param counts The counts for the spikeins with sense and anti columns
#' #'@param prefix The sample id
#' calc_ratio_talk <- function(counts, prefix)
#' {
#'   if (sum(counts$sense) == 0)
#'   {
#'     #no sense counts (!) just output a warning
#'     print(paste("0 sense counts in sample:",prefix))
#'     ratio = -1
#'   }
#'   else
#'   {
#'     plot_spikeins_talk(counts, NULL, prefix, 2e+06, 2e+04)
#'     par(new=TRUE)
#'     
#'     if ((prefix=="Col-03-"))
#'     {
#'       par(new=FALSE)      
#'     }
#'     
#'     #plot_spikeins(counts, NULL, prefix, 3e+06, 7e+03, aslog=FALSE)
#'     m <- lm(counts$anti[counts$anti!=0] ~ 0 + counts$sense[counts$anti!=0])
#'     ratio = (m$coefficients)
#'   }
#'   return(ratio)
#' }
