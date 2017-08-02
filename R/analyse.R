# This file is part of RoSA.
# 
# RoSA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# RoSA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with RoSA  If not, see <http://www.gnu.org/licenses/>.

#================================================================================
# Constants
GENEID = "Geneid"
SENSE  = "sense"
ANTI   = "anti"
LENGTH = "length"
RATIO = "ratio"
TPMSENSE = "tpmsense"
TPMANTI = "tpmanti"
TPMRATIO = "tpmratio"

##############################################################################################################
#' Analyse RNA-Seq counts data for spurious antisense
#' 
#' RoSA calculates the incidence of spurious antisense counts in RNA-Seq data, using either (or both) antisense and
#' sense counts from spike-ins and spliced reads at known splice junctions.
#' 
#' This function was written with the intention of obtaining count and group parameters from an 
#' edgeR DGEList object d, via d$counts and d$samples$group. Input data from other sources may need to be adjusted
#' to match the required format.
#'
#' @param data A list of dataframes (1 for each replicate) containing sense and antisense counts by gene. 
#' Each dataframe should have columns: 
#' \itemize{
#' \item "Geneid" (the id of each gene)
#' \item "sense" (raw sense counts for each gene)
#' \item "anti" (raw antisense counts for each gene)
#' \item "tpmsense" (normalised sense counts (TPM))
#' \item "tpmanti" (normalised antisense counts (TPM))
#' \item "length" (length of gene)
#' }
#' @param spikein_sense A matrix of integers, where each column c corresponds to a replicate, 
#' and each row r corresponds to a spike in. Each entry is the sense counts for the r-th spike-in for replicate c.
#' The c-th column corresponds to the c-th entry in data. If there is no spike-in data, set to NULL.
#' @param spikein_anti A matrix of integers, where each column c corresponds to a replicate, 
#' and each row r corresponds to a spike in. Each entry is the antisense counts for the r-th spike-in for replicate c.
#' The c-th column corresponds to the c-th entry in data. If there is no spike-in data, set to NULL.
#' @param splice_sense A matrix of integers, where each column c corresponds to a replicate, 
#' and each row r corresponds to a gene. Each entry is the spliced sense counts (may be 0 or NA) for the r-th 
#' gene for replicate c. The c-th column corresponds to the c-th entry in data. If there is no spliced reads data, set to NULL.
#' @param splice_anti A matrix of integers, where each column c corresponds to a replicate, 
#' and each row r corresponds to a gene. Each entry is the spliced sense counts (may be 0 or NA) for the r-th 
#' gene for replicate c. The c-th column corresponds to the c-th entry in data. If there is no spliced reads data, set to NULL.
#' @param spike_ids A dataframe containing spike-in ids in column 'Geneid'
#' @param splice_ids A dataframe containing gene ids in column 'Geneid'
#' @param spike_lengths A dataframe containing spike-in lengths in column 'length'
#' @param splice_lengths A dataframe containing gene lengths in column 'length'
#' @param groups A character array listing each replicate by its position in data, and the spikein and 
#' spliced reads matrices. E.g. if replicates in data are in the order WT1,WT2,WT3,Mutant1,Mutant2,Mutant3 then 
#' groups should be set to c("WT1","WT2","WT3","Mutant1","Mutant2","Mutant3"). In future this parameter
#' will allow replicates to be combined in the analysis.
#' @param resultdir Full path to directory where plots should be output
#' @param global Use only spike-in ratios to calculate correction (Default \code{FALSE})
#' @param xmin Minimum x-value for plots of antisense vssense counts (Default \code{1e-12})
#' @param xmax Maximum x-value for plots of antisense vs sense counts (Default \code{1e-2})
#' @param ymin Minimum y-value for plots of antisense vs sense counts (Default \code{1e-12})
#' @param ymax Maximum y-value for plots of antisense vs sense counts (Default \code{1e-4})
#' @param legendx Location of left side of legend (Default \code{-12})
#' @param legendy Location of top of legend (Default \code{0})
#' @return Returns a list with entries:
#' \itemize{
#' \item corrected_data (original data with corrected counts values),
#' \item spikeratios (ratios for the spike-in data), and 
#' \item spliceratios (ratios for the spliced reads data).
#' }
#' Each of spikeratios and spliceratios consists of a list with two entries: the first entry is a dataframe 
#' containing the global spike-in/spliced ratio of antisense:sense, and the standard deviation of the residuals 
#' for the associated linear model, for each condition listed in groups; the second entry is a list of dataframes,
#' one for each condition listed in groups, containing counts and ratios for each spike-in/gene.
#' @export
rosa <- function(data, spikein_sense, spikein_anti, splice_sense, splice_anti, spike_ids, splice_ids, 
                 spike_lengths, splice_lengths, groups, resultdir, global=FALSE, xmin=1e-12, xmax=1e-2, ymin=1e-12, ymax=1e-4,
                 legendx=-12, legendy=0)
{
  # validate
  message("Validating parameters")
  paramcheck = validate_ratio_data(spikein_sense, spikein_anti, spike_ids, spike_lengths, groups, "Spike-ins")
  if (paramcheck$error)
  {
    stop(paramcheck$message)
  }
  else
  {
    paramcheck = validate_ratio_data(splice_sense, splice_anti, splice_ids, splice_lengths, groups, "Spliced reads")
  }
  if (paramcheck$error)
  {
    stop(paramcheck$message)
  }
  else
  {
    paramcheck = validate_data(data, groups)
  }
  
  if (paramcheck$error)
  {
    stop(paramcheck$message)
  }
  else
  {
    if (!dir.exists(resultdir))
    {
      paramcheck = list("error"=TRUE, "message"=paste(resultdir, 
                          " does not exist on the filesystem. Please supply *the full path* of an existing directory for results."))
    }
    if (is.null(splice_sense) & is.null(spikein_sense))
    {
      paramcheck = list("error"=TRUE, "message"="Both splice counts and spike-in counts are NULL. At least one must be non-NULL.")
    }
  }
  
  if (paramcheck$error)
  {
    stop(paramcheck$message)
  }
  
  message("Calculating ratios")
  
  
  if (is.null(spikein_sense))
  {
    # total counts
    totalcounts = colSums(rbind(colSums(splice_sense,na.rm=TRUE),colSums(splice_anti,na.rm=TRUE)))
    
    # calc ratios
    spikeratios = NULL
    spliceratios = calculate_ratios(splice_sense, splice_anti, splice_ids, splice_lengths, groups, totalcounts)
  }
  else if (is.null(splice_sense))
  {
    # total counts
    totalcounts = colSums(rbind(colSums(spikein_sense),colSums(spikein_anti)))
    
    # calc ratios
    spliceratios = NULL
    spikeratios = calculate_ratios(spikein_sense, spikein_anti, spike_ids, spike_lengths, groups, totalcounts)
  }
  else
  {
    # total counts
    totalcounts = colSums(rbind(colSums(spikein_sense),colSums(spikein_anti),
                                colSums(splice_sense,na.rm=TRUE),colSums(splice_anti,na.rm=TRUE)))
    # calc ratios
    spikeratios = calculate_ratios(spikein_sense, spikein_anti, spike_ids, spike_lengths, groups, totalcounts)
    spliceratios = calculate_ratios(splice_sense, splice_anti, splice_ids, splice_lengths, groups, totalcounts)
  }

  grouplist = split(seq_along(groups), groups)
  newdata = data[order(unlist(grouplist))]
  newtotalcounts = totalcounts[order(unlist(grouplist))]
  
  # calc correction
  message("Calculating antisense correction")
  corrected_data = make_correction(newdata, spikeratios, spliceratios, totalcounts, global)
  
  # make plots
  message("Making plots...")
  grouplist = split(seq_along(groups), groups)
  newdata = data[order(unlist(grouplist))]
  
  if (is.null(spikein_sense))
  {
    useratios = spliceratios 
  }
  else
  {
    useratios = spikeratios
  }
  plot_data_and_spikein_ratios(useratios, list(1,newdata), title="Original", "original data.pdf", resultdir,
                               xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, legx=legendx, legy=legendy, aslog=TRUE, show_spikeins=FALSE)
  
  plot_data_and_spikein_ratios(useratios, list(1,corrected_data), title="Corrected", "corrected data.pdf", resultdir,
                               xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, legx=legendx, legy=legendy, aslog=TRUE, show_spikeins=FALSE)
  
  # plot spikeins sense v antisense and draw best fit line
  # not sure how useful this is
  # plot_spikein_ratios(spikeratios, "Spike-in.pdf", xmax=2e6, ymax=1e4, aslog=TRUE)
  
  # Plot spliced sense v antisense data + spike-ins
  if (!is.null(spikein_sense))
  {
    plot_data_and_spikein_ratios(spikeratios, spliceratios, title="Spike-ins overlaid on spliced\n", "data and spikeins.pdf", resultdir,
                               xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, legx=legendx, legy=legendy, aslog=TRUE)
  }
  
  message("Finished")
  return(list(corrected_data,spikeratios,spliceratios))
}


##############################################################################################################
#' Calculate antisense:sense ratios
#' 
#' Count and group parameters for this function can be obtained from an edgeR DGEList object d, via d$counts and d$samples$group
#' 
#' @param sensecounts Sense strand counts as a matrix where each row corresponds to a location - a gene or spike-in - and 
#' each column corresponds to a source - lane, replicate, condition or experiment.
#' @param anticounts Antisense strand counts as a matrix. Locations and sources must match sensecounts.  
#' @param ids Id of each row in sensecounts/anticounts e.g. gene ids or spike-in ids
#' @param groups An array indicating the group of each column in sensecounts and anticounts.
#' @param totalcounts The total read counts for each source
#' @return the antisense:sense ratios for each group
#' 
calculate_ratios<-function(sensecounts, anticounts, ids, lengths, groups, totalcounts)
{
  # merge across groups - note not currently working as intended, ordering is not 
  # being changed to correspond to grouping
  grouplist = split(seq_along(groups), groups)
  sensegroups = sapply(grouplist, function(x) rowSums(sensecounts[, x, drop = FALSE]))
  antigroups = sapply(grouplist, function(x) rowSums(anticounts[, x, drop = FALSE]))
  
  # merge sense and antisense counts
  allcounts = lapply(1:ncol(sensegroups), function(x) 
                makecountsdata(sensegroups[,x],antigroups[,x],ids, lengths, totalcounts[x]))
  
  # calculate ratios
  # set up holder for ratio results
  ratios = data.frame("rep"=as.character(names(grouplist)), ratio="0", rmse="0", stringsAsFactors=FALSE)
  
  # calc ratios of sense vs antisense
  ratios$ratio = mapply(calc_ratio, allcounts, as.list(names(grouplist)), SIMPLIFY=FALSE)
  ratios$rmse = mapply(calc_rmse, allcounts, as.list(names(grouplist)), SIMPLIFY=FALSE)
  
  # remove any stray NAs
  allcounts = lapply(1:length(allcounts), function(x) allcounts[[x]][complete.cases(allcounts[[x]]),])
  
  return(list(ratios, allcounts))
}

##############################################################################################################
#' Make a dataset with corrected antisense counts
#'
make_correction <- function(data, spikeratios, spliceratios, totalcounts, global=FALSE)
{
  # if we don't have spike-in data, use the splices data
  if (is.null(spikeratios))
  {
    by_rep_ratios = spliceratios[[1]]$ratio
  }
  else
  {
    by_rep_ratios = spikeratios[[1]]$ratio
  }
  
  # if we don't have splices data, just apply the spike-in ratios globally
  if (global | is.null(spliceratios))
  {
    result <- lapply(1:length(data), function(i) insert_global_ratio(data[[i]], by_rep_ratios[[i]]))
  }
  else
  {
    by_gene_ratios = lapply(spliceratios[[2]], function(x) data.frame(x$Geneid,x$tpmratio))
    by_gene_ratios = lapply(by_gene_ratios, setNames, c(GENEID,RATIO))
    
    # do correction with individual per replicate ratios
    result <- lapply(1:length(data), function(i) insert_indiv_ratio(data[[i]], by_gene_ratios[[i]], by_rep_ratios[[i]]))
  }
  
  # sort out column names
  result <- lapply(result, setNames, c(GENEID, LENGTH, TPMSENSE, TPMANTI, SENSE, ANTI))
  
  return(result)
}

##############################################################################################################
#' Plot spike-in ratios for each of a set of groups
#'
#' @param ratiodata a list (such as that returned by analyse) containing the elements:
#' ratios: a data frame with columns rep and ratio. Rep = group id, ratio = antisense:sense ratio for that replicate
#' counts: list of data frames, 1 for each group in ratios$rep where each data frame has columns Geneid, sense, anti, ratio
#' Geneid = id of gene/spikein, sense = sense counts, anti = antisense counts, ratio = ratio of antisense to sense counts
#' @param name name of file to output to, must end in .pdf. For each plot the name will be appended to the group id for the plot.
#' @param xmax max extent of x axis
#' @param ymax max extent of y axis
#' @param aslog TRUE if plotting log/log else FALSE
plot_spikein_ratios <- function(ratiodata, name, xmax, ymax, aslog=TRUE)
{
  groups = ratiodata[[1]]
  counts = ratiodata[[2]]
  
  # for each group, plot the spike-in ratios to the pdf file <rep><space><name>
  ignore = lapply(groups$rep, function(x) plot_single_spikein_group(
    counts[[match(x,groups$rep)]], paste(x, name, sep=" "), x, xmax, ymax, aslog))
}

plot_single_spikein_group<- function(sp, name, label, xmax, ymax, aslog=TRUE)
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
    graphics.off()
  }
}

##############################################################################################################
#' Plot data and spike-in ratios for each of a set of groups
#'
#' @param ratiodata a list (such as that returned by analyse) containing the elements:
#' ratios: a data frame with columns rep and ratio. Rep = group id, ratio = antisense:sense ratio for that replicate
#' counts: list of data frames, 1 for each group in ratios$rep where each data frame has columns Geneid, sense, anti, ratio
#' Geneid = id of gene/spikein, sense = sense counts, anti = antisense counts, ratio = ratio of antisense to sense counts
#' @param data a list of dataframes where each dataframe contains sense and antisense counts for one group, in columns
#' Geneid, sense, anti, ratio. Entry x in the list must correspond to entry x in ratiodata$counts
#' @param name name of file to output to, must end in .pdf. For each plot the name will be appended 
#' to the group id for the plot.
#' @param resultdir Full path to an existing directory where plots will be output to
#' @param xmax max extent of x axis
#' @param ymax max extent of y axis
#' @param aslog TRUE if plotting log/log else FALSE
#' @param show_spikeins TRUE if spikeins are also to be plotted
plot_data_and_spikein_ratios <- function(ratiodata, data, title, name, resultdir, xmin = 1e-12, ymin=1e-12, xmax=0, ymax=0,
                                         legx=0, legy=0, aslog=TRUE, show_spikeins=TRUE)
{
  groups = ratiodata[[1]]
  spikeins = ratiodata[[2]]
  splices = data[[2]]
  
  # for each group, plot the spike-in ratios to the pdf file <rep><space><name>
  ignore = lapply(groups$rep, function(x) plot_single_data_and_spikeins(
    spikeins[[match(x,groups$rep)]], 
    splices[[match(x,groups$rep)]], 
    x, 
    title,
    file.path(resultdir,paste(x, name, sep=" ")),
    xmin, 
    ymin, 
    xmax,
    ymax,
    legx, 
    legy, 
    aslog,
    show_spikeins,
    list(ratiodata[[1]]$ratio[[match(x,groups$rep)]], data[[1]]$ratio[[match(x,groups$rep)]]),
    list(ratiodata[[1]]$rmse[[match(x,groups$rep)]], data[[1]]$rmse[[match(x,groups$rep)]])))
}

##############################################################################################################
#' plot_antisensedata_and_spikeins: plot data and spikeins on same graph, sense vs antisense
#' @param sp spike in antisense and sense counts
#' @param data spliced antisense and sense counts
#' @param label condition or replicate name for graph labels
#' @param name name of pdf file to output to, must end in .pdf
#' @param aslog plot as log/log if TRUE, else linear
#' @param show_spikeins TRUE if spikeins are also to be plotted
#' @param summary_ratios list containing spikeratio (global spike-in ratio) and spliceratio (global splice ratio)
plot_single_data_and_spikeins <- function(sp, 
                                          data, 
                                          label, 
                                          title = "",
                                          name, 
                                          xmin = 1e-12, 
                                          ymin = 1e-12, 
                                          xmax = 0,
                                          ymax = 0,
                                          legx=0, 
                                          legy=0, 
                                          aslog=TRUE,
                                          show_spikeins=TRUE,
                                          summary_ratios=NULL,
                                          summary_rmse=NULL)
{
  # deal with log of zero values
  if (aslog & (xmax==0))
  {
    xmax = 1
  }
  if (aslog & (ymax==0))
  {
    ymax = 1
  }
  
  plot_title <- paste(label, ": ", title, " ", "normalised antisense vs sense counts", sep="")
  
  # symbols
  data_pch <- 21  # symbol for data points - circle
  sp_pch <- 21    # symbol for spike-in points
  
  legend_labels <- c(paste(label, " spike-ins", sep=""), paste(label, " spike-in ratio", sep=""))
  legend_colours <- c("black", "black")
  legend_bg <- c("black", NA)
  legend_lty <- c(NA, 1)
  legend_pch <- c(sp_pch, NA)
  
  if (!aslog & legx==0 & legy==0)
  {
    legx = 0.6 * xmin
    legy = 0.4 * ymin
  }
  else if (legx==0 & legy==0)
  {
    legx = 0.6 * log10(xmin)
    legy = 0.4 * log10(ymin)
  }
  
  pdf(file=name,width=7,height=6)
  if (!aslog)
  {
    plot(data$sense, data$anti, pch=data_pch, cex=0.5, col="black", 
         xlab='Sense counts', ylab='Antisense counts', xlim = c(0,xmax), ylim = c(0,ymax))
    
    if (show_spikeins)
    {
      par(new=TRUE)
      plot(sp$sense, sp$anti, pch=sp_pch, cex=0.8, col='blue', bg='cyan', axes = FALSE, xlab = '', ylab = '', 
           xlim = c(xmin,xmax), ylim = c(ymin,ymax))
      abline(lm(sp$anti ~ 0 + sp$sense), col='blue')
      legend(legx, legy, legend_labels, col = legend_colours, pt.bg = legend_bg, lty = legend_lty, pch = legend_pch)
    }
  }
  else
  {
    heatscatter(log10(data$tpmsense),log10(data$tpmanti), #colpal="spectral",
                cor=FALSE,cex.main=1, xlim=c(log10(xmin),log10(xmax)), ylim=c(log10(ymin),log10(ymax)),
                xlab=expression('log'[10]*'(Sense TPM)'), ylab=expression('log'[10]*'(Antisense TPM)'), main='')
    
    use <- is.finite(log10(data$tpmanti)) & is.finite(log10(data$tpmsense))
    m = lm(log10(data$tpmanti[use])~log10(data$tpmsense[use]))
    
    if (show_spikeins)
    {
      legwidth = strwidth("Spike-in ratio: 0.0000")
      temp <-legend("bottomright", bty="n", legend = c(" "," ", " "), text.width=legwidth)
      text(temp$rect$left + temp$rect$w, temp$text$y, 
           c(parse(text=paste("R", "^2", ": ",format(round(summary(m)$adj.r.squared,digits=2),nsmall=2))),
             paste("Spike-in ratio: ",format(round(summary_ratios[[1]],digits=4),nsmall=4), "SD: ", format(round(summary_rmse[[1]],digits=8),nsmall=8)),
             paste("Splices ratio: ",format(round(summary_ratios[[2]],digits=4),nsmall=4), "SD: ", format(round(summary_rmse[[2]],digits=8),nsmall=8))), 
             pos=2, cex=1)
    }
    else
    {
      legend("bottomright", bty="n", legend=parse(text=paste("R", "^2", ": ",format(round(summary(m)$adj.r.squared, digits=2),nsmall=2))), cex=1)
    }

    if (show_spikeins)
    {
      par(new=TRUE)
  
      plot(log10(sp$tpmsense), log10(sp$tpmanti), pch=sp_pch, cex=1, col='white', bg='black', axes = FALSE, xlab = '', ylab = '',
           xlim=c(log10(xmin),log10(xmax)), ylim=c(log10(ymin),log10(ymax)))
      use <- is.finite(log10(sp$tpmanti)) & is.finite(log10(sp$tpmsense))
      abline(lm(log10(sp$tpmanti[use]) ~ log10(sp$tpmsense[use])), col='black')
      legend(legx, legy, legend_labels, col = legend_colours, pt.bg = legend_bg, lty = legend_lty, pch = legend_pch)
    }
    
  }
  if (show_spikeins)
  {
    title(main=plot_title, line=-2.2)
  }
  else
  {
    title(main=plot_title, line=-1.2)
  }
  dev.off()
}

#' Plot sense/antisense expression boxplots of a set of reads
#' @param data list of dataframes which must each have a column called colname
#' @param colname name of column to plot
boxplot_counts <- function(data, groups, sensename=SENSE, antiname=ANTI, datasetname)
{

  plotdata = sapply(data, function(x) x[[sensename]][x[[sensename]] != 0])
  title = paste("Sense expression levels (", datasetname, ")")
  boxplot(plotdata, log="y",names=groups, main=title, 
          xlab="Replicate", ylab="Number of reads", col="blue")
  plotdata = sapply(data, function(x) x[[antiname]][x[[sensename]] != 0]+1)
  title = paste("Antisense expression levels (", datasetname, ")")
  boxplot(plotdata, log="y",names=groups, main=title, 
          xlab="Replicate", ylab="Number of reads", col="blue")

}

##############################################################################################################
#' Given separate sense, anti and ids columns, make a dataframe of data with columns Geneid, sense and anti
#'
#'@param sense The counts for the sense column
#'@param anti The counts for the anti columns
#'@param ids The ids
makecountsdata <- function(sense, anti, ids, lengths, totalcounts)
{
  ratio = anti/sense
  ratio[sense == 0] = 1 # keep sense and antisense counts the same when sense counts = 0
  rep = data.frame(ids, sense, anti, lengths, ratio)
  names(rep) <- c(GENEID, SENSE, ANTI, LENGTH, RATIO)
  tpmsense = calc_tpm(rep, SENSE, totalcounts)
  tpmanti = calc_tpm(rep, ANTI, totalcounts)
  tpmratio = tpmanti/tpmsense
  tpmratio[sense == 0] = 1
  rep = data.frame(ids, sense, anti, lengths, ratio, tpmsense, tpmanti, tpmratio)
  names(rep) <- c(GENEID, SENSE, ANTI, LENGTH, RATIO, TPMSENSE, TPMANTI, TPMRATIO)
  return(rep)
}

##############################################################################################################
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
    use <- is.finite(log10(counts$tpmanti)) & is.finite(log10(counts$tpmsense))
    m <- lm(counts$tpmanti[use] ~ 0 + counts$tpmsense[use])
    ratio = (m$coefficients)[[1]]
  }
  return(ratio)
}

#' Calculate the standard deviation of the residuals
#' @param counts dataframe containing TPM sense and antisense counts
#' @param prefix sample name
#' @return RMSE for the linear regression model calculated for antisense vs sense counts
calc_rmse <- function(counts, prefix)
{
  if (sum(counts$sense, na.rm=TRUE) == 0)
  {
    #no sense counts (!) just output a warning
    print(paste("0 sense counts in sample:",prefix)) 
    ratio = -1
  }
  else
  {
    use <- is.finite(log10(counts$tpmanti)) & is.finite(log10(counts$tpmsense))
    m <- lm(counts$tpmanti[use] ~ 0 + counts$tpmsense[use])
    stderr = coef(summary(m))[,2]
    RSS <- c(crossprod(m$residuals))
    MSE <- RSS / length(m$residuals)
    RMSE <- sqrt(MSE)
  }
  return(RMSE)
}

calc_tpm <- function(counts, name, totalcounts)
{
  # divide by a million
  scalefactor = totalcounts / 10e6
  
  # calc reads per kilobase / scalefactor = tpm
  tpm = ((counts[[name]] / counts[[LENGTH]]) / 10e4) / scalefactor
  
  return(tpm)
}

invert_tpm <- function(newcounts, totalcounts)
{
  scalefactor = totalcounts / 10e6
  
  newcounts$sense = round(newcounts$tpmsense * scalefactor * 10e4 * newcounts[[LENGTH]])
  newcounts$anti = round(newcounts$tpmanti * scalefactor * 10e4 * newcounts[[LENGTH]])
  return(newcounts)
}

############################################################################################################
# Apply antisense correction ratio on a per gene basis
insert_indiv_ratio <- function(d, ratio, globalratio)
{
  d$ratio <- rep(globalratio, length(d$Geneid))
  d$ratio[na.omit(match(ratio$Geneid,d$Geneid))] <- ratio$ratio[which(ratio$Geneid %in% d$Geneid)]
  d <- apply_correction(d, "ratio")
  
  d = d[c(GENEID,LENGTH,"adjtpmsense","adjtpmanti","adjsense","adjanti")]
  return(d)
}

# Apply antisense correction ratio globally across all genes
insert_global_ratio <- function(d, ratio)
{
  d$correction <- rep(ratio, length(d$Geneid))
  d <- apply_correction(d, "correction")
  d = d[c(GENEID,LENGTH,"adjtpmsense","adjtpmanti","adjsense","adjanti")]
  return(d)
}

# Apply correction to both strands
apply_correction <- function(d, correction)
{
  # calculate antisense as a proportion of sense, proportional to correction ratio
  d$adjtpmanti[!is.na(d[correction])] <- pmax(0, d$tpmanti[!is.na(d[correction])] - 
                                                d$tpmsense[!is.na(d[correction])] * d[correction][!is.na(d[correction])])
  
  # take the difference between the new and old antisense and add it to the old sense to get an adjusted sense value
  # these are reads which were assigned to antisense but should have been sense reads
  d$adjtpmsense[!is.na(d[correction])] <- d$tpmsense[!is.na(d[correction])] + 
                                          (d$tpmanti[!is.na(d[correction])] - d$adjtpmanti[!is.na(d[correction])])
  
  # calculate antisense as a proportion of sense, proportional to correction ratio
  d$adjanti[!is.na(d[correction])] <- pmax(0, round(d$anti[!is.na(d[correction])] - 
                                                d$sense[!is.na(d[correction])] * d[correction][!is.na(d[correction])]))
  
  # take the difference between the new and old antisense and add it to the old sense to get an adjusted sense value
  # these are reads which were assigned to antisense but should have been sense reads
  d$adjsense[!is.na(d[correction])] <- d$sense[!is.na(d[correction])] + 
                                        (d$anti[!is.na(d[correction])] - d$adjanti[!is.na(d[correction])])
  
  return(d)
}

##############################################################################################################
#' Check that:
#' - sensecounts and anticounts have same dimensions
#' - group size = number of columns
#' - ids size = number of rows
#' - lengths size = number of rows; lengths are all non-zero
#' - totalcounts columns = number of groups
#'
#'@param counts The counts for the spikeins with sense and anti columns
#'@param prefix The sample id
validate_ratio_data<-function(sensecounts, anticounts, ids, lengths, groups, dataset)
{
  warnmsg <- list()
  
  if (is.null(sensecounts) | (is.null(anticounts)))
  {
    return(list("error"=FALSE,"message"="!"))
  }

  if (nrow(sensecounts) != nrow(anticounts))
  {
    return(list("error"=TRUE, "message"=paste(dataset, ": Sensecounts and anticounts have different dimensions. Their dimensions must be the same.")))
  }
  else if (ncol(sensecounts) != ncol(anticounts))
  {
    return(list("error"=TRUE, "message"=paste(dataset, ": Sensecounts and anticounts have different dimensions. Their dimensions must be the same.")))
  }

  if (nrow(lengths) != nrow(sensecounts))
  {
    return(list("error"=TRUE, "message"=paste(dataset, ": Length vector has a different number of rows to sensecounts and anticounts.")))
  }

  if (any(lengths==0))
  {
    return(list("error"=TRUE, "message"=paste(dataset, ": Length vector has at least one zero value. Lengths must be non-zero.")))
  }

  if (nrow(ids) != nrow(sensecounts))
  {
    return(list("error"=TRUE, "message"=paste(dataset, ": Ids has a different number of rows than sensecounts and anticounts.")))
  }
  if (length(groups) != ncol(sensecounts))
  {
    return(list("error"=TRUE, "message"=paste(dataset, ": The length of the groups parameter must match the number of columns in sensecounts or anticounts.")))
  }

  return(list("error"=FALSE, "message"="!"))
}

validate_data <- function(data, groups)
{
  # data list has as many elements as there are groups
  if (length(data) != length(groups))
  {
    return(list("error"=TRUE, "message"="The number of elements in the data list must be equal to the number of groups."))
  }

  return(list("error"=FALSE, "message"="!"))
}