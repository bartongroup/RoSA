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
#' Analyse antisense data
#' 
#' Count and group parameters for this function can be obtained from an edgeR DGEList object d, via d$counts and d$samples$group
#'
#' export
rosa <- function(data, spikein_sense, spikein_anti, splice_sense, splice_anti, spike_ids, splice_ids, 
                 spike_lengths, splice_lengths, groups, global=FALSE)
{
  # validate
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
    paramcheck = validate_data(data, groups, splice_ids)
  }
  if (paramcheck$error)
  {
    stop(paramcheck$message)
  }
  
  # total counts
  totalcounts = colSums(rbind(colSums(spikein_sense),colSums(spikein_anti),colSums(splice_sense,na.rm=TRUE),colSums(splice_anti,na.rm=TRUE)))
  
  # calc ratios
  # sensesplices end up with NAs at e.g. position 18719 which then screws up everything
  spikeratios = calculate_ratios(spikein_sense, spikein_anti, spike_ids, spike_lengths, groups, totalcounts)
  spliceratios = calculate_ratios(splice_sense, splice_anti, splice_ids, splice_lengths, groups, totalcounts)
  
  grouplist = split(seq_along(groups), groups)
  newdata = data[order(unlist(grouplist))]
  newtotalcounts = totalcounts[order(unlist(grouplist))]
  
  # calc correction
  corrected_data = make_correction(newdata, spikeratios, spliceratios, totalcounts, global)
  
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
#' export
calculate_ratios<-function(sensecounts, anticounts, ids, lengths, groups, totalcounts)
{
  # merge across groups
  # TODO this is not working - split + seq_along is changing the group order when mutant/wt
  # but we need something like this to group replicates together when we want to
  grouplist = split(seq_along(groups), groups)
  sensegroups = sapply(grouplist, function(x) rowSums(sensecounts[, x, drop = FALSE]))
  antigroups = sapply(grouplist, function(x) rowSums(anticounts[, x, drop = FALSE]))
  
  # merge sense and antisense counts
  allcounts = lapply(1:ncol(sensegroups), function(x) makecountsdata(sensegroups[,x],antigroups[,x],ids, lengths, totalcounts[x]))
  
  # calculate ratios
  # set up holder for ratio results
  ratios = data.frame("rep"=as.character(names(grouplist)), "ratio"=0, stringsAsFactors=FALSE)
  
  # calc ratios of sense vs antisense
  ratios$ratio = mapply(calc_ratio, allcounts, as.list(names(grouplist)), SIMPLIFY=FALSE)
  
  # remove any stray NAs
  allcounts = lapply(1:length(allcounts), function(x) allcounts[[x]][complete.cases(allcounts[[x]]),])
  
  return(list(ratios, allcounts))
}

##############################################################################################################
#' Make a dataset with corrected antisense counts
#'
make_correction <- function(data, spikeratios, spliceratios, totalcounts, global=FALSE)
{
  by_rep_ratios = spikeratios[[1]]$ratio
  if (global)
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
#' @param name name of file to output to, must end in .pdf. For each plot the name will be appended to the group id for the plot.
#' @param xmax max extent of x axis
#' @param ymax max extent of y axis
#' @param aslog TRUE if plotting log/log else FALSE
#' @param show_spikeins TRUE if spikeins are also to be plotted
plot_data_and_spikein_ratios <- function(ratiodata, data, title, name, xmin = 1e-12, ymin=1e-12, xmax=0, ymax=0,
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
    paste(x, name, sep=" "), 
    xmin, 
    ymin, 
    xmax,
    ymax,
    legx, 
    legy, 
    aslog,
    show_spikeins,
    list(ratiodata[[1]]$ratio[[match(x,groups$rep)]], data[[1]]$ratio[[match(x,groups$rep)]])))
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
                                          summary_ratios=NULL)
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
             paste("Spike-in ratio: ",format(round(summary_ratios[[1]],digits=4),nsmall=4), "SD: ", 0.001),
             paste("Splices ratio: ",format(round(summary_ratios[[2]],digits=4),nsmall=4), "SD: ", 0.001)), 
             pos=2, cex=0.85)
    }
    else
    {
      legend("bottomright", bty="n", legend=parse(text=paste("R", "^2", ": ",format(round(summary(m)$adj.r.squared, digits=2),nsmall=2))), cex=0.85)
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

validate_data <- function(data, groups, ids)
{
  # data list has as many elements as there are groups
  if (length(data) != length(groups))
  {
    return(list("error"=TRUE, "message"="The number of elements in the data list must be equal to the number of groups."))
  }

  return(list("error"=FALSE, "message"="!"))
}