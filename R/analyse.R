
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
#' @return the antisense:sense ratios for each group
#' 
#' export
analysenew<-function(sensecounts, anticounts, ids, groups)
{
  # validate
  paramcheck = validate_counts_data(sensecounts, anticounts, ids, groups)
  if (paramcheck$error)
  {
    stop(paramcheck$message)
  }
  
  # merge across groups
  sensegroups = sapply(split(seq_along(groups), groups), function(x) rowSums(sensecounts[, x, drop = FALSE]))
  antigroups = sapply(split(seq_along(groups), groups), function(x) rowSums(anticounts[, x, drop = FALSE]))
  
  # merge sense and antisense counts
  allcounts = lapply(1:ncol(sensegroups), function(x) makecountsdata(sensegroups[,x],anticounts[,x],ids))
  
  # calculate ratios
  # set up holder for ratio results
  ratios = data.frame("rep"=as.character(unique(groups)), "ratio"=0, stringsAsFactors=FALSE)
  
  # calc ratios of sense vs antisense
  ratios$ratio = mapply(calc_ratio, allcounts, as.list(unique(groups)), SIMPLIFY=FALSE)
  
  return(list(ratios, allcounts))
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
    counts[[seq_along(x)]], paste(x, name, sep=" "), x, xmax, ymax, aslog))
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
#' Given separate sense, anti and ids columns, make a dataframe of data with columns Geneid, sense and anti
#'
#'@param sense The counts for the sense column
#'@param anti The counts for the anti columns
#'@param ids The ids
makecountsdata <- function(sense, anti, ids)
{
  ratio = anti/sense
  rep = data.frame(ids, sense, anti, ratio)
  names(rep) <- c("Geneid", "sense", "anti", "ratio")
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
  if (sum(counts$sense) == 0)
  {
    #no sense counts (!) just output a warning
    print(paste("0 sense counts in sample:",prefix)) 
    ratio = -1
  }
  else
  {
    m <- lm(counts$anti[counts$anti!=0] ~ 0 + counts$sense[counts$anti!=0])
    ratio = (m$coefficients)[[1]]
  }
  return(ratio)
}

##############################################################################################################
#' Check that:
#' - sensecounts and anticounts have same dimensions
#' - group size = number of columns
#' - ids size = number of rows
#'
#'@param counts The counts for the spikeins with sense and anti columns
#'@param prefix The sample id
validate_counts_data<-function(sensecounts, anticounts, ids, groups)
{
  warnmsg <- list()
  
  if (nrow(sensecounts) != nrow(anticounts))
  {
    return(list("error"=TRUE, "message"="Sensecounts and anticounts have different dimensions. Their dimensions must be the same."))
  }
  else if (ncol(sensecounts) != ncol(anticounts))
  {
    return(list("error"=TRUE, "message"="Sensecounts and anticounts have different dimensions. Their dimensions must be the same."))
  }
  
  if (nrow(ids) != nrow(sensecounts))
  {
    return(list("error"=TRUE, "message"="Ids has a different number of rows than sensecounts and anticounts."))
  }
  if (length(groups) != ncol(sensecounts))
  {
    return(list("error"=TRUE, "message"="The length of the groups parameter must match the number of columns in sensecounts or anticounts."))
  }
  return(list("error"=FALSE, "message"="!"))
}