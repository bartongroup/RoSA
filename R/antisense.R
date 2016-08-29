## Functions to support analysis, modelling and removal of antisense noise in RNA-Seq data

# data: list of dataframes with sense/antisense counts by fwd/rev strand
# names: name of each data item in data
# groups: group (replicate for lanes; condition for replicates etc) for each data item in data
# adjustment: correction factor for antisense counts, either a list of ratios, one for each data item;
# or a list of vectors of ratios, one vector for each data item, with one entry per gene in data
# return data with antisense counts corrected as per adjust values
analyse_antisense <- function(data, names, groups, adjustment)
{
  
}




## Expects data in the form of a dataframe with columns containing:
# Geneid, fwd sense counts, fwd antisense counts, rev sense counts, rev antisense counts
# e.g. create by:
# vir <- read.csv("~/Documents/Arabidopsis_RNAMeth_polyA/Vir.counts.txt", sep="\t")
# viranti <- read.csv("~/Documents/Arabidopsis_RNAMeth_polyA/Vir.anticounts.txt", sep="\t")
# names(vir)[c(7,8)] = c("fwd.sense", "rev.sense")
# names(viranti)[c(7,8)] = c("fwd.anti", "rev.anti")
# vir <- merge(vir,viranti,by="Geneid")
# data <- vir[c("Geneid","fwd.sense","rev.sense","fwd.anti","rev.anti","Length.x")]

# data: dataframe
# id_col: gene id column index
# fwdsense_col: fwd sense counts column index
# fwdanti_col: fwd antisense counts column index
# revsense_col: rev sense counts column index
# revanti_col: rev antisense counts column index
analyse_antisense <- function(data, id_col, fwdsense_col, fwdanti_col, revsense_col, revanti_col, length)
{
  names(data)[c(id_col, fwdsense_col, fwdanti_col, revsense_col, revanti_col, length)] <-
    c("Geneid", "fwd.sense", "fwd.anti", "rev.sense", "rev.anti", "length")
  
  # calc forward and reverse antisense ratios
  data$revratio <- data$fwd.anti/data$rev.sense
  data$fwdratio <- data$rev.anti/data$fwd.sense
  
  hist(log10(data$revratio[is.finite(data$revratio)]), breaks=100, col='cornflowerblue', 
       xlab='log(antisense:sense ratio)', ylab='Number of genes', main='Distribution of antisense:sense ratios by gene - reverse strand')
  hist(log10(data$fwdratio[is.finite(data$fwdratio)]), breaks=100, col='cornflowerblue', 
       xlab='log(antisense:sense ratio)', ylab='Number of genes', main='Distribution of antisense:sense ratios by gene - forward strand')

  plot(log10(vir$sensecounts), log10(vir$anticounts),  col='red', xlab='log(Sense counts)', ylab='log(Antisense counts)', xlim =c(0,7) , ylim = c(0,6))
  par(new=TRUE)
  plot(log10(data$fwd.sense), log10(data$rev.anti), xlab='Log(Sense counts)', ylab='Log(Antisense counts)', 
       main='Antisense:sense - fwd strand', xlim =c(0,7) , ylim = c(0,6))
  abline(log10(5),0,col="red")
  
  anti <- stack(data, select = c("rev.anti", "fwd.anti"))
  anti$len <- rep(data$length,2)
  par(mfrow = c(2, 2))
  
  hist(log10(anti$values), breaks=150, col='cornflowerblue')
  hist(log10(anti$values/anti$len), breaks=150, col='cornflowerblue')
  
  sense <- stack(data, select = c("rev.sense", "fwd.sense"))
  sense$len <- rep(data$length,2)
  
  hist(log10(sense$values), breaks=150, col='cornflowerblue')
  hist(log10(sense$values/sense$len), breaks=150, col='cornflowerblue')
  
  scalefactor <- 0.0011
  antisense_thresh <- 5
  scale_antisense(scalefactor, antisense_thresh, data)
  
  adjanti <- stack(data, select=c("adjfwd","adjrev"))
  adjanti$len <- rep(data$length,2)
  hist(log10(adjanti$values), breaks=150, col='cornflowerblue')
  hist(log10(adjanti$values/adjanti$len), breaks=150, col='cornflowerblue')
  
  par(mfrow=c(1,1))
  hist(log10(anti$values), breaks=150, col=rgb(1,1,0,0.7), xlim=c(0,5), ylim=c(0,2500))
  par(new=TRUE)
  hist(log10(adjanti$values), breaks=150, col=rgb(0,1,1,0.4), xlim=c(0,5), ylim=c(0,2500), xlab="", ylab="", main="")
  

}

scale_antisense <- function(scalefactor, antisense_thresh, data)
{
  data$scaledrevanti <- pmax(0, round(data$rev.anti - scalefactor * data$fwd.sense))
  data$scaledfwdanti <- pmax(0, round(data$fwd.anti - scalefactor * data$rev.sense))
  plot(log10(data$fwd.sense), log10(data$scaledrevanti), xlab='Log(Sense counts)', ylab='Log(Scaled antisense counts)', 
       main='Antisense:sense - fwd strand', xlim =c(0,7) , ylim = c(0,6))
  abline(log10(antisense_thresh),0,col="red")
  
  plot(log10(data$rev.sense), log10(data$scaledfwdanti), xlab='Log(Sense counts)', ylab='Log(Scaled antisense counts)', 
       main='Antisense:sense - rev strand', xlim =c(0,7) , ylim = c(0,6))
  abline(log10(antisense_thresh),0,col="red")
}
