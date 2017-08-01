##############################################################################################################
# Plot graphs of spikeins and data from Arabidopsis_RNAMeth_polyA data / 
# analyse with reference to other Arabidopsis experiments
# A version of the code in this file was used to generate plots for my Barton-Simpson group talk of 25/7/16
##############################################################################################################

source("R/analyse.R")

#' Read sense and antisense counts data
#' export
read_data <- function(sensename, antiname)
{
  sense <- read.csv(sensename, sep="\t")
  anti <- read.csv(antiname, sep="\t")
  names(sense)[c(7,8)] = c("fwd.sense", "rev.sense")
  names(anti)[c(7,8)] = c("fwd.anti", "rev.anti")
  data <- merge(sense,anti,by="Geneid")
  
  data$sense = data$fwd.sense + data$rev.sense
  data$anti = data$fwd.anti + data$rev.anti
  names(data)[c(6)] = "length"
  
  totalcounts = sum(data$sense, na.rm=TRUE) + sum(data$anti, na.rm=TRUE)
  data$tpmsense = calc_tpm(data, "sense", totalcounts)
  data$tpmanti = calc_tpm(data, "anti", totalcounts)
  
  data <- data[c("Geneid","sense", "anti","tpmsense","tpmanti","length")]
  return(data)
}

get_splices <- function(start_dir, suffix="-antisplicedoutput.csv")
{
  # get spliced counts files
  filenames = list.files(path=start_dir, pattern=paste("*", suffix, sep=""), full.names=TRUE)
  allfiles = lapply(filenames, read.csv, stringsAsFactors=FALSE)
  
  return(allfiles)
}

#================================================================================
#' Merge forward and reverse spike-in counts into one dataframe
#' 
#' @param fwdfile The file containing the forward strand counts
#' @param revfile The file containing the reverse strand counts
#' @param genecol Column number of gene id column
#' @param fwdcol Column number of (fwd) sense counts column
#' @param revcol Column number of (rev) antisense counts column
#'
merge_fwd_rev <- function(fwdfile, revfile, genecol, lengthcol, fwdcol, revcol)
{
  counts = data.frame(fwdfile[genecol], fwdfile[lengthcol], fwdfile[fwdcol], revfile[revcol])
  names(counts) = c(GENEID, LENGTH, SENSE, ANTI)
  return(counts)
}

get_spikeins <- function(start_dir, 
                         fwdsuffix="fwdcounts.txt", 
                         revsuffix="revcounts.txt", 
                         skiplines=0,
                         genecol=1,
                         fwdcol=7,
                         revcol=8,
                         lengthcol=6)
{
  # get files containing all the fwd and rev counts
  fwdnames = list.files(path=start_dir, pattern=paste("*", fwdsuffix, sep=""), full.names=TRUE)
  revnames = list.files(path=start_dir, pattern=paste("*", revsuffix, sep=""), full.names=TRUE)
  
  # return early if there are no counts files in this directory
  if (length(fwdnames) == 0)
  {
    print(paste("No counts files found in directory", start_dir))
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
  allcounts = lapply(allcounts, function(x) transform(x, ratio = anti/sense))
  return(allcounts)
}


plot_boxplots <- function(results)
{
  pdf(file="RatiosBoxplotAll.pdf",width=12,height=6)
  boxplot(as.numeric(results$ratio)~results$experiment,data=results, main="Antisense:sense ratios", 
          xlab="Experiment", ylab="Antisense:sense ratios", col="blue")
  dev.off()
  
  pdf(file="RatiosBoxplotNopolyA.pdf",width=11,height=6)
  boxplot(as.numeric(results$ratio[results$experiment!="Arabidopsis_RNAMeth_polyA"])~results$experiment[results$experiment!="Arabidopsis_RNAMeth_polyA"],data=results, main="Antisense:sense ratios", 
          xlab="Experiment", ylab="Antisense:sense ratios", col="blue")
  dev.off()
  
  pdf(file="ExoBoxplot.pdf",width=7,height=6)
  exo <- results[results$experiment=="Arabidopsis_Exosome",]
  exo$con <- c("Col","Col","Col","Col","Col","Col","Col",
               "hen","hen","hen","hen","hen","hen","hen",
               "mtr","mtr","mtr","mtr","mtr","mtr","mtr")
  boxplot(as.numeric(exo$ratio)~con,data=exo, main="Antisense:sense ratios", xlab="Condition", 
          ylab="Antisense:sense ratios (Exosome)", col="blue")
  dev.off()
  
  pdf(file="RiboBoxplot.pdf",width=7,height=6)
  ribo <- results[results$experiment=="Arabidopsis_RNAMeth_ribo",]
  ribo$con <- c("Col","Col","Col","Col","Col","Col","Col","Vir","Vir","Vir","Vir","Vir","Vir","Vir")
  boxplot(as.numeric(ribo$ratio)~con,data=ribo, main="Antisense:sense ratios", xlab="Condition", 
          ylab="Antisense:sense ratios (Ribominus)", col="blue")
  dev.off()
  
  polya <- results[results$experiment=="Arabidopsis_RNAMeth_polyA",]
  polya$con <- c("Col","Col","Col","Vir","Vir","Vir")
  boxplot(as.numeric(polya$ratio)~con,data=polya, main="Antisense:sense ratios", xlab="Condition", 
          ylab="Antisense:sense ratios (polyA)", col="blue")
}

# load splice counts data by replicate
path <- "~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/antisense-260816/"
replicates <- c("Col-01_S1", "Col-02_S1", "Col-03_S2", "Vir-1_S2", "Vir-2_S2", "Vir-3_S1")
conditions <- c("Col","Col","Col","Vir","Vir","Vir")
data <- lapply(replicates, function(r) read_data(paste(path,r,"-sensecounts.txt",sep=""), paste(path,r,"-anticounts.txt", sep="")))

# add lengths to this data too!

# Note all data files generated from featureCounts need to have the first line removed, as it messes up
# the columns. Use e.g. sed -i '' '1d' <filename> to remove first line (on a Mac, elsewhere omit the '')

# load this counts file just to get lengths of genes, and a complete list of gene ids
Col.counts.genes <- read.table("~/Documents/Arabidopsis_RNAMeth_polyA/Col_counts.txt", header=TRUE, quote="\"", stringsAsFactors = FALSE)
lengths = Col.counts.genes[c(1,6)]

# load spike-ins
sps <- get_spikeins("/Users/kmourao/Documents/ENCODE/fortalk2/polyA/")

# load splices
splices <- get_splices("~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/spliced_antisense/")

# add the lengths and full set of gene ids via outer join - this adds a load of NAs but ensures the dataframes are all the same size
splices = lapply(1:6, function(z) merge(x=splices[[z]], y=lengths, by.x="gene_id", by.y="Geneid", all.y=TRUE))
groups = c("WT1","WT2","WT3","Mutant1","Mutant2","Mutant3")

# sort
splices = lapply(1:6, function(x) splices[[x]][order(splices[[x]]$gene_id),])

# munge data into format expected by rosa TODO add functions to do this if data is different
sensesps = matrix(unlist(lapply(sps, `[`, 3)),ncol=6,byrow=FALSE)
antisps = matrix(unlist(lapply(sps, `[`, 4)),ncol=6,byrow=FALSE)
sensesplices = matrix(unlist(lapply(splices, `[`, 3)),ncol=6,byrow=FALSE)
antisplices = matrix(unlist(lapply(splices, `[`, 2)),ncol=6,byrow=FALSE)

rosa_result = rosa(data,sensesps, antisps, sensesplices, antisplices,
                  sps[[1]][1], splices[[1]][1], sps[[1]][2], splices[[1]][5], groups, 
                  resultdir="/Users/kmourao/Documents/kmourao/km-antisense/results", global=FALSE)

rosa_result = rosa(data,NULL, NULL, sensesplices, antisplices,
                   NULL, splices[[1]][1], NULL, splices[[1]][5], groups, 
                   resultdir="/Users/kmourao/Documents/kmourao/km-antisense/results", global=FALSE)

####################################################################################################
## Plot Col spikeins by replicate, aggregated from the lanes data (because that's what I've got)

plots <- function(counts, prefix)
{
  plot_spikeins_talk(counts, NULL, prefix, 2e+06, 2e+04)
  par(new=TRUE)
  
  if ((prefix=="WT3"))
  {
    par(new=FALSE)      
  }
}

pdf(file="Colspikeinsbyrep.pdf",width=7,height=6)
result <- lapply(1:3, function(x) plots(sps[[x]], groups[[x]]))
dev.off()

pdf(file="Virspikeinsbyrep.pdf",width=7,height=6)
result <- lapply(4:6, function(x) plots(sps[[x]], groups[[x]]))
dev.off()



# ####################################################################################################
# ## Plot Ribo-minus spikeins by replicate (in this case directly from the replicates data)  
# path <- "~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Arabidopsis_RNAMeth_ribo/spikeins/"
# 
# condition <- c("Col-1", "Col-2", "Col-3", "Col-4", "Col-5", "Col-6", "Col-7",
#                "Vir-1", "Vir-2", "Vir-3", "Vir-4", "Vir-5", "Vir-6", "Vir-7")
# symbol <- c(8,8,8,8,8,8,8, 21,21,21,21,21,21,21)
# colour <- c('black', 'blue', 'red', 'green', 'pink', 'yellow', 'magenta')
# isfirst <- c('TRUE', 'FALSE', 'FALSE', 'FALSE', 'FALSE', 'FALSE', 'FALSE',
#              'FALSE', 'FALSE', 'FALSE', 'FALSE', 'FALSE', 'FALSE', 'FALSE')
# details <- data.frame(condition,symbol,colour,isfirst)
# details$rep <- lapply(details$condition, function(c) load_spike_in_rep_direct(c, path))
# 
# plot_spikeins_by_rep(details, "AllRiboMinus.pdf", 'Spike-ins by replicate (Ribo-minus)',10000, 2e+06)
# 
# ####################################################################################################
# ## Analyse spike-in results across all the Arabidopsis experiments and compare to the polyA results
# 
# bpath <- "~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/"
# sp <- "/spikeins/"
# paths <- c("~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Arabidopsis_Exosome/spikeins/",
#            "~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Arabidopsis_xrn3-8/spikeins/",
#            "~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Arabidopsis_RNAMeth_ribo/spikeins/",
#            "~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Arabidopsis_RNAMeth_polyA/spikeins/")
# conditions <- list( c("Col-1", "Col-2", "Col-3", "Col-4", "Col-5", "Col-6", "Col-7",
#                       "hen-1", "hen-2", "hen-3", "hen-4", "hen-5", "hen-6", "hen-7",
#                       "mtr-1", "mtr-2", "mtr-3", "mtr-4", "mtr-5", "mtr-6", "mtr-7"),
#                     c("Col-1", "Col-2", "Col-3", "xrn-1", "xrn-2", "xrn-3"),
#                     c("Col-1", "Col-2", "Col-3", "Col-4", "Col-5", "Col-6", "Col-7",
#                       "Vir-1", "Vir-2", "Vir-3", "Vir-4", "Vir-5", "Vir-6", "Vir-7"),
#                     c("Col-01", "Col-02", "Col-03", "Vir-1", "Vir-2", "Vir-3"))
# 
# # expects counts files to be in path/condition-<fwd/rev>counts.txt
# 
# # calc experiment labels for each condition & replicate
# labels <- sapply(1:length(paths), function(p) sapply(conditions[[p]], function(c) 
#   rbind(sub(sp,"",sub(bpath,"",paths[p])))))
# 
# results <- calcratios(paths, conditions, labels)
# 
# plot_boxplots(results)
# 
# ## plot histogram of ratios for non-polyA (but include outlying polyA point to get scale right for comparisons)
# r <- as.numeric(results$ratio[results$experiment!="Arabidopsis_RNAMeth_polyA" | results$condition == "Col-03"])
# plot_histogram(r, 0.012, 20, "RatioDist.pdf")
# 
# ## plot with curve fitted to to non-polyA data
# rcurve <- as.numeric(results$ratio[results$experiment!="Arabidopsis_RNAMeth_polyA"])
# plot_histogram(r, 0.012, 20, "RatioDistWithCurve.pdf", rcurve)
# 
# ## plot all data + curve fitted to non-polyA data
# r <- as.numeric(results$ratio)
# plot_histogram(r, 0.012, 20, "RatioDistAll.pdf", rcurve)
# 
# 
# # calc probabilities of Col1 and Col3 results
# mean=mean(as.numeric(results$ratio[results$experiment!="Arabidopsis_RNAMeth_polyA"]))
# sd=sd(as.numeric(results$ratio[results$experiment!="Arabidopsis_RNAMeth_polyA"]))
# 
# probcol3 <- pnorm(0.011263, mean, sd, lower.tail=FALSE)
# probcol1 <- pnorm(0.0050597, mean, sd, lower.tail=FALSE)
# 

# Ribominus genes identification
# gtf = read.csv("/Users/kmourao/Documents/Arabidopsis_RNAMeth_polyA/genes_fixed.gtf", sep="\t", header=FALSE)
# rrnas = strsplit(as.character(gtf), " ")
# rrnas = sapply(rrnas, "[", 4)
# rrnas = strtrim(rrnas,9)
# rrnas = data.frame(rrnas, stringsAsFactors = FALSE)
# names(rrnas)[1] = c("Geneid")
# ribocounts = lapply(data, function(x) merge(x, rrnas, by="Geneid", all.y=TRUE)[c("Geneid","sense","anti")])



# Read sense and antisense spliced counts data
read_spliced_data <- function(path, repname)
{
  data <- read.csv(paste(path,repname,"-antisplicedoutput.csv",sep=""))
  names(data)[1] <- "Geneid"
  return(data)
}

# make data for ERCCdashboard
make_ERCC_dash_data <- function()
{
  path <- "~/Documents/Arabidopsis_RNAMeth_polyA/spikeins_by_lane/"
  numlanes <- 4
  condition <- c("Col-01_S1_", "Vir-1_S2_", "Vir-3_S1_", "Col-02_S1_", "Col-03_S2_", "Vir-2_S2_")
  details <- data.frame(condition)
  details$rep <- lapply(details$condition, function(c) load_spike_in_rep(c, path, numlanes))

  colnames <- c("Features", "Mix1_1", "Mix1_2", "Mix1_3", "Mix2_1", "Mix2_2", "Mix2_3")
  
  allsense <- do.call("cbind",lapply(details$rep, function(x) x[2]))
  allsense <- do.call("cbind",list(details$rep[[1]][1], allsense))
  names(allsense) <- colnames

  allanti <- do.call("cbind",lapply(details$rep, function(x) x[3]))
  allanti <- do.call("cbind",list(details$rep[[1]][1], allanti))
  names(allanti) <- colnames
  
  temp <- cbind(allsense[2:7], allanti[2:7])
  allspikeins <- sapply(unique(colnames(temp)), function(x) rowSums(temp[, colnames(temp) == x, drop = FALSE]))
  allspikeins <- do.call("cbind", list(details$rep[[1]][1], allspikeins))
  names(allspikeins) <- colnames
  
  exdat <- runDashboard(datType="count", isNorm = FALSE,
                                              exTable=allspikeins,
                                              filenameRoot="test", sample1Name="Mix1",
                                              sample2Name="Mix2", erccmix="RatioPair",
                                              erccdilution=1/100, spikeVol=8,
                                              totalRNAmass=4, choseFDR=0.1)
  
  exdat <- dynRangePlot(exdat, allPoints="TRUE", labelReps ="TRUE")
  
  return (list(allsense, allanti, allspikeins))
}

