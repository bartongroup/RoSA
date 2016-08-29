##############################################################################################################
# Plot graphs of spikeins and data from Arabidopsis_RNAMeth_polyA data / 
# analyse with reference to other Arabidopsis experiments
# A version of the code in this file was used to generate plots for my Barton-Simpson group talk of 25/7/16
##############################################################################################################

source("analyse-spikeins.R") # main plotting code


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




# Note all data files generated from featureCounts need to have the first line removed, as it messes up
# the columns. Use e.g. sed -i '' '1d' <filename> to remove first line (on a Mac, elsewhere omit the '')
# load spike-ins
sp <- read.csv("~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Col_spikein_counts.csv")
vsp <- read.csv("~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Vir_spikein_counts.csv")

# load data in form of spliced reads which match intron-exon structure on sense strand
col <- read.csv("~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Col-countsofsplicedantisense.csv")
vir <- read.csv("~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Vir-countsofsplicedantisense.csv")

# rename columns appropriately
names(vsp)[7] <- "sense"
names(vsp)[8] <- "anti"
names(sp)[7] <- "sense"
names(sp)[8] <- "anti"

# plot Col and Vir spikeins sense v antisense and draw best fit line
plot_spikeins(sp, "Colspikeins.pdf", "Col", 10000, 2e+06)
plot_spikeins(vsp, "Virspikeins.pdf", "Vir", 2000, 2e+06)

# Plot Col and Vir sense v antisense data + spike-ins
plot_antisensedata_and_spikeins(col, sp, "Col", "Col.pdf", aslog=FALSE)
plot_antisensedata_and_spikeins(vir, vsp, "Vir", "Vir.pdf", aslog=FALSE)

# Plot all on one plot
plot_all(col,vir,sp,vsp,"Col","Vir","ColVir.pdf",60000,5e+06)

# Plot sense and antisense reads by spikein ordered by sense reads
plot_spikeins_by_id(sp,vsp,"Col","Vir")

####################################################################################################
## Plot Col spikeins by replicate, aggregated from the lanes data (because that's what I've got)
path <- "~/Documents/Arabidopsis_RNAMeth_polyA/spikeins_by_lane/"

numlanes <- 4
condition <- c("Col-01_S1_", "Col-02_S1_", "Col-03_S2_")
symbol <- c(8,8,8)
colour <- c('black', 'blue', 'red')
isfirst <- c('TRUE', 'FALSE', 'FALSE')

details <- data.frame(condition,symbol,colour,isfirst)
details$rep <- lapply(details$condition, function(c) load_spike_in_rep(c, path, numlanes))
plot_spikeins_by_rep(details, "Colspikeinsbyrep.pdf", 'Spike-ins by replicate (Col)',10000, 2e+06)

## Plot Vir spikeins by replicate, aggregated from the lanes data
condition <- c("Vir-1_S2_", "Vir-2_S2_", "Vir-3_S1_")
details <- data.frame(condition,symbol,colour,isfirst)
details$rep <- lapply(details$condition, function(c) load_spike_in_rep(c, path, numlanes))

plot_spikeins_by_rep(details, "Virspikeinsbyrep.pdf", 'Spike-ins by replicate (Vir)',10000, 2e+06)

####################################################################################################
## Plot Col with Vir overlaid in grey
path <- "~/Documents/Arabidopsis_RNAMeth_polyA/spikeins_by_lane/"

condition <- c("Col-01_S1_", "Col-02_S1_", "Col-03_S2_", "Vir-1_S2_", "Vir-2_S2_", "Vir-3_S1_")
symbol <- c(8,8,8,8,8,8)
colour <- c('black', 'blue', 'red', 'gray', 'gray', 'gray')
isfirst <- c('TRUE', 'FALSE', 'FALSE', 'FALSE', 'FALSE', 'FALSE')
details <- data.frame(condition,symbol,colour,isfirst)
details$rep <- lapply(details$condition, function(c) load_spike_in_rep(c, path, numlanes))

plot_spikeins_by_rep(details, "ColVirSpikeins.pdf", 'Spike-ins by replicate (Col with Vir overlaid)',10000, 2e+06)

####################################################################################################
## Plot Ribo-minus spikeins by replicate (in this case directly from the replicates data)  
path <- "~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Arabidopsis_RNAMeth_ribo/spikeins/"

condition <- c("Col-1", "Col-2", "Col-3", "Col-4", "Col-5", "Col-6", "Col-7",
               "Vir-1", "Vir-2", "Vir-3", "Vir-4", "Vir-5", "Vir-6", "Vir-7")
symbol <- c(8,8,8,8,8,8,8, 21,21,21,21,21,21,21)
colour <- c('black', 'blue', 'red', 'green', 'pink', 'yellow', 'magenta')
isfirst <- c('TRUE', 'FALSE', 'FALSE', 'FALSE', 'FALSE', 'FALSE', 'FALSE',
             'FALSE', 'FALSE', 'FALSE', 'FALSE', 'FALSE', 'FALSE', 'FALSE')
details <- data.frame(condition,symbol,colour,isfirst)
details$rep <- lapply(details$condition, function(c) load_spike_in_rep_direct(c, path))

plot_spikeins_by_rep(details, "AllRiboMinus.pdf", 'Spike-ins by replicate (Ribo-minus)',10000, 2e+06)

####################################################################################################
## Analyse spike-in results across all the Arabidopsis experiments and compare to the polyA results

bpath <- "~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/"
sp <- "/spikeins/"
paths <- c("~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Arabidopsis_Exosome/spikeins/",
           "~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Arabidopsis_xrn3-8/spikeins/",
           "~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Arabidopsis_RNAMeth_ribo/spikeins/",
           "~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Arabidopsis_RNAMeth_polyA/spikeins/")
conditions <- list( c("Col-1", "Col-2", "Col-3", "Col-4", "Col-5", "Col-6", "Col-7",
                      "hen-1", "hen-2", "hen-3", "hen-4", "hen-5", "hen-6", "hen-7",
                      "mtr-1", "mtr-2", "mtr-3", "mtr-4", "mtr-5", "mtr-6", "mtr-7"),
                    c("Col-1", "Col-2", "Col-3", "xrn-1", "xrn-2", "xrn-3"),
                    c("Col-1", "Col-2", "Col-3", "Col-4", "Col-5", "Col-6", "Col-7",
                      "Vir-1", "Vir-2", "Vir-3", "Vir-4", "Vir-5", "Vir-6", "Vir-7"),
                    c("Col-01", "Col-02", "Col-03", "Vir-1", "Vir-2", "Vir-3"))

# expects counts files to be in path/condition-<fwd/rev>counts.txt

# calc experiment labels for each condition & replicate
labels <- sapply(1:length(paths), function(p) sapply(conditions[[p]], function(c) 
  rbind(sub(sp,"",sub(bpath,"",paths[p])))))

results <- calcratios(paths, conditions, labels)

plot_boxplots(results)

## plot histogram of ratios for non-polyA (but include outlying polyA point to get scale right for comparisons)
r <- as.numeric(results$ratio[results$experiment!="Arabidopsis_RNAMeth_polyA" | results$condition == "Col-03"])
plot_histogram(r, 0.012, 14, "RatioDist.pdf")

## plot with curve fitted to to non-polyA data
rcurve <- as.numeric(results$ratio[results$experiment!="Arabidopsis_RNAMeth_polyA"])
plot_histogram(r, 0.012, 14, "RatioDistWithCurve.pdf", rcurve)

## plot all data + curve fitted to non-polyA data
r <- as.numeric(results$ratio)
plot_histogram(r, 0.012, 14, "RatioDistAll.pdf", rcurve)


# calc probabilities of Col1 and Col3 results
mean=mean(as.numeric(results$ratio[results$experiment!="Arabidopsis_RNAMeth_polyA"]))
sd=sd(as.numeric(results$ratio[results$experiment!="Arabidopsis_RNAMeth_polyA"]))

probcol3 <- pnorm(0.011263, mean, sd, lower.tail=FALSE)
probcol1 <- pnorm(0.0050597, mean, sd, lower.tail=FALSE)



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

