setup <- function()
{
# load read counts for genes, on antisense and sense strands - from featurecounts
Col.anti.counts.genes <- read.table("~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Col-anti-counts-genes.out", header=TRUE, quote="\"")
Col.sense.counts.genes <- read.table("~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Col-sense-counts-genes.out", header=TRUE, quote="\"")

# rename dataframes and columns
anti <- Col.anti.counts.genes
rm(Col.anti.counts.genes)
sense <- Col.sense.counts.genes
rm(Col.sense.counts.genes)

colnames(anti)[7] <- "fwd"
colnames(anti)[8] <- "rev"
colnames(sense)[7] <- "fwd"
colnames(sense)[8] <- "rev"

# calculate raw ratios
ratios <- anti[,c("Geneid","fwd","rev")]
ratios <- merge(x = ratios, y = sense, by = "Geneid", all = TRUE)
ratios <- ratios[,c("Geneid","fwd.x", "fwd.y", "rev.x", "rev.y")]
colnames(ratios)[2] <- "fwd.anti"
colnames(ratios)[3] <- "fwd.sense"
colnames(ratios)[4] <- "rev.anti"
colnames(ratios)[5] <- "rev.sense"

ratios["fwdratio"] <- ratios["rev.anti"]/ratios["fwd.sense"]
ratios["revratio"] <- ratios["fwd.anti"]/ratios["rev.sense"]

# draw histograms of raw ratios
hist(log(ratios["fwdratio"][sapply(ratios["fwdratio"], is.finite)]), xlim=c(-10,5),xlab="Log(fwd antisense:sense ratio)", 
     ylab="Number of genes", col="blue", main="All genes antisense ratios: forward strand")
hist(log(ratios["revratio"][sapply(ratios["revratio"], is.finite)]), xlim=c(-10,5),xlab="Log(rev antisense:sense ratio)", 
     ylab="Number of genes", col="blue", main="All genes antisense ratios: reverse strand")

# load list of genes identified as having antisense-sense structure (actually now just have full dataset, 
# filters below take out the genes with specific structure)
result <- read.csv("~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/result.csv")
colnames(result)[1] <- "Geneid"

# inner join with ratios to get ratios for only genes which have same sense-antisense structure (in theory)
ratios2 <- merge(ratios,result["Geneid"],by="Geneid")

# draw histograms of these ratios
hist(log(ratios2["fwdratio"][sapply(ratios2["fwdratio"], is.finite)]), xlim=c(-10,5),xlab="Log(fwd antisense:sense ratio)", ylab="Number of genes", col="blue", main="Genes with reads>100, k=0.2: forward strand")
hist(log(ratios2["revratio"][sapply(ratios2["revratio"], is.finite)]), xlim=c(-10,5),xlab="Log(rev antisense:sense ratio)", ylab="Number of genes", col="blue", main="Genes selected by k=1: Antisense ratios: reverse strand")

# various filters on the ratios
filter <- result[["in_norm_counts"]] < 0.2 * result[["ex_norm_counts"]] & result[["exonsum"]] > 100
filter <- result[["in_norm_counts"]] < 0.2 * result[["ex_norm_counts"]]
filter <- result[["exonsum"]] > 10
result2 <- result[filter,]
ratios2 <- merge(ratios,result2["Geneid"],by="Geneid")

}

spikeins <- function()
{
# read in spike ins data and calculate ratios - plot histogram  
spikeins <- read.csv("~/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/spikein_counts.csv")
spike_ratios <- spikeins[,c("Geneid","fwd.sense", "fwd.anti", "rev.sense", "rev.anti")]
spike_ratios["fwdratio"] <- spike_ratios["rev.anti"]/spike_ratios["fwd.sense"]
hist(log(spike_ratios["fwdratio"][sapply(spike_ratios["fwdratio"], is.finite)]), xlim=c(-10,5), xlab="Log(antisense:sense ratio)", ylab="Number of genes", col="blue", main="Spike-in antisense ratios")

}

# getting annotation data
t <- read.table("~/Documents/Arabidopsis_RNAMeth_polyA/genes_fixed.gtf", sep = "\t", quote="\"")
# extract geneids
t$gene_id <- sapply(strsplit(sapply(strsplit(as.character(t$V9), "gene_id "), `[`, 2), ";"), "[", 1)
# group by geneid to get strand per geneid and inner join to get strands onto ratios2 table
t2 <- unique(t[c("V7", "gene_id")])
colnames(t2) = c("strand","Geneid")
ratios2 <- merge(ratios2, t2, by="Geneid")

outlierpos <- ratios2[ratios2$fwdratio > 1,] # but this produces far fewer than in earlier ofinterest calcs??

# take ratios before merge with results to do this calc
outlierpos <- ratios[ratios$fwdratio > 1 & ratios$rev.anti > 100,]
outlierpos <- outlierpos[order(-outlierpos$fwdratio),]


# to get full set of what we want:
extranti <- ratios2[ratios2$fwdratio > 1 | ratios2$revratio > 1,]
extranti<- extranti[complete.cases(extranti[,8]),] # get rid of nans
wrongfstrand <- extranti[extranti$fwdratio > 1 & extranti$strand=="+",]
wrongrstrand <- extranti[extranti$revratio > 1 & extranti$strand=="-",]

