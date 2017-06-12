# RoSA: a tool for the Removal of Spurious Antisense

In stranded RNA-Seq experiments we have the opportunity to detect and measure antisense transcription, important since antisense transcripts impact gene transcription in several different ways. Stranded RNA-Seq determines the strand from which an RNA fragment originates, and so can be used to identify where antisense transcription may be implicated in gene regulation. 

However, spurious antisense reads are often present in experiments, and can manifest at levels greater than 1% of sense transcript levels. This is enough to disrupt analyses by causing false antisense counts to dominate the set of genes with high antisense transcription levels.   

The RoSA (Removal of Spurious Antisense) tool detects the presence of high levels of spurious antisense transcripts, by analysing ERCC spike-in data to find the ratio of antisense:sense transcripts in the spike-ins. Similarly, RoSA will calculate a correction to the antisense counts based on either the spike-in antisense:sense ratio, or, where possible, using antisense and sense counts around splice sites to provide a gene-specific correction. 

# Pre-requisites to running RoSA
## Data
RoSA depends on having a stranded RNA-Seq dataset with spike-ins. Several different read counts are needed:

- Full read counts
- Antisense counts (via a python script to create an antisense annotation)
- Spike-in counts (by aligning the reads data to the spike-ins, and generating read counts for each strand)
- Spliced read counts (via a python script to filter “sense-strand” spliced reads)

Name each read counts file according to its strand and library as follows: 
\<library name\>fwdcounts.txt and \<library name>revcounts.txt.
The counts files should contain the following tab-delimited columns:
- Geneid
- Chromosome
- Start
- End
- Strand
- Length
- Forward counts
- Reverse counts

Each first row of the file should be a header row, with column names (the names themselves do not matter). Note that if featureCounts is used, it outputs an additional line at the top of the file which must be removed before inputting to RoSA.

## Software
RoSA is an R package with no additional dependencies.

# Installation

# Running RoSA
