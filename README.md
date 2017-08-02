# RoSA: a tool for the Removal of Spurious Antisense

In stranded RNA-Seq experiments we have the opportunity to detect and measure antisense transcription, important since antisense transcripts impact gene transcription in several different ways. Stranded RNA-Seq determines the strand from which an RNA fragment originates, and so can be used to identify where antisense transcription may be implicated in gene regulation. 

However, spurious antisense reads are often present in experiments, and can manifest at levels greater than 1% of sense transcript levels. This is enough to disrupt analyses by causing false antisense counts to dominate the set of genes with high antisense transcription levels.   

The RoSA (Removal of Spurious Antisense) tool detects the presence of high levels of spurious antisense transcripts, by:
* analysing ERCC spike-in data to find the ratio of antisense:sense transcripts in the spike-ins; or
* using antisense and sense counts around splice sites to provide a set of gene-specific estimates; or
* both.

Once RoSA has an estimate of the spurious antisense, expressed as a ratio of antisense:sense counts, RoSA will calculate a correction to the antisense counts based on the ratio. Where a gene-specific estimate is available for a gene, it will be used in preference to the global estimate obtained from either spike-ins or spliced reads.

## Pre-requisites to running RoSA
### Data
RoSA depends on having a stranded RNA-Seq dataset. Ideally the data will also have spikeins, but these are not essential as RoSA can also estimate the spurious antisense in your data from looking at spliced reads at known splice junctions.

### Annotation
The data should have already been aligned to an annotation. RoSA will need the gff or gtf file for the annotation, in order to create an antisense annotation which can then be used to generate antisense counts, and to determine where known splice junctions occur.

### Dependencies

The package depends on one third-party package, LSD, which you may need to install first, 
if it is not present already:

```
install.packages("LSD")
```

## Installation

## How to use RoSA

RoSA (currently) is a combination of an R package and some python scripts. The R package takes as input datasets containing several different read counts:

- Full read counts by gene
- Antisense counts by gene (via a python script to create an antisense annotation)
At least one of:
- Spike-in sense and antisense counts (by aligning the reads data to the spike-ins, and generating read counts for each strand)
- Spliced sense and antisense counts (via a python script to filter “sense-strand” spliced reads)

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

## Contact information

The `RoSA` R package was developed within [The Barton Group](http://www.compbio.dundee.ac.uk) at [The University of Dundee](http://www.dundee.ac.uk)
by Dr. Kira Mourão and Dr. Nick Schurch.

To **report problems** or ask for **assistance**, please raise a new issue [on the project's support forum](https://github.com/bartongroup/RoSA/issues).
Providing a *reproducible working example* that demonstrates your issue is strongly encouraged.  Please also browse/search
the support forum before posting a new issue, in case your question is already answered there.
