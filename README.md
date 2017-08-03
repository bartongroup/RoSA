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

The R package depends on one third-party package, LSD, which you may need to install first, 
if it is not present already:

```
install.packages("LSD")
```

The python scripts are targetted for [python 2.7](https://www.python.org/download/releases/2.7/) and depend on:
- scipy (version 0.17.1 - will not work with latest version)
- numpy
- pandas (version 0.18 - 0.19.2 - will not work with latest version)
- pysam
- pyyaml
- six
- (optionally) drmaa

The python script to find and count spliced antisense and sense reads also depends on [sambamba](http://lomereiter.github.io/sambamba/) being installed.

## Installation

## How to use RoSA

RoSA (currently) is a combination of an R package and some python scripts. The R package takes as input datasets containing several different read counts:

- Full read counts by gene
- Antisense counts by gene (via RoSA's python script to create an antisense annotation, and then read counting as usual)
At least one of:
- Spike-in sense and antisense counts (by aligning the reads data to the spike-ins, and generating read counts for each strand)
- Spliced sense and antisense counts (via RoSA's python script to filter “sense-strand” spliced reads)

The python scripts can be used to 
* create an antisense annotation (as gtf) from a standard annotation (as gff or gtf), which can then be used to generate antisense read counts via your favourite read counting tool (e.g. featureCounts):
```
from viewseq.gffParser import GffParser
p = GffParser("inputannotation.gff")
p.build_antisense_gtf_gene_only([(p.featureType.exon,"outfile.gtf")])
```
or
```
from viewseq.gtfParser import GtfParser
p = GtfParser("inputannotation.gtf")
p.build_antisense_gtf_gene_only([(p.featureType.exon,"outfile.gtf")])
```
* generate sense and antisense counts of reads at splice junctions. The script takes the usual annotation (as gtf/gff) and corresponding alignment (as bam) and outputs counts of spliced sense and antisense reads to a designated output file. Because this script must process an entire bam file of reads, this is very slow.

## Known issues

Installing pysam via pip can result in a [broken installation](https://github.com/pysam-developers/pysam/issues/475)
Workaround: use bioconda instead.

## Contact information

The `RoSA` R package was developed within [The Barton Group](http://www.compbio.dundee.ac.uk) at [The University of Dundee](http://www.dundee.ac.uk)
by Dr. Kira Mourão and Dr. Nick Schurch.

To **report problems** or ask for **assistance**, please raise a new issue [on the project's support forum](https://github.com/bartongroup/RoSA/issues).
Providing a *reproducible working example* that demonstrates your issue is strongly encouraged.  Please also browse/search
the support forum before posting a new issue, in case your question is already answered there.
