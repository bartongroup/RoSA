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

RoSA (currently) is a combination of an R package and some python scripts for preprocessing.
The R package depends on one third-party package, LSD, which you may need to install first, 
if it is not present already:

```
install.packages("LSD")
```

The python scripts are targetted for [python 2.7](https://www.python.org/download/releases/2.7/) and depend on:
- scipy (version 0.16.1 - 0.17.1 - see known issues)
- numpy
- pandas (not 0.20.1 - see known issues)
- six
- (optionally) drmaa

The python script to find and count spliced antisense and sense reads also depends on [sambamba](http://lomereiter.github.io/sambamba/) being installed.

## Installation

To install the R rosa package from github, while RoSA is still a private repository, use `devtools::install_github`. You will need to install devtools if you have not already done so, and to obtain a personal access token through your github account.
```
devtools::install_github("bartongroup/km-rosa", auth_token="<your PAT>")
```
Load the package to get access to RoSA:
```
library("rosa", lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")
```

RoSA's python scripts are set up in a python package which can be installed via:
```
pip install git+https://github.com/bartongroup/km-rosa.git#subdirectory=python
```
and uninstalled via:
```
pip uninstall rosa
```
The pip installer will handle installing the python dependencies with the correct versions.


## How to use RoSA

The R package takes as input datasets containing several different read counts:

1. Full read counts by gene
2. Antisense counts by gene (via RoSA's python script (*makeannotation.py*) to create an antisense annotation, and then read counting as usual)
3. At least one of:
     1. Spike-in sense and antisense counts (by aligning the reads data to the spike-ins, and generating read counts for each strand)
     2. Spliced sense and antisense counts (via RoSA's python script (*antisense.py*) to filter spliced reads which occur at known splice junctions)

Help for the R rosa functionality can be found by typing `help(rosa)` after installing and loading the RoSA R package. Since some of RoSA's inputs are non-standard, the python preprocessing scripts are supplied to make it easier to generate the inputs required by the package (specifically, input (2) and input (3(i))).

The *make_annotation* script creates an antisense annotation (as gtf) from a standard annotation (as gff or gtf), which can then be used to generate antisense read counts (input 2) via your favourite read counting tool (e.g. [featureCounts](http://subread.sourceforge.net)):
```
make_annotation -a <annotation file as gff or gtf> -o <output filename without file extension>
```
The annotation produced by *make_annotation* only contains antisense features and so cannot be used in place of a standard annotation. The *make_annotation* script can also be called from R:
```
make_annotation(<annotation file>, <output filename>)
```

The *count_spliced* script generates sense and antisense counts of reads at splice junctions (input 3(i)). The script takes a standard annotation (as gtf/gff) and corresponding alignment (as bam) and outputs counts of spliced sense and antisense reads to a designated output file. An index file (.bai file) should also have been pre-generated and be in the same directory as the bam file. Because the script must process an entire bam file of reads, it is very slow. The script is set up to break the bam file into chunks and process each chunk separately using sambamba and some custom filtering. On a cluster with drmaa installed, the script will use drmaa to submit each chunk as a separate job. On a single machine, the script will spawn a new process to run each chunk separately. Once all of the jobs have run, the script collates the results to give a count of the spliced reads. The script is run as follows:
```
count_spliced -a <annotation file as gff or gtf> -l <alignment file as bam> -o <output file>
```
An additional -i option allows the user to input a file containing the extents of all introns, which is otherwise calculated as part of splice counting process. This is primarily useful for testing on the same annotation multiple times, as the intron file can be calculated once and re-used, saving some time. The location of the generated intron file is output by the script.

The *count_spliced* script can also be called from R:
```
count_spliced(<annotation file>,<alignment file>,<output filename>)
```

Both scripts output timestamped logfiles recording progress.

## Known issues

* RoSA is incompatible with pandas version 0.20.1, due to [bug #16921](https://github.com/pandas-dev/pandas/issues/16921)
* "command not found" errors when calling the preprocessing scripts from RStudio can be caused on Macs by a bug where RStudio uses an incorrect PATH variable. Running RStudio from the command line avoids the problem: 
`open /Applications/RStudio.app`

## Contact information

The `RoSA` R and python tool was developed within [The Barton Group](http://www.compbio.dundee.ac.uk) at [The University of Dundee](http://www.dundee.ac.uk)
by Dr. Kira Mourão and Dr. Nick Schurch.

To **report problems** or ask for **assistance**, please raise a new issue [on the project's support forum](https://github.com/bartongroup/RoSA/issues).
Providing a *reproducible working example* that demonstrates your issue is strongly encouraged.  Please also browse/search
the support forum before posting a new issue, in case your question is already answered there.
