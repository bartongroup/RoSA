# ViewSeq
ViewSeq is a tool to summarise and visualise RNA-Seq data. It gives a swift high-level overview of the main characteristics of RNA-Seq datasets. The tool automatically runs read summarisation tools, and presents the results via an interactive graphical web interface, enabling quick identification of potential genes of interest, and including a visual overview of the alignment quality. It aims to help with questions such as:

•	How good is the quality of the alignment?

•	Which genes are most highly expressed, or not expressed at all?

•	What are the ratios of the read counts in introns / untranslated regions (UTRs) / antisense expression?

•	What is the range of exon lengths and numbers of exons per gene?

## Installation

ViewSeq can be installed directly from this GitHub repository via pip:  
`pip install git+https://github.com/bartongroup/km-ViewSeq.git`

## Pre-requisites

[featureCounts](http://subread.sourceforge.net) should be installed so that ViewSeq can perform counts of reads assigned to different features. If featureCounts is not installed, ViewSeq will still run, but only annotation statistics will be calculated and displayed.

## Quick start

Download the example data from the example directory in the repository. Here you will find an annotation file (*annotation.gtf*), a corresponding alignment file (*example.bam*) and a configuration file (*config.yaml*). The data comes from an RNA-Seq experiment on *Arabidopsis thaliana*, but is restricted to Chromosome 1, due to file size limitations on GitHub. The configuration file assumes that featureCounts is in your PATH, otherwise you should edit the *featurecounts_exe* line in *config.yaml* to give the location of the featureCounts executable.

In the example directory, run `viewseq config.yaml`

ViewSeq will output timestamped progress updates as it processes the annotation and reads files:

```
14:57:21 Initialising summary  
14:57:21 Adding annotation file: annotation.gtf  
14:57:21 Adding read file: example.bam  
14:57:21 Building Summary  
14:57:21 Calculating statistics...  
14:57:21 Reading in annotation file annotation.gtf  
14:57:23 Extracting gtf attributes to columns  
14:57:24 Checking types exist in gff file: exon, CDS  
14:57:24 Checking types exist in gff file: three_prime_UTR, five_prime_UTR, CDS  
14:57:24 Extracting regions for different features  
14:57:26 Processing alignments: example.bam  
14:57:26 Checking chromosome names in annotation and reads files match  
14:57:27 Running featurecounts 
  < featureCounts output >
14:57:30 Reading featurecounts output  
14:57:32 Running mismatch counts  
14:57:32 Counting mismatches in reads file: example.bam  
14:57:32 No mismatch tag (NM) was found in the reads file example.bam  
14:57:32 Mismatch counting completed  
14:57:32 Outputting visualisation data  
14:57:32 Outputting first exons JSON file  
14:57:32 Outputting last exons JSON file  
14:57:32 Outputting annotation features JSON files  
14:57:32 Outputting totals JSON file  
14:57:32 Outputting metadata JSON file  
14:57:32 Outputting matches JSON file  
14:57:32 Scraping STAR logs  
14:57:32 Finished scraping STAR logs  
14:57:32 Summary complete. Results can be browsed in output/summary.html.
```

Open output/summary.html in a web browser to see the results.


## The configuration file *config.yaml*
