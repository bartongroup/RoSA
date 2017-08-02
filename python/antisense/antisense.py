#!/usr/bin/env python

""" Antisense functions """

import datetime
import argparse
import logging
import pandas as pd
import pandas.io.common
import os
import tempfile
import shutil

from summary_routines.runcommand import RunCommand
from summary_routines.gtfparser import GtfParser
from summary_routines.constants import *
from anti_reads_as_job import AntisenseJob

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"

# constants
GENE_ID = 'gene_id'
EXON_COUNTS = 'exoncounts'
INTRON_COUNTS = 'introncounts'
EXON_SUM = 'exonsum'


class Antisense:

    def __init__(self):
        self.__logger = logging.getLogger(__name__)

        # constants defined by featureCounts
        self.__start_col = "Start"  # name of Start column in featureCounts output
        self.__stop_col = "End"    # name of Stop column in featureCounts output
        self.__geneid = "Geneid"   # name of Gene id column in featureCounts output

    def id_structured_antisense(self, exondata, introndata, output, k=0.5, min_reads=100):
        """ Find genes where the antisense expression is structured as the sense expression is
        Note, limited to genes which have at least one intron
        :param exondata: ftrcounts output containing exon reads counts data
        :param introndata: ftrcounts output containing intron reads counts data
        :param output: path and name of output file
        :param k: fraction of exon counts which intron counts must be less than (defaults to 0.5)
        :param min_reads: min number of reads on antisense strand to consider the gene (defaults to 100)
        :return: list of genes
        """

        self.__logger.info("Identifying structured antisense genes from featureCounts data")

        exons = pd.read_csv(exondata, header=1, sep='\t', engine='python')
        introns = pd.read_csv(introndata, header=1, sep='\t', engine='python')
        exons = exons.rename(columns={self.__geneid: GENE_ID})
        introns = introns.rename(columns={self.__geneid: GENE_ID})

        # get reads columns
        cols = [col for col in exons.columns if ".bam" in col]
        exons[EXON_COUNTS] = exons[cols].sum(axis=1)

        cols = [col for col in introns.columns if ".bam" in col]
        introns[INTRON_COUNTS] = introns[cols].sum(axis=1)

        self.__logger.info("Calculate mean exon reads and mean intron reads per gene")
        # need to normalise by exon length
        exons["length"] = abs(exons["Start"] - exons["End"])
        exons["ex_norm_counts"] = exons[EXON_COUNTS] / exons["length"]
        exonmean = exons[[GENE_ID, "ex_norm_counts"]].groupby(GENE_ID).mean()
        exonsum = exons[[GENE_ID, EXON_COUNTS]].groupby(GENE_ID).sum()
        exonsum = exonsum.rename(columns={EXON_COUNTS: EXON_SUM})

        introns["length"] = abs(introns["Start"] - introns["End"])
        introns["in_norm_counts"] = introns[INTRON_COUNTS] / introns["length"]
        intronmean = introns[[GENE_ID, "in_norm_counts"]].groupby(GENE_ID).mean()

        # merge the results by gene
        exonmean = exonmean.join(exonsum)
        combined = exonmean.join(intronmean, how="inner")
        combined["in_norm_counts"] = combined["in_norm_counts"].fillna(0)

        combined.to_csv(output)
        # filter for only genes where intronmean < k x exonmean
        result = combined[(combined["in_norm_counts"] < k * combined["ex_norm_counts"]) & (combined[EXON_SUM] > min_reads)]

        #result.to_csv(output)

        return result

    def get_diff_strand_overlaps(self, overlapping_file, annotation_file):

        # given list of overlapping genes, find genes which overlap on same strand and export
        overlaps = pd.read_csv(overlapping_file, sep='\t', engine='python')
        a = GtfParser(annotation_file)

        annotation = a.featuredict['exon'].groupby(GENE_ID).first().reset_index()
        overlaps.columns = [GENE_ID, "overlap", "gene2"]
        overlaps = pd.merge(overlaps, annotation[[GENE_ID, STRAND]], on=GENE_ID)
        overlaps.columns = ["gene1", "overlap", GENE_ID, "strand1"]
        overlaps = pd.merge(overlaps, annotation[[GENE_ID, STRAND]], on=GENE_ID)
        overlaps.columns = ["gene1", "overlap", "gene2", "strand1", "strand2"]
        overlaps = overlaps[overlaps["strand1"] != overlaps["strand2"]]
        result = pd.concat([overlaps["gene1"], overlaps["gene2"]])
        result = result.drop_duplicates()

        return result

    def filter_gtf(self, filter_file, gtf_file):

        # output gtf file which only has genes in filter_file in it
        filterlist = pd.read_csv(filter_file, sep='\t', engine='python')
        to_filter = pd.read_csv(gtf_file, sep='\t', engine='python')

        # extract all rows from to_filter which contain an entry in filter
        # then output as gtf
        filterlist.columns = [PARENT]

        to_filter.columns = [CHROMOSOME, SOURCE, TYPE, START, STOP, SCORE, STRAND, PHASE, ATTRIBUTES]
        to_filter[PARENT] = to_filter[ATTRIBUTES].str.extract('[Gg]ene_id \"([^\"]*)', expand=False)

        filtered = pd.merge(to_filter, filterlist, on=PARENT)

        filtered.drop([PARENT], inplace=True, axis=1)

        return filtered

    def find_anti_reads_in_sense(self, anti_file, sense_file):

        anti = pd.read_csv(anti_file, r"\s*", engine="python")
        anti.columns = ["read_id", "assigned", "gene_id", "details"]
        sense_chunks = pd.read_csv(sense_file, r"\s*", iterator=True, chunksize=1000, engine="python")

        duplicates = []
        for chunk in sense_chunks:
            chunk.columns = ["read_id", "assigned", "gene_id", "details"]
            duplicates.append(chunk[chunk["read_id"].isin(anti["read_id"])])

        all_duplicates = pd.concat(duplicates)

        return all_duplicates

    def find_anti_with_sense_structure(self, annotation_file, alignment_file, outfile, intron_file=""):

        logger = logging.getLogger(__name__)

        counts = pd.DataFrame()
        r = RunCommand()

        if intron_file == "":
            logger.info("Calculating introns")
            intronlist = self._get_introns(annotation_file)

            hometempdir = "{}/temp_anti".format(os.environ['HOME'])
            try:
                os.mkdir(hometempdir)
            except OSError:
                pass

            intron_file = os.path.join(hometempdir, "intronfile.csv")
            intronlist.to_csv(path_or_buf=intron_file, sep="\t", index=False)
        else:
            intronlist = pd.read_csv(intron_file, header=0, sep="\t")

        lastchro = -1
        poslist = intronlist[["Chromosome", "Stop"]].groupby("Chromosome").max()

        logger.info("Building cmds")
        cmds = []
        fnames = []
        for (chro, maxpos) in poslist.to_records(index=True):

            logger.info("Chromosome {}".format(chro))
            logger.info("Max pos {}".format(maxpos))
            step = 100000
            for lastpos in range(1, maxpos, step):  # slice out a range of base pair positions

                if chro != lastchro:
                    lastchro = chro

                pos = lastpos + step - 1
                cmd = ["python", "../anti_reads_as_job.py", "-l", alignment_file, "-i", intron_file, "-s",
                       str(lastpos), "-e", str(pos), "-o", chro]
                filename = "anti-{}.{}-{}".format(chro, lastpos, pos)

                self.__logger.info("Appending command: {}".format(cmd))
                cmds.append(cmd)
                fnames.append(filename)

        self.__logger.info("Running command list")
        result_files = r.run_commands(cmds, fnames, output_file=True)

        for result_file in result_files:

            try:
                result = pd.read_csv(result_file)

                if counts.empty:
                    counts = result
                else:
                    counts = counts.append(result, ignore_index=True)

            except pandas.io.common.EmptyDataError as e:
                # no data for read_chunks at all
                self.__logger.info("No data reading in results from {}: {}".format(result_file, e))
                pass

        logger.info("Completed processing reads")

        if not counts.empty:
            # aggregate counts in case we have multiple cases of same gene: need to aggregate into single count
            counts = counts.groupby("gene_id")["anticounts", "sensecounts"].sum()

            counts.fillna(0, inplace=True)

            # identify genes with intron structure on opposite strand
            counts["ratio"] = counts["anticounts"] / counts["sensecounts"]
            counts.sort_values("ratio", inplace=True)

        # write counts out to a temporary directory
        logger.info("Writing results to temp directory")
        tempdir = tempfile.gettempdir()
        loc = os.path.join(tempdir, "counts.csv")
        counts.to_csv(path_or_buf=loc, sep=',')

        # finally move the output to specified location
        logger.info("Moving output to {}".format(outfile))
        shutil.move(loc, outfile)

    def _get_introns(self, annotation_file):

        a = GtfParser(annotation_file)
        intronlist = a.export_to_gtf2("","intron", out_to_file=False)
        intronlist = intronlist[["Chromosome", "Strand", "Start", "Stop", "gene_id"]]

        return intronlist

def main():

    logname = "antisense-log-{}.txt".format(datetime.datetime.now().isoformat())
    logging.basicConfig(filename=logname, filemode='w', level=logging.DEBUG,
                        format='%(asctime)s %(message)s')  # include timestamp

    # parse arguments, which should just be a config file
    argparser = argparse.ArgumentParser(description="A script to calculate antisense counts")
    required_named = argparser.add_argument_group('required named arguments')
    required_named.add_argument('-a', action='store', dest='annot_file',
                                help="Specify the location of the annotation file.", required=True)

    required_named.add_argument('-l', action='store', dest='align_file',
                                help="Specify the location of the alignment file.", required=True)

    required_named.add_argument('-o', action='store', dest='output_file',
                                help="Specify the name of the output file.", required=True)

    argparser.add_argument('-i', action='store', dest='input_file',
                                help="Specify the name of the input file.")

    args = argparser.parse_args()

    a = Antisense()
    if args.input_file:
        a.find_anti_with_sense_structure(args.annot_file, args.align_file, args.output_file, args.input_file)
    else:
        a.find_anti_with_sense_structure(args.annot_file, args.align_file, args.output_file)

if __name__ == "__main__":
    main()
