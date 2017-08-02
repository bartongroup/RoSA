#!/usr/bin/env python

""" An array job which counts antisense reads """

import datetime
import argparse
import logging
import pandas as pd
import pandas.io.common
import numpy as np
import sys

from summary_routines.runcommand import RunCommand

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"

# constants
GENE_ID = 'gene_id'
EXON_COUNTS = 'exoncounts'
INTRON_COUNTS = 'introncounts'
EXON_SUM = 'exonsum'

class AntisenseJob:

    def __init__(self):
        self.__logger = logging.getLogger(__name__)

        # constants defined by featureCounts
        self.__start_col = "Start"  # name of Start column in featureCounts output
        self.__stop_col = "End"    # name of Stop column in featureCounts output
        self.__geneid = "Geneid"   # name of Gene id column in featureCounts output

    def calc_anti_for_slice(self, alignment_file, intron_file, chro, lastpos, pos):

        counts = pd.DataFrame()
        self.__logger.info("Slice: {} - {}".format(lastpos, pos))
        self.__logger.info("Chromosome: {}".format(chro))
        self.__logger.info("Intron file: {}".format(intron_file))
        self.__logger.info("Alignment file: {}".format(alignment_file))

        intronlist = pd.read_csv(intron_file, header=0, sep="\t")
        intronlist["Chromosome"] = intronlist["Chromosome"].astype(str)

        try:
            read_chunks = self._get_reads(alignment_file, chro, lastpos, pos)  # chunk reads in the selected range

            for reads_df in read_chunks:

                reads_df.columns = ["qname", "flag", "chromosome", "pos", "cigar"]

                # convert chromosome to number so that we can compare with intron list chromosomes
                reads_df["chromosome"] = reads_df["chromosome"].astype(str)

                if not reads_df.empty:

                    reads_df = self._add_strand_col(reads_df)
                    reads_df = self._create_splice_points(reads_df)

                    if counts.empty:
                        counts = self._count_reads_by_sense(intronlist, reads_df)
                    else:
                        counts = counts.append(self._count_reads_by_sense(intronlist, reads_df), ignore_index=True)
        except IOError as e:
            # non-existent file if there are no reads in this region
            self.__logger.info("Error at read_chunks: {}".format(e))
            pass
        except pandas.io.common.EmptyDataError as e:
            # no data for read_chunks at all
            self.__logger.info("Error at read_chunks: {}".format(e))
            pass

        counts.to_csv(sys.stdout, index=False)

    def _get_reads(self, alignment_file, chro, last_pos, position):

        r = RunCommand()

        # because we are expanding into a sam file and then reading, we don't want to create giant files...
        # so filter out small bits of the bam file at a time
        # use position rather than a a region to select reads, so we don't get the same read twice
        # don't put -F argument in quotes if using serial call - will break subprocess
        cmd = ["sambamba", "view", "-F",
               "position > {} and position <= {} and cigar =~ /[0-9A-Z]*N[0-9A-Z]*/ and proper_pair and mapping_quality>50 and [NH]==1".format(last_pos, position),
               alignment_file, "{}:{}-{}".format(chro, last_pos, position)]

        self.__logger.info("Running command: ".format(cmd))
        result = r.run_commands([cmd], output_file=True, forceSerial=True)
        # for highIO we need position in cmd string in quotes, o/w should be without quotes

        reads_chunks = pd.read_csv(result[0], sep="\\t", usecols=[0, 1, 2, 3, 5],
                                   header=None, iterator=True, chunksize=10000, engine="python")

        return reads_chunks

    def _add_strand_col(self, reads_df):

        self.__logger.info("Adding strand column")

        # add a strand column
        # 0x10 - seq on reverse strand + 0x40 - 1st in pair : +
        # or 0x20 - mate on reverse strand + 0x80 - 2nd in pair: +
        # 0x10 - seq on reverse strand + 0x80 - 2nd in pair: -
        # or 0x20 - mate on reverse strand + 0x40 - 2nd in pair: -
        # require 0x02 and 0x01 as well, then we can filter out invalid reads
        reads_df["strand"] = ""
        reads_df["strand"] = np.where((((reads_df["flag"] & 0x10) > 0) & ((reads_df["flag"] & 0x40) > 0)) |
                                      (((reads_df["flag"] & 0x20) > 0) & ((reads_df["flag"] & 0x80) > 0)), "+", "")

        reads_df["strand"] = np.where((((reads_df["flag"] & 0x10) > 0) & ((reads_df["flag"] & 0x80) > 0)) |
                                      (((reads_df["flag"] & 0x20) > 0) & ((reads_df["flag"] & 0x40) > 0)), "-",
                                      reads_df["strand"])

        reads_df = reads_df[reads_df["strand"] != ""]

        return reads_df

    def _create_splice_points(self, reads_df):

        self.__logger.info("Creating splice points")

        # make a row for each splice point
        # add a column for the splice position
        reads_df["split_cigar"] = reads_df["cigar"].str.findall(r'([0-9]+)([MIDNSHPX=])')

        # split out into columns and then pivot into rows
        s = reads_df["split_cigar"].apply(pd.Series).stack()

        # parse the cigar string by index to sum all values up to this row's position
        # so if index = 1 and cigar = 23M41N10M then we want 23 + 41
        s = s.apply(pd.Series)
        s.columns = ["cigar_value", "cigar_type"]
        s["cigar_value"] = pd.to_numeric(s["cigar_value"])  # make numeric for summing

        s["cumsums"] = s.groupby(level=[0]).cumsum()["cigar_value"]

        s.index = s.index.droplevel(-1)  # to line up with df's index
        s.name = "split_cigar"  # s needs a name in order to join with df
        del reads_df["split_cigar"]  # delete the current parent column
        reads_df = reads_df.join(s)  # join the new parent column
        reads_df.reset_index(inplace=True)

        # keep only cigar type N, these are where the splices occur
        reads_df = reads_df[reads_df["cigar_type"] == "N"]
        reads_df["end_intron"] = reads_df["pos"] + reads_df["cumsums"] - 1
        reads_df["start_intron"] = reads_df["end_intron"] - reads_df["cigar_value"] + 1

        return reads_df[["index", "qname", "chromosome", "cigar", "strand", "start_intron", "end_intron"]]

    def _count_reads_by_sense(self, intronlist, reads_df):

        self.__logger.info("Counting spliced reads on sense and antisense strands")

        # for each intron find reads which match splicing on opposite strand
        introns_and_reads = pd.merge(reads_df, intronlist, how="inner",
                                     left_on=["start_intron", "end_intron", "chromosome"],
                                     right_on=["Start", "Stop", "Chromosome"])

        # extract the read ids
        # filter those reads out of the main bam file, and then call featureCounts on the filtered bam

        anti_i_r = introns_and_reads[introns_and_reads["Strand"] != introns_and_reads["strand"]]
        sense_i_r = introns_and_reads[introns_and_reads["Strand"] == introns_and_reads["strand"]]

        # now count number of spliced reads on sense and antisense strand, per gene
        sensecounts = sense_i_r[["gene_id", "index"]].groupby("index").first().reset_index().groupby("gene_id").count()
        anticounts = anti_i_r[["gene_id", "index"]].groupby("index").first().reset_index().groupby("gene_id").count()

        try:
            counts = pd.merge(sensecounts, anticounts, how='outer', left_index=True, right_index=True)
            counts = counts.rename(columns={"index_x": "sensecounts", "index_y": "anticounts"})

        except KeyError as e:
            # one (or both) dataframe(s) is empty
            if anticounts.empty:
                # stop now, no antisense to deal with!
                e.message = "No antisense counts present in data"
                raise e

        return counts.reset_index()


def main():

    logname = "antisense-log-{}.txt".format(datetime.datetime.now().isoformat())
    logging.basicConfig(filename=logname, filemode='w', level=logging.DEBUG,
                        format='%(asctime)s %(message)s')  # include timestamp

    # parse arguments, which should just be a config file
    argparser = argparse.ArgumentParser(description="A script to calculate antisense counts")
    required_named = argparser.add_argument_group('required named arguments')
    required_named.add_argument('-i', action='store', dest='intron_file',
                                help="Specify the location of the intron file.", required=True)

    required_named.add_argument('-l', action='store', dest='align_file',
                                help="Specify the location of the alignment file.", required=True)

    required_named.add_argument('-s', action='store', dest='start',
                                help="Specify the start position.", required=True)

    required_named.add_argument('-e', action='store', dest='end',
                                help="Specify the end position.", required=True)

    required_named.add_argument('-o', action='store', dest='chromosome',
                                help="Specify the chromosome.", required=True)

    args = argparser.parse_args()

    a = AntisenseJob()
    a.calc_anti_for_slice(args.align_file, args.intron_file, args.chromosome, args.start, args.end)


if __name__ == "__main__":
    main()
