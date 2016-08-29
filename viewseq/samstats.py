#!/usr/bin/env python

import logging
import pysam
import locale

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"


class SamStats:
    def __init__(self, reads_filenames):
        """ Get read counts for the reads files.
        :param reads_filenames: list of reads files
        """

        self.__logger = logging.getLogger(__name__)

        locale.setlocale(locale.LC_ALL, 'en_GB.UTF-8')

        self.options = ["Total reads",
                        "0 mismatches",
                        "1 mismatch",
                        "2 mismatches",
                        "3 mismatches"]

        # run counts
        self.__logger.info("Running mismatch counts")
        results = self._get_mismatches(reads_filenames)

        self.stats = {}

        # munge the result into the expected tuples, with formatted count values
        for strand, totals in results.items():
            self.stats[strand] = []
            format_totals = [locale.format("%d", t, grouping=True) for t in totals]  # format locale of each total
            self.stats[strand].extend(zip(self.options, format_totals))

    def _get_mismatches(self, reads_filenames):
        """
        Get mismatch counts from each reads file
        :param reads_filenames: list of reads files
        :return: counts of reads with 0, 1, 2 and 3 mismatches, and total number of reads, by strand
        """

        fwdtotal = 0                            # total number of reads on forward strand
        revtotal = 0                            # total number of reads on reverse strand
        fwdtotals = {0: 0, 1: 0, 2: 0, 3: 0}    # reads with 0,1,2 and 3 mismatches on forward strand
        revtotals = {0: 0, 1: 0, 2: 0, 3: 0}    # reads with 0,1,2 and 3 mismatches on reverse strand

        MAX_MISMATCH = 3                        # max number of mismatches we will count
        TAG = "NM"                              # SAM tag which indicates number of mismatches (in this case edit dist)

        reads_file = ""
        try:
            for reads_file in reads_filenames:

                self.__logger.info("Counting mismatches in reads file: {}".format(reads_file))

                # read the bam file and iterate over each read
                samfile = pysam.AlignmentFile(reads_file, "rb")
                iteration = samfile.fetch(until_eof=True)

                for x in iteration:
                    # determine strand of read and get mismatch count
                    if (x.is_read1 and x.is_reverse) or (x.is_read2 and x.mate_is_reverse):
                        value = x.get_tag(TAG)
                        if value <= MAX_MISMATCH:
                            fwdtotals[value] += 1
                        fwdtotal += 1
                    elif (x.is_read1 and x.mate_is_reverse) or (x.is_read2 and x.is_reverse):
                        value = x.get_tag(TAG)
                        if value <= MAX_MISMATCH:
                            revtotals[value] += 1
                        revtotal += 1
                    else:   # unstranded - count everything as forward strand
                        value = x.get_tag(TAG)
                        if value <= MAX_MISMATCH:
                            fwdtotals[value] += 1
                        fwdtotal += 1

        except KeyError:
            self.__logger.info("*** Mismatch counts will not be reported *** No mismatch tag ({}) was found in the reads file {}".format(TAG, reads_file))

            # still need to count fwd and reverse reads though!
            for reads_file in reads_filenames:

                self.__logger.info("Counting reads only: {}".format(reads_file))

                # read the bam file and iterate over each read
                samfile = pysam.AlignmentFile(reads_file, "rb")
                iteration = samfile.fetch(until_eof=True)

                for x in iteration:
                    if (x.is_read1 and x.is_reverse) or (x.is_read2 and x.mate_is_reverse):
                        fwdtotal += 1
                    elif (x.is_read1 and x.mate_is_reverse) or (x.is_read2 and x.is_reverse):
                        revtotal += 1
                    else:  # unstranded - count everything as forward strand
                        fwdtotal += 1

        self.__logger.info("Mismatch counting completed")

        # put the results into a labelled structure
        ftotals = [fwdtotal]
        ftotals.extend(fwdtotals.values())
        rtotals = [revtotal]
        rtotals.extend(revtotals.values())
        results = {"forward strand": ftotals, "reverse strand": rtotals}

        return results

    def get_stats(self):
        """ Output the stats info.
        :return: an object containing stats for the reads files """

        return self.stats
