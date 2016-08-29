#!/usr/bin/env python

import pandas as pd
import math
import logging

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"


class Normaliser:
    def __init__(self, annotation, total_mapped_reads):
        """
        Initialise normalisation
        :param annotation: the annotation to take feature lengths and transcript counts from
        :param total_mapped_reads: total number of reads mapped in the experiment
        :return:
        """

        self.__logger = logging.getLogger(__name__)

        self._annotation = annotation
        self._total_mapped_reads = total_mapped_reads
        self._transcript_counts = annotation.get_transcripts_per_gene()

    def normalise_fpkm(self, s, ftr, fpkm_col_name):
        """
        Normalise series s corresponding to ftr, using FPKM.
        Assumes that s is ordered by feature, one line per feature, in the same order
        as the annotation
        :param s: series data to normalise
        :param ftr: ftr this series corresponds to in annotation
        :param fpkm_col_name: name to assign to normalised data
        :return: dataframe with normalised data
        """

        # for FPKM calc
        # C = Number of reads mapped to a ftr i.e. the variable
        # N = Total mapped reads in the experiment
        # L = ftr length in base-pairs for a gene - and what does that mean if there are multiple transcripts?
        # here we use length averaged over the transcripts
        # FPKM = (10^9 * C)/(N * L)

        self.__logger.debug("Normalising counts for {}".format(ftr))

        # constant name for temp column
        fpkm_div = "fpkm_div"

        # calculate average length of ftr across all its transcripts
        ftr_lengths = self._annotation.get_ftr_lengths(ftr, bygene=True)
        ave_lengths = ftr_lengths.div(self._transcript_counts.transcript_counts, axis='index')

        # calculate FPKM divisors as total experiment reads * length of each ftr
        fpkm_divisor = self._total_mapped_reads * ave_lengths
        fpkm_divisor.columns = [fpkm_div]

        # calculate normalised column
        fpkm_col = 10**9 * s.div(fpkm_divisor[fpkm_div], axis='index')
        fpkm_col = pd.DataFrame(fpkm_col)  # convert to dataframe
        fpkm_col.columns = [fpkm_col_name]

        # adjust normalisation in the face of high total reads
        total_order = math.log10(self._total_mapped_reads)
        if total_order > 6:
            fpkm_col = fpkm_col * (10**(total_order-6))

        return fpkm_col
