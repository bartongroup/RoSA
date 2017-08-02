#!/usr/bin/python
"""
Unit tests for antisense module
"""

import unittest
import csv
import pandas as pd
from antisense import Antisense

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"


class TestAntisense(unittest.TestCase):

    """
    Test functions in antisense module.
    """

    def setUp(self):

        self.exondata = "Data/exon-anti-counts-km.out"
        self.introndata = "Data/intron-anti-counts-km.out"

        self.overlaps = "Data/overlapping_genes_Arabidopsis-3.txt"
        self.gtf = "Data/genes_fixed.gtf"

    def test_id_antisense(self):

        a = Antisense()
        result = a.id_structured_antisense(self.exondata, self.introndata, "result.csv")

    def test_get_diff_strand_overlaps(self):

        a = Antisense()
        diff_overlaps = a.get_diff_strand_overlaps(self.overlaps, self.gtf)
        diff_overlaps.to_csv(path="Data/diff_overlaps.csv", sep='\t', index=False)

    def test_filter(self):

        # needs import csv

        a = Antisense()
        filtered = a.filter_gtf("Data/diff_overlaps.csv", "Data/genes_fixed.gtf")
        filtered.to_csv(path_or_buf="Data/filtered.gtf", sep="\t", index=False, quoting=csv.QUOTE_NONE)

    def test_find_reads(self):

        a = Antisense()
        duplicates = a.find_anti_reads_in_sense("/Users/kmourao/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Col.fwd.assigned.filtered.anti.reads",
                                                "/Users/kmourao/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Col.fwd.assigned.filtered.sense.reads")
        duplicates.to_csv(path_or_buf="/Users/kmourao/Documents/Arabidopsis_RNAMeth_polyA/antisense_counting/Col_duplicates.csv", sep="\t", index=False)

    def test_find_anti_with_sense_structure(self):

        a = Antisense()
        a.find_anti_with_sense_structure("Data/mini_genes_fixed.gtf",
                                         "Data/mini-anti.bam",
                                         "Data/anti-results.csv")

        results = pd.read_csv("Data/anti-results.csv")
        self.assertListEqual(results["anticounts"].tolist(), [1, 3, 2])
        self.assertListEqual(results["sensecounts"].tolist(), [871, 123, 0])

    def test_find_anti_with_sense_structure_overlapping(self):
        a = Antisense()
        a.find_anti_with_sense_structure("Data/genes_fixed.gtf",
                                         "Data/AT1G68568.bam",
                                         "Data/anti-results.csv")

        results = pd.read_csv("Data/anti-results.csv")
        self.assertListEqual(results["anticounts"].tolist(),
                             [0.0, 0.0, 0.0, 4.0, 12.0, 10.0, 32.0, 18.0, 28.0, 44.0, 39.0, 283.0, 75.0, 30.0,
                              128.0, 8.0, 20.0, 87.0, 965.0, 673.0, 6.0, 366.0])
        self.assertListEqual(results["sensecounts"].tolist(),
                             [34, 101, 6, 2729, 5384, 3517, 9426, 4735, 7155, 10958, 9069, 65159, 17238, 6646,
                              21193, 1306, 3068, 12570, 78262, 51969, 445, 25253])
