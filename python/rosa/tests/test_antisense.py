#!/usr/bin/python

#  Copyright 2017 Kira Mourao, Nick Schurch
#
#  This file is part of RoSA.
#
#  RoSA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  RoSA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with RoSA.  If not, see <http://www.gnu.org/licenses/>.

"""
Unit tests for antisense module
"""

import unittest
import csv
import pandas as pd
from rosa.antisense import Antisense

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
        self.assertListEqual(results["anticounts"].tolist(), [1])
        self.assertListEqual(results["sensecounts"].tolist(), [871])

    def test_find_anti_with_sense_structure_overlapping(self):
         a = Antisense()
         a.find_anti_with_sense_structure("Data/genes_fixed.gtf",
                                          "Data/AT1G68568.bam",
                                          "Data/anti-results.csv")

         results = pd.read_csv("Data/anti-results.csv")

         # beware, these results may change with different versions of sambamba etc.
         # These results with sambamba 0.6.6
         self.assertListEqual(results["anticounts"].tolist(),
                              [0.0, 0.0, 0.0, 0.0, 4.0, 12.0, 10.0, 32.0, 18.0, 28.0, 44.0, 39.0, 283.0, 30.0,
                               128.0, 8.0, 136.0, 20.0, 87.0, 964.0, 673.0, 6.0, 364.0])
         self.assertListEqual(results["sensecounts"].tolist(),
                              [34, 101, 6, 10, 2729, 5381, 3512, 9419, 4731, 7136, 10947, 9060, 65106, 6644,
                               21173, 1306, 21373, 3068, 12554, 78167, 519631, 445, 25226])
