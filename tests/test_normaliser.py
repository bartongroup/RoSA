#!/usr/bin/python
"""
Unit tests for normaliser
"""

import unittest
import pandas as pd

from viewseq.normaliser import Normaliser
from viewseq.gffparser import GffParser

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"


class TestNormaliser(unittest.TestCase):

    """
    Test normalisation functions.
    """

    def setUp(self):
        self.test_data_gff = "Data/mini.gff"

    def test_normalise_fpkm(self):
        """ Test FPKM normalisation. """

        # set up annotation
        p = GffParser(self.test_data_gff)

        # set up data series
        data = {"AT1G01010": 217, "AT1G01020": 0, "AT1G01030": 0, "AT1G01040": 173, "AT1G01046": 1}
        s = pd.Series(data)
        s.name = p.featureType.exon

        # initialise normaliser
        nm = Normaliser(p, 1000)

        normed = nm.normalise_fpkm(s, p.featureType.exon, "exon-FPKM")

        self.assertListEqual(normed["exon-FPKM"].tolist(), [128554.50236966823, 0.0, 0.0,
                                                            45104.940685699388, 4830.9178743961356])

    def test_normalise_high_reads(self):
        """ Test FPKM normalisation when total reads is very high. """

        # set up annotation
        p = GffParser(self.test_data_gff)

        # set up data series
        data = {"AT1G01010": 42, "AT1G01020": 66, "AT1G01030": 69, "AT1G01040": 461, "AT1G01046": 1}
        s = pd.Series(data)
        s.name = p.featureType.exon

        # initialise normaliser
        nm = Normaliser(p, 1694457307)

        normed = nm.normalise_fpkm(s, p.featureType.exon, "exon-FPKM")

        self.assertListEqual(normed["exon-FPKM"].tolist(), [24.881516587677716, 60.081929904415084,
                                                            36.22047244094486, 120.19293442836653,
                                                            4.8309178743961327])

