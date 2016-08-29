#!/usr/bin/python
"""
Unit tests for samstats module
"""

import unittest
import locale
from viewseq.samstats import SamStats

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"


class TestSamStats(unittest.TestCase):

    """
    Test functions in samstats module.
    """

    def setUp(self):

        """Initialize the framework for testing.

        Define location of bam file for testing
        """

        self.test_data_bam = "Data/mini.bam"
        self.test_data_nM_bam = "Data/mini-nM.bam"
        self.test_data_unstranded = "Data/unstranded.bam"
        self.test_data_nMunstranded = "Data/mini-nM-nostrand.bam"
        locale.setlocale(locale.LC_ALL, 'en_GB.UTF-8')

    def tearDown(self):
        """Tidy up files created for testing."""

    def test_samstats(self):

        read_files = [self.test_data_bam]
        s = SamStats(read_files)
        self.assertIsNotNone(s)

    def test_readcounts(self):

        read_files = [self.test_data_bam]
        s = SamStats(read_files)
        stats = s.get_stats()

        self.assertItemsEqual(stats["forward strand"],
            [("Total reads", locale.format("%d", 5515, grouping=True)),
            ("0 mismatches", locale.format("%d", 4854, grouping=True)),
            ("1 mismatch", locale.format("%d", 407, grouping=True)),
            ("2 mismatches", locale.format("%d", 130, grouping=True)),
            ("3 mismatches", locale.format("%d", 40, grouping=True))])
        self.assertItemsEqual(stats["reverse strand"],
            [("Total reads", locale.format("%d", 2, grouping=True)),
            ("0 mismatches", locale.format("%d", 2, grouping=True)),
            ("1 mismatch", locale.format("%d", 0, grouping=True)),
            ("2 mismatches", locale.format("%d", 0, grouping=True)),
            ("3 mismatches", locale.format("%d", 0, grouping=True))])

    def test_missing_readcounts(self):
        read_files = [self.test_data_nM_bam]
        s = SamStats(read_files)
        stats = s.get_stats()

        self.assertItemsEqual(stats["forward strand"],
            [("Total reads", locale.format("%d", 5516, grouping=True)),
            ("0 mismatches", locale.format("%d", 0, grouping=True)),
            ("1 mismatch", locale.format("%d", 0, grouping=True)),
            ("2 mismatches", locale.format("%d", 0, grouping=True)),
            ("3 mismatches", locale.format("%d", 0, grouping=True))])
        self.assertItemsEqual(stats["reverse strand"],
            [("Total reads", locale.format("%d", 1, grouping=True)),
            ("0 mismatches", locale.format("%d", 0, grouping=True)),
            ("1 mismatch", locale.format("%d", 0, grouping=True)),
            ("2 mismatches", locale.format("%d", 0, grouping=True)),
            ("3 mismatches", locale.format("%d", 0, grouping=True))])

    def test_unstranded(self):
        """ Test counts are made for unstranded data, and put into forward strand counts. """

        read_files = [self.test_data_unstranded]
        s = SamStats(read_files)
        stats = s.get_stats()

        self.assertItemsEqual(stats["forward strand"],
                              [("Total reads", locale.format("%d", 2159, grouping=True)),
                               ("0 mismatches", locale.format("%d", 1660, grouping=True)),
                               ("1 mismatch", locale.format("%d", 209, grouping=True)),
                               ("2 mismatches", locale.format("%d", 290, grouping=True)),
                               ("3 mismatches", locale.format("%d", 0, grouping=True))])
        self.assertItemsEqual(stats["reverse strand"],
                              [("Total reads", locale.format("%d", 0, grouping=True)),
                               ("0 mismatches", locale.format("%d", 0, grouping=True)),
                               ("1 mismatch", locale.format("%d", 0, grouping=True)),
                               ("2 mismatches", locale.format("%d", 0, grouping=True)),
                               ("3 mismatches", locale.format("%d", 0, grouping=True))])

        read_files = [self.test_data_nMunstranded]
        s = SamStats(read_files)
        stats = s.get_stats()

        self.assertItemsEqual(stats["forward strand"],
                              [("Total reads", locale.format("%d", 5517, grouping=True)),
                               ("0 mismatches", locale.format("%d", 0, grouping=True)),
                               ("1 mismatch", locale.format("%d", 0, grouping=True)),
                               ("2 mismatches", locale.format("%d", 0, grouping=True)),
                               ("3 mismatches", locale.format("%d", 0, grouping=True))])
        self.assertItemsEqual(stats["reverse strand"],
                              [("Total reads", locale.format("%d", 0, grouping=True)),
                               ("0 mismatches", locale.format("%d", 0, grouping=True)),
                               ("1 mismatch", locale.format("%d", 0, grouping=True)),
                               ("2 mismatches", locale.format("%d", 0, grouping=True)),
                               ("3 mismatches", locale.format("%d", 0, grouping=True))])


