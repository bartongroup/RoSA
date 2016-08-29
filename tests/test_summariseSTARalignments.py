#!/usr/bin/python
"""
Unit tests for summariser module
"""

import unittest
import json
import os
from viewseq.summariseSTARAlignments import STARLogScraper

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"


class TestStarLogScraper(unittest.TestCase):

    """
    Test STARLogScraper class
    """

    def setUp(self):

        self.outfile = "test_STAR_output.json"

    def tearDown(self):
        """Tidy up files created for testing."""

        try:
            os.remove(self.outfile)
        except OSError as e:
            if e.errno == 2:
                # no such file, assume already deleted
                pass

    def test_scrapeSTARLog(self):

        """ Test typical behaviour of STAR log scraping with default names. """

        s = STARLogScraper()
        s.scrape("Data/STARlogtest", self.outfile)

        with open(self.outfile, 'r') as fp:

            logdata = json.load(fp)
            self.assertAlmostEqual(logdata["data"]["con1/rep1/Log.final.out"]["Uniquely-mapping (%)"], 89.53, 2)
            self.assertAlmostEqual(logdata["data"]["con1/rep2/Log.final.out"]["No. Mreads"], 78.729, 2)
            self.assertAlmostEqual(logdata["data"]["con1/rep3/Log.final.out"]["Unmapped - too short (%)"], 8.76, 2)
            self.assertAlmostEqual(logdata["data"]["test2/rep1/Log.final.out"]["Multi-mapping (%)"], 0.79, 2)
            self.assertItemsEqual(logdata["limits"].keys(),
                                  ["Multi-mapping (%)", "No. Mreads", "Uniquely-mapping (%)",
                                   "Unmapped - multi-mapping (%)", "Unmapped - other (%)", "Unmapped - too short (%)"])
            self.assertItemsEqual(logdata["conditions"],
                                  ["con1/rep1/Log.final.out", "con1/rep2/Log.final.out", "con1/rep3/Log.final.out",
                                   "test2/rep1/Log.final.out", "other_cons/reps/rep1/Log.final.out"])

            # check that limit was adjusted correctly
            self.assertAlmostEqual(logdata["limits"]["Unmapped - too short (%)"], 20, 1)

    def test_scrapeSTARLog_diff_log_name(self):

        """ Test STAR log scraping with non-standard log name. """

        s = STARLogScraper()
        s.scrape("Data/STARlogtest", self.outfile, logfilename="a_different_log_name")

        with open(self.outfile, 'r') as fp:
            logdata = json.load(fp)

            self.assertAlmostEqual(logdata["data"]["test2/rep1/a_different_log_name"]["Uniquely-mapping (%)"], 92.27, 2)
            self.assertItemsEqual(logdata["limits"].keys(),
                                  ["Multi-mapping (%)", "No. Mreads", "Uniquely-mapping (%)",
                                   "Unmapped - multi-mapping (%)", "Unmapped - other (%)", "Unmapped - too short (%)"])
            self.assertItemsEqual(logdata["conditions"], ["test2/rep1/a_different_log_name"])
