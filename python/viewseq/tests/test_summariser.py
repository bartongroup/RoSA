#!/usr/bin/python
"""
Unit tests for summariser module
"""

import unittest
import tempfile
import pandas as pd
import math
import shutil

import os.path
from viewseq.summariser import InputData, Summary
import viewseq
from viewseq.constants import *

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"


class TestSummariser(unittest.TestCase):

    """
    Test InputData and Summary classes, and functions in summariser module.
    """

    def setUp(self):

        """Initialize the framework for testing.

        Define and create a new system temporary directory with the
        `tempfile <http://docs.python.org/2/library/tempfile.html>`_ package
        for use as an existing directory.
        Creates empty gff, bam and wig files in this directory to test input data routine.
        """

        # setup a new temp directory
        self.existing_path = tempfile.mkdtemp()

        # setup temp output directory
        self.output_dir = "{}/summary_output".format(self.existing_path)

        # make an empty gff v3 file
        self.gff_filename = "{}/empty.gff".format(self.existing_path)
        temp_file = open(self.gff_filename, "w")
        temp_file.close()

        # make an empty bam file
        self.bam_filename = "{}/empty.bam".format(self.existing_path)
        temp_file = open(self.bam_filename, "w")
        temp_file.close()

        # make an empty wig file
        self.wig_filename = "{}/empty.bigwig".format(self.existing_path)
        temp_file = open(self.wig_filename, "w")
        temp_file.close()

        # make an empty unknown file
        self.unknown_filename = "{}/empty.unk".format(self.existing_path)
        temp_file = open(self.unknown_filename, "w")
        temp_file.close()

        # ensure a non-existent file
        self.non_filename = "{}/this_is_not_a_file.gff".format(self.existing_path)
        if os.path.exists(self.non_filename):  # being paranoid
            os.remove(self.non_filename)

        # make a very simple gff v3 file, borrowed from Nick's code
        self.simple_gff = "{}/simple_gff_v3.gff".format(self.existing_path)
        temp_file = open(self.simple_gff, "w")
        temp_file.write("##gff-version 3\n# These are \n# header\n# lines\n")
        temp_file.write("chrI\t.\tchromosome\t1\t10000\t.\t.\t.\t.\n")
        temp_file.write("chrI\t.\texon\t1000\t2000\t.\t-\t.\t.\n")
        temp_file.write("chrI\t.\tgene\t900\t2200\t.\t-\t.\t.\n")
        temp_file.write("chrI\t.\tmRNA\t1000\t2000\t.\t.-\t.\t.\n")
        temp_file.write("chr2\t.\tchromosome\t1\t20000\t.\t.\t.\t.\n")
        temp_file.write("chr2\t.\texon\t2000\t4000\t.\t+\t.\t.\n")
        temp_file.write("chr2\t.\tgene\t1800\t4400\t.\t+\t.\t.\n")
        temp_file.write("chr2\t.\texon\t4100\t4300\t.\t+\t.\t.\n")
        temp_file.write("chr3\t.\tchromosome\t1\t15000\t.\t.\t.\t.\n")
        temp_file.write("chr3\t.\texon\t5000\t5450\t.\t+\t.\t.\n")
        temp_file.write("chr3\t.\tgene\t5100\t5300\t.\t+\t.\t.\n")
        temp_file.close()

        # make an even simpler simple gff v3 file
        self.too_simple_gff = "{}/too_simple_gff_v3.gff".format(self.existing_path)
        temp_file = open(self.too_simple_gff, "w")
        temp_file.write("##gff-version 3\n# These are \n# header\n# lines\n")
        temp_file.write("chrI\t.\tchromosome\t1\t10000\t.\t.\t.\t.\n")
        temp_file.write("chrI\t.\texon\t1000\t2000\t.\t-\t.\t.\n")
        temp_file.write("chrI\t.\tgene\t900\t2200\t.\t-\t.\t.\n")
        temp_file.write("chr2\t.\tchromosome\t1\t20000\t.\t.\t.\t.\n")
        temp_file.write("chr2\t.\texon\t2000\t4000\t.\t+\t.\t.\n")
        temp_file.write("chr2\t.\tgene\t1800\t4400\t.\t+\t.\t.\n")
        temp_file.write("chr2\t.\texon\t4100\t4300\t.\t+\t.\t.\n")
        temp_file.write("chr3\t.\tchromosome\t1\t15000\t.\t.\t.\t.\n")
        temp_file.write("chr3\t.\texon\t5000\t5450\t.\t+\t.\t.\n")
        temp_file.write("chr3\t.\tgene\t5100\t5300\t.\t+\t.\t.\n")
        temp_file.close()

        # make a less simple gff v3 file
        self.simple2_gff = "{}/simple2_gff_v3.gff".format(self.existing_path)
        temp_file = open(self.simple2_gff, "w")
        temp_file.write("##gff-version 3\n# These are \n# header\n# lines\n")
        temp_file.write("chrI\t.\tchromosome\t1\t10000\t.\t.\t.\t.\n")
        temp_file.write("chrI\t.\texon\t1000\t2000\t.\t-\t.\t.\n")
        temp_file.write("chrI\t.\tgene\t900\t2200\t.\t-\t.\t.\n")
        temp_file.write("chr2\t.\tchromosome\t1\t20000\t.\t.\t.\t.\n")
        temp_file.write("chr2\t.\texon\t2000\t4000\t.\t+\t.\t.\n")
        temp_file.write("chr2\t.\tgene\t1800\t4400\t.\t+\t.\t.\n")
        temp_file.write("chr2\t.\texon\t4100\t4300\t.\t+\t.\t.\n")
        temp_file.write("chr3\t.\tchromosome\t1\t15000\t.\t.\t.\t.\n")
        temp_file.write("chr3\t.\texon\t5000\t5450\t.\t+\t.\t.\n")
        temp_file.write("chr3\t.\tgene\t5100\t5300\t.\t+\t.\t.\n")
        temp_file.close()

        # make the sequence ontology example gff v3 file
        self.seqont_gff = "{}/seqont_gff_v3.gff".format(self.existing_path)
        temp_file = open(self.seqont_gff, "w")
        temp_file.write("##gff-version 3\n##sequence-region		ctg123	1	1497228\n")
        temp_file.write("ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001;Name=EDEN\n")
        temp_file.write("ctg123\t.\tTF_binding_site\t1000\t1012\t.\t+\t.\tID=tfbs00001;Parent=gene00001\n")
        temp_file.write("ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00001;Parent=gene00001;Name=EDEN.1\n")
        temp_file.write("ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00002;Parent=gene00001;Name=EDEN.2\n")
        temp_file.write("ctg123\t.\tmRNA\t1300\t9000\t.\t+\t.\tID=mRNA00003;Parent=gene00001;Name=EDEN.3\n")
        temp_file.write("ctg123\t.\texon\t1300\t1500\t.\t+\t.\tID=exon00001;Parent=mRNA00003\n")
        temp_file.write("ctg123\t.\texon\t1050\t1500\t.\t+\t.\tID=exon00002;Parent=mRNA00001,mRNA00002\n")
        temp_file.write("ctg123\t.\texon\t3000\t3902\t.\t+\t.\tID=exon00003;Parent=mRNA00001,mRNA00003\n")
        temp_file.write("ctg123\t.\texon\t5000\t5500\t.\t+\t.\tID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003\n")
        temp_file.write("ctg123\t.\texon\t7000\t9000\t.\t+\t.\tID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003\n")
        temp_file.write("ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00001;Parent=mRNA00001;Name=edenprotein.1\n")
        temp_file.write("ctg123\t.\tCDS\t3000\t3902\t.\t+\t0\tID=cds00001;Parent=mRNA00001;Name=edenprotein.1\n")
        temp_file.write("ctg123\t.\tCDS\t5000\t5500\t.\t+\t0\tID=cds00001;Parent=mRNA00001;Name=edenprotein.1\n")
        temp_file.write("ctg123\t.\tCDS\t7000\t7600\t.\t+\t0\tID=cds00001;Parent=mRNA00001;Name=edenprotein.1\n")
        temp_file.write("ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00002;Parent=mRNA00002;Name=edenprotein.2\n")
        temp_file.write("ctg123\t.\tCDS\t5000\t5500\t.\t+\t0\tID=cds00002;Parent=mRNA00002;Name=edenprotein.2\n")
        temp_file.write("ctg123\t.\tCDS\t7000\t7600\t.\t+\t0\tID=cds00002;Parent=mRNA00002;Name=edenprotein.2\n")
        temp_file.write("ctg123\t.\tCDS\t3301\t3902\t.\t+\t0\tID=cds00003;Parent=mRNA00003;Name=edenprotein.3\n")
        temp_file.write("ctg123\t.\tCDS\t5000\t5500\t.\t+\t1\tID=cds00003;Parent=mRNA00003;Name=edenprotein.3\n")
        temp_file.write("ctg123\t.\tCDS\t7000\t7600\t.\t+\t1\tID=cds00003;Parent=mRNA00003;Name=edenprotein.3\n")
        temp_file.write("ctg123\t.\tCDS\t3391\t3902\t.\t+\t0\tID=cds00004;Parent=mRNA00003;Name=edenprotein.4\n")
        temp_file.write("ctg123\t.\tCDS\t5000\t5500\t.\t+\t1\tID=cds00004;Parent=mRNA00003;Name=edenprotein.4\n")
        temp_file.write("ctg123\t.\tCDS\t7000\t7600\t.\t+\t1\tID=cds00004;Parent=mRNA00003;Name=edenprotein.4\n")
        temp_file.close()

        # a similar gff file with 3' and 5' UTRs
        self.utrseqont_gff = "{}/utrseqont_gff_v3.gff".format(self.existing_path)
        temp_file = open(self.utrseqont_gff, "w")
        temp_file.write("##gff-version 3\n##sequence-region		ctg123	1	1497228\n")
        temp_file.write("ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001;Name=EDEN\n")
        temp_file.write("ctg123\t.\tTF_binding_site\t1000\t1012\t.\t+\t.\tID=tfbs00001;Parent=gene00001\n")
        temp_file.write("ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00001;Parent=gene00001;Name=EDEN.1\n")
        temp_file.write("ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00002;Parent=gene00001;Name=EDEN.2\n")
        temp_file.write("ctg123\t.\tmRNA\t1300\t9000\t.\t+\t.\tID=mRNA00003;Parent=gene00001;Name=EDEN.3\n")
        temp_file.write("ctg123\t.\texon\t1300\t1500\t.\t+\t.\tID=exon00001;Parent=mRNA00003\n")
        temp_file.write("ctg123\t.\texon\t1050\t1500\t.\t+\t.\tID=exon00002;Parent=mRNA00001,mRNA00002\n")
        temp_file.write("ctg123\t.\texon\t3000\t3902\t.\t+\t.\tID=exon00003;Parent=mRNA00001,mRNA00003\n")
        temp_file.write("ctg123\t.\texon\t5000\t5500\t.\t+\t.\tID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003\n")
        temp_file.write("ctg123\t.\texon\t7000\t9000\t.\t+\t.\tID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003\n")
        temp_file.write("ctg123\t.\tfive_prime_UTR\t1050\t1200\t.\t+\t0\tParent=mRNA00001,mRNA00002\n")
        temp_file.write("ctg123\t.\tfive_prime_UTR\t1300\t1499\t.\t+\t0\tParent=mRNA00003\n")
        temp_file.write("ctg123\t.\tfive_prime_UTR\t3000\t3300\t.\t+\t0\tParent=mRNA00003\n")
        temp_file.write("ctg123\t.\tfive_prime_UTR\t3000\t3390\t.\t+\t0\tParent=mRNA00003\n")
        temp_file.write("ctg123\t.\tthree_prime_UTR\t7601\t9000\t.\t+\t0\tParent=mRNA00001,mRNA00002,mRNA00003\n")
        temp_file.close()

        self.config_yaml = "{}/config.yaml".format(self.existing_path)
        temp_file = open(self.config_yaml, "w")
        temp_file.write("annotation: Data/mini.gff\n")
        temp_file.write("readsfiles: \n")
        temp_file.write(" - Data/mini.bam\n")
        temp_file.write("outputdir: {}\n".format(self.output_dir))
        temp_file.write("featureCounts_exe: featureCounts\n")
        temp_file.write("featureCounts_options: -s 2\n")
        temp_file.write("STAR_rootpath: Data/STARlogtest\n")
        temp_file.write("STAR_logfilename: Log.final.out\n")

        temp_file.close()

    def tearDown(self):
        """Tidy up files created for testing."""

        os.remove(self.gff_filename)
        os.remove(self.bam_filename)
        os.remove(self.wig_filename)
        os.remove(self.unknown_filename)
        os.remove(self.config_yaml)

        try:
            shutil.rmtree("{}".format(self.output_dir))
        except OSError:
            pass      # assume directory not found

    def test_add_gff(self):

        """ Test the adding of a new gff file to InputData."""
        d = InputData()

        # successfully locate and identify gff file
        d.add_input_file(self.gff_filename)

    def test_add_bam(self):

        """ Test the adding of a new bam file to InputData."""
        d = InputData()

        # successfully locate and identify bam file
        d.add_input_file(self.bam_filename)

    def test_add_wig(self):

        """ Test the adding of a new wig file to InputData."""
        d = InputData()

        # successfully locate and identify wig file
        d.add_input_file(self.wig_filename)

    def test_add_all(self):

        """ Test the adding of a new gff file to InputData."""
        d = InputData()

        # successfully locate and identify gff file
        d.add_input_file(self.gff_filename)

        # successfully locate and identify bam file
        d.add_input_file(self.bam_filename)

        # successfully locate and identify wig file
        d.add_input_file(self.wig_filename)

    def test_add_unknown(self):

        """ Test the adding of a file with unknown extension to InputData."""
        d = InputData()

        # reject file
        self.assertRaises(InputData.UnknownFileExtError, d.add_input_file, self.unknown_filename)

    def test_add_input_no_file(self):

        """ Test trying to add a non-existent input file to InputData."""
        d = InputData()

        # fail to locate non-existent file
        self.assertRaises(IOError, d.add_input_file, self.non_filename)

    def test_get_annotation_none(self):

        """ Test no annotation is returned if it wasn't initialised yet. """
        d = InputData()

        self.assertIsNone(d.get_annotation())

    def test_get_annotation(self):

        """ Test annotation object correctly returned. """
        d = InputData()
        d.add_input_file(self.simple_gff)
        p = d.get_annotation()

        self.assertEqual(type(p).__name__, "GffParser")

    def test_get_no_feature_counts(self):

        """ Test no counts returned if no reads files were supplied. """
        d = InputData()
        self.assertIsNone(d.get_feature_counts(None, "", "", ""))

    def test_get_feature_counts(self):

        """ Test feature counts object correctly returned. """
        d = InputData()
        d.add_input_file("Data/mini.gff")
        d.add_input_file("Data/mini.bam")
        a = d.get_annotation()

        f = d.get_feature_counts(a, "featureCounts", "-p -s 2")
        self.assertIsInstance(f, viewseq.fcparser.FeatureCountParser)

    def test_get_feature_counts_mismatch(self):

        """ Test IOError caught if chromosome names don't match. """
        d = InputData()
        d.add_input_file(self.simple_gff)
        d.add_input_file("Data/mini.bam")
        a = d.get_annotation()

        self.assertIsNone(d.get_feature_counts(a, "featureCounts", "-p -s 2"))

    def test_get_noreads_data(self):

        """ Test reads data is None when there are no reads files. """
        d = InputData()
        self.assertIsNone(d.get_reads_data())

    def test_get_reads_data(self):

        """ Test reads data is returned correctly. """
        d = InputData()
        d.add_input_file("Data/mini.bam")
        r = d.get_reads_data()
        self.assertIsInstance(r, viewseq.samstats.SamStats)

    def test_vis_output(self):

        """ Test visualisation output. """

        d = InputData()
        d.add_input_file("Data/mini.gff")
        d.add_input_file("Data/mini.bam")
        s = Summary(d, "featureCounts", "-s 2")

        s.output_for_visualisation(self.output_dir)

        self._check_vis_output()

    def _check_vis_output(self):

        # Check each file is output and has expected contents
        visfile = os.path.join(self.output_dir, "genes.json")
        self.assertTrue(os.path.exists(visfile), "JSON file for genes data does not exist")
        testdata = pd.read_json(visfile, orient="split")
        self.assertListEqual(testdata["exon"].tolist(),
                             [181, 0, 1, 102, 0])
        self.assertListEqual(testdata["intron"].tolist()[:4],
                             [3, 0, 0, 0])
        self.assertTrue(math.isnan(testdata["intron"].tolist()[4]))
        self.assertListEqual(testdata[EXON_COUNTS].tolist(),
                             [6, 12, 2, 23, 1])
        self.assertListEqual(testdata[TRANSCRIPT_COUNTS].tolist(),
                             [1, 2, 1, 2, 1])

        visfile = os.path.join(self.output_dir, "exon-data.json")
        self.assertTrue(os.path.exists(visfile), "JSON file for exon lengths does not exist")
        testdata = pd.read_json(visfile, orient="split")
        self.assertItemsEqual(testdata[LENGTH].tolist(),
                             [283, 281, 120, 390, 153, 461, 167, 48, 90, 46, 74, 86, 67, 76, 633, 336, 294, 280, 380,
                              1525, 1306, 114, 211, 395, 220, 173, 123, 161, 234, 151, 183, 162, 96, 629, 98, 191, 906,
                              165, 407, 326, 1036, 165, 219, 207])

        visfile = os.path.join(self.output_dir, "intron-data.json")
        self.assertTrue(os.path.exists(visfile), "JSON file for intron lengths does not exist")
        testdata = pd.read_json(visfile, orient="split")
        self.assertItemsEqual(testdata[LENGTH].tolist(),
                             [82, 209, 100, 78, 112, 90, 96, 78, 88, 81, 83, 88, 90, 85, 86, 90, 84, 89, 276, 84, 79,
                              81, 98, 85, 81, 173, 87, 151, 113, 112, 106, 248, 91, 106, 161])

        visfile = os.path.join(self.output_dir, "firstexons.json")
        self.assertTrue(os.path.exists(visfile), "JSON file for first exons does not exist")
        testdata = pd.read_json(visfile, orient="split")
        self.assertEqual(len(testdata), 7)

        visfile = os.path.join(self.output_dir, "lastexons.json")
        self.assertTrue(os.path.exists(visfile), "JSON file for last exons does not exist")
        testdata = pd.read_json(visfile, orient="split")
        self.assertEqual(len(testdata), 7)

    def test_vis_output_with_gtf(self):

        """ Test visualisation output. """

        d = InputData()
        d.add_input_file("Data/cds3.gtf")
        d.add_input_file("Data/mini.bam")
        s = Summary(d, "featureCounts", "-s 2 -S fr")

        s.output_for_visualisation(self.output_dir)

        # Check each file is output and has expected contents
        visfile = os.path.join(self.output_dir, "genes.json")
        self.assertTrue(os.path.exists(visfile), "JSON file for genes data does not exist")
        testdata = pd.read_json(visfile, orient="split")

        # Here just check we actually got 3 genes output
        self.assertEqual(len(testdata), 3)

        self.assertListEqual(testdata[EXON_COUNTS].tolist(), [6, 4, 1])
        self.assertListEqual(testdata[LENGTH].tolist(), [804, 1266, 9903])

    def test_vis_output_with_gtf_and_hash(self):

        """ Test visualisation output. """

        d = InputData()
        d.add_input_file("Data/cds4.gtf")
        d.add_input_file("Data/mini.bam")
        s = Summary(d, "featureCounts", "-s 2 -S fr")

        s.output_for_visualisation(self.output_dir)

        # Check each file is output and has expected contents
        visfile = os.path.join(self.output_dir, "genes.json")
        self.assertTrue(os.path.exists(visfile), "JSON file for genes data does not exist")
        testdata = pd.read_json(visfile, orient="split")

        # Here just check we actually got 4 genes output
        self.assertEqual(len(testdata), 4)

        self.assertListEqual(testdata[EXON_COUNTS].tolist(), [1, 6, 4, 1])
        self.assertListEqual(testdata[LENGTH].tolist(), [601, 804, 1266, 9903])

    def test_build_summary(self):

        viewseq.summariser.build_summary("Data/mini.gff", ["Data/mini.bam"], "", output=self.output_dir,
                                            fc_options="-s 2")
        self._check_vis_output()

    def test_main(self):

        viewseq.summariser.run_main([self.config_yaml])

        self._check_vis_output()


if __name__ == '__main__':
    unittest.main(verbosity=2)
