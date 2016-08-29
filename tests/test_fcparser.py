#!/usr/bin/python
"""
Unit tests for fcparser module
"""

import unittest
import tempfile
import os.path
import math
from viewseq.fcparser import FeatureCountParser
from viewseq.gffparser import GffParser
from viewseq.gtfparser import GtfParser

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"


class TestFCParser(unittest.TestCase):

    """
    Test functions in parser module.
    """

    def setUp(self):

        """Initialize the framework for testing.

        Define and create a new system temporary directory with the
        `tempfile <http://docs.python.org/2/library/tempfile.html>`_ package
        for use as an existing directory.
        Creates featurecounts files in this directory to test parsing functionality.
        """

        self.test_data_gff = "Data/mini.gff"
        self.test_data_bam = "Data/mini.bam"
        self.test_unstranded_bam = "Data/unstranded.bam"
        self.test_unstranded_gtf = "Data/mm9.gtf"
        self.test_nofive_gff = "Data/nofiveUTR.gff"

        # setup a new temp directory
        self.existing_path = tempfile.mkdtemp()

        # make a featurecounts file
        self.fcout = "{}/fcout.out".format(self.existing_path)
        temp_file = open(self.fcout, "w")
        temp_file.write("# Program:featureCounts v1.4.0; Command:\"/sw/opt/subread-1.4.0/bin/featureCounts\" \"-p\" "
                        "\"-P\" \"-b\" \"-T\" \"8\" \"-s\" \"2\" \"-S\" \"-a\" \"/something.gtf\" \"-o\" "
                        "\"/some_file_path/featureCounts.out\" \"/some_file_path.bam\"\n")
        temp_file.write("Geneid\tChr\tStart\tEnd\tStrand\tLength\tA_bam_file_path.bam\n")
        temp_file.write("AT1G01010\t1\t3631\t3913\t+\t100\t31\n")
        temp_file.write("AT1G01010\t1\t3996\t4276\t+\t65\t4\n")
        temp_file.write("AT1G01010\t1\t4486\t4605\t+\t103\t37\n")
        temp_file.write("AT1G01010\t1\t4706\t5095\t+\t154\t66\n")
        temp_file.write("AT1G01010\t1\t5174\t5326\t+\t450\t35\n")
        temp_file.write("AT1G01010\t1\t5439\t5899\t+\t211\t55\n")

        temp_file.write("AT1G01020\t1\t5928\t6263\t-\t98\t22\n")
        temp_file.write("AT1G01020\t1\t6437\t7069\t-\t154\t24\n")
        temp_file.write("AT1G01020\t1\t7157\t7450\t-\t79\t35\n")
        temp_file.write("AT1G01020\t1\t7564\t7649\t-\t311\t111\n")
        temp_file.write("AT1G01020\t1\t7762\t7835\t-\t498\t136\n")
        temp_file.write("AT1G01020\t1\t7942\t7987\t-\t76\t5\n")
        temp_file.write("AT1G01020\t1\t8236\t8325\t-\t93\t43\n")
        temp_file.write("AT1G01020\t1\t8417\t8464\t-\t14\t1\n")
        temp_file.write("AT1G01020\t1\t8571\t8737\t-\t65\t21\n")

        temp_file.close()

        # a gff with overlapping +ve and -ve strand exons, so we can create annotations
        self.overlapping_gff = "{}/overlapping_gff_v3.gff".format(self.existing_path)
        temp_file = open(self.overlapping_gff, "w")
        temp_file.write("##gff-version 3\n##sequence-region		ctg123	1	1497228\n")
        temp_file.write("Chr1\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001;Name=EDEN\n")
        temp_file.write("Chr1\t.\tTF_binding_site\t1000\t1012\t.\t+\t.\tID=tfbs00001;Parent=gene00001\n")
        temp_file.write("Chr1\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00001;Parent=gene00001;Name=EDEN.1\n")
        temp_file.write("Chr1\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00002;Parent=gene00001;Name=EDEN.2\n")
        temp_file.write("Chr1\t.\tmRNA\t1300\t9000\t.\t+\t.\tID=mRNA00003;Parent=gene00001;Name=EDEN.3\n")
        temp_file.write("Chr1\t.\texon\t1300\t1500\t.\t+\t.\tID=exon00001;Parent=mRNA00003\n")
        temp_file.write("Chr1\t.\texon\t1050\t1500\t.\t+\t.\tID=exon00002;Parent=mRNA00001,mRNA00002\n")
        temp_file.write("Chr1\t.\texon\t3000\t3902\t.\t+\t.\tID=exon00003;Parent=mRNA00001,mRNA00003\n")
        temp_file.write("Chr1\t.\texon\t5000\t5500\t.\t+\t.\tID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003\n")
        temp_file.write("Chr1\t.\texon\t7000\t9000\t.\t+\t.\tID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003\n")
        temp_file.write("Chr1\t.\tfive_prime_UTR\t1050\t1200\t.\t+\t0\tParent=mRNA00001,mRNA00002\n")
        temp_file.write("Chr1\t.\tfive_prime_UTR\t1300\t1499\t.\t+\t0\tParent=mRNA00003\n")
        temp_file.write("Chr1\t.\tfive_prime_UTR\t3000\t3300\t.\t+\t0\tParent=mRNA00003\n")
        temp_file.write("Chr1\t.\tfive_prime_UTR\t3000\t3390\t.\t+\t0\tParent=mRNA00003\n")
        temp_file.write("Chr1\t.\tthree_prime_UTR\t7601\t9000\t.\t+\t0\tParent=mRNA00001,mRNA00002,mRNA00003\n")

        temp_file.write("Chr1\t.\tgene\t5550\t12000\t.\t-\t.\tID=gene00002;Name=OPP\n")
        temp_file.write("Chr1\t.\tmRNA\t5590\t12000\t.\t-\t.\tID=mRNA00006;Parent=gene00002\n")
        temp_file.write("Chr1\t.\tmRNA\t9000\t12000\t.\t-\t.\tID=mRNA00007;Parent=gene00002\n")
        temp_file.write("Chr1\t.\texon\t5590\t7000\t.\t-\t.\tID=exon00006;Parent=mRNA00006\n")
        temp_file.write("Chr1\t.\texon\t9000\t12000\t.\t-\t.\tID=exon00007;Parent=mRNA00006,mRNA00007\n")
        temp_file.write("Chr1\t.\tthree_prime_UTR\t5590\t6000\t.\t-\t.\tParent=mRNA00006\n")
        temp_file.write("Chr1\t.\tthree_prime_UTR\t9000\t9200\t.\t-\t.\tParent=mRNA00007\n")
        temp_file.write("Chr1\t.\tfive_prime_UTR\t11801\t12000\t.\t-\t.\tParent=mRNA00006,mRNA00007\n")

        temp_file.write("Chr2\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00003;Name=EDEN_AGAIN\n")
        temp_file.write("Chr2\t.\tTF_binding_site\t1000\t1012\t.\t+\t.\tID=tfbs00002;Parent=gene00003\n")
        temp_file.write("Chr2\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00008;Parent=gene00003;Name=EDEN_AGAIN.1\n")
        temp_file.write("Chr2\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00009;Parent=gene00003;Name=EDEN_AGAIN.2\n")
        temp_file.write("Chr2\t.\tmRNA\t1300\t9000\t.\t+\t.\tID=mRNA00010;Parent=gene00003;Name=EDEN_AGAIN.3\n")
        temp_file.write("Chr2\t.\texon\t1300\t1500\t.\t+\t.\tID=exon00001;Parent=mRNA00010\n")
        temp_file.write("Chr2\t.\texon\t1050\t1500\t.\t+\t.\tID=exon00002;Parent=mRNA00008,mRNA00009\n")
        temp_file.write("Chr2\t.\texon\t3000\t3902\t.\t+\t.\tID=exon00003;Parent=mRNA00008,mRNA00010\n")
        temp_file.write("Chr2\t.\texon\t5000\t5500\t.\t+\t.\tID=exon00004;Parent=mRNA00008,mRNA00009,mRNA00010\n")
        temp_file.write("Chr2\t.\texon\t7000\t9000\t.\t+\t.\tID=exon00005;Parent=mRNA00008,mRNA00009,mRNA00010\n")
        temp_file.write("Chr2\t.\tfive_prime_UTR\t1050\t1200\t.\t+\t0\tParent=mRNA00008,mRNA00009\n")
        temp_file.write("Chr2\t.\tfive_prime_UTR\t1300\t1499\t.\t+\t0\tParent=mRNA00010\n")
        temp_file.write("Chr2\t.\tfive_prime_UTR\t3000\t3300\t.\t+\t0\tParent=mRNA00010\n")
        temp_file.write("Chr2\t.\tfive_prime_UTR\t3000\t3390\t.\t+\t0\tParent=mRNA00010\n")
        temp_file.write("Chr2\t.\tthree_prime_UTR\t7601\t9000\t.\t+\t0\tParent=mRNA00008,mRNA00009,mRNA00010\n")
        temp_file.close()

    def tearDown(self):
        """Tidy up files created for testing."""

        os.remove(self.fcout)
        os.remove(self.overlapping_gff)
        os.removedirs(self.existing_path)

    def test_no_reads_file(self):
        """ Test parser throws IOError when there is no reads file. """

        read_files = ["Data/does_not_exist.bam"]
        a = GffParser(self.overlapping_gff)
        self.assertRaises(IOError, FeatureCountParser, a, "-p -s 2", "featureCounts", read_files)

    def test_duff_reads_file(self):
        """ Test parser throws IOError when the reads file is not in the right format. """

        read_files = [self.fcout]
        a = GffParser(self.overlapping_gff)
        self.assertRaises(IOError, FeatureCountParser, a, "featureCounts", "-p -s 2", read_files)

    def test_unmatched_chromosomes(self):
        """ Test parser throws IOError when the chromosomes in the annotation and reads file do not match. """

        read_files = [self.test_data_bam]
        a = GffParser(self.overlapping_gff)
        self.assertRaises(IOError, FeatureCountParser, a, "featureCounts", "-p -s 2", read_files)

    def test_bam_without_chromosomes(self):

        """ Test with wacky bam file which only has chromosomes listed in binary part, not in text header """
        read_files = [self.test_unstranded_bam]
        a = GtfParser(self.test_unstranded_gtf)
        fp = FeatureCountParser(a, "featureCounts", "-p -s 2", read_files)

        self.assertIsNotNone(fp)

    def test_parser(self):
        """ Test parser builds counts file. """

        read_files = [self.test_data_bam]
        a = GffParser(self.test_data_gff)
        fp = FeatureCountParser(a, "featureCounts", "-p -s 2", read_files)

        self.assertIsNotNone(fp)

    def test_counts_per_feature(self):
        """ Test counts per feature calculated correctly. """

        read_files = [self.test_data_bam]
        a = GffParser(self.test_data_gff)
        fp = FeatureCountParser(a, "featureCounts", "-p -s 2", read_files)

        counts = fp.counts_per_feature(fp.featureType.exon)
        self.assertItemsEqual(counts,
            [32, 51, 11, 49, 12, 49, 173, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

        counts = fp.counts_per_feature(fp.featureType.intron)
        self.assertItemsEqual(counts,
                             [0, 6, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0])

        counts = fp.counts_per_feature(fp.featureType.three_prime_UTR)
        self.assertItemsEqual(counts, [16, 0, 0, 0, 0, 0, 0, 0])

        counts = fp.counts_per_feature(fp.featureType.five_prime_UTR)
        self.assertItemsEqual(counts, [8, 0, 0, 0, 170, 0])

    def test_reads_in_first_feature(self):
        """ Test number of reads in first feature calculated correctly. """

        read_files = [self.test_data_bam]
        a = GffParser(self.test_data_gff)
        fp = FeatureCountParser(a, "featureCounts", "-p -s 2", read_files)

        counts = fp.reads_in_first_feature(fp.featureType.exon)
        self.assertListEqual(counts["counts"].tolist(), [32, 173, 0, 0, 0, 0, 0])

        counts = fp.reads_in_first_feature(fp.featureType.intron)
        self.assertListEqual(counts["counts"].tolist(), [0, 0, 0, 0, 0, 0])

    def test_reads_in_feature_wrong_feature(self):
        """ Test number of reads in first feature calculated correctly. """

        read_files = [self.test_data_bam]
        a = GffParser(self.test_data_gff)
        fp = FeatureCountParser(a, "featureCounts", "-p -s 2", read_files)

        with self.assertRaises(ValueError):
            fp.reads_in_first_feature("rubbish")

        with self.assertRaises(ValueError):
            fp.reads_in_last_feature("rubbish")

    def test_reads_in_last_feature(self):
        """ Test number of reads in first feature calculated correctly. """

        read_files = [self.test_data_bam]
        a = GffParser(self.test_data_gff)
        fp = FeatureCountParser(a, "featureCounts", "-p -s 2", read_files)

        counts = fp.reads_in_last_feature(fp.featureType.exon)
        self.assertListEqual(counts["counts"].tolist(), [0, 0, 0, 49, 0, 0, 0])

        counts = fp.reads_in_last_feature(fp.featureType.intron)
        self.assertListEqual(counts["counts"].tolist(), [0, 0, 0, 0, 0, 0])

    def test_create_data_by_gene(self):

        """ Test creation of dataframes by gene. """

        read_files = [self.test_data_bam]
        a = GffParser(self.test_data_gff)
        fp = FeatureCountParser(a, "featureCounts", "-p -s 2", read_files)

        exon_lengths = a.get_ftr_lengths(a.featureType.exon, bygene=True)

        total_mapped_reads = 1000

        result = fp.create_data_by_gene([fp.featureType.exon, fp.featureType.intron], total_mapped_reads,
                                        [exon_lengths])

        self.assertListEqual(result[fp.featureType.exon].tolist(), [204, 0, 0, 173, 0])
        self.assertListEqual(result[fp.featureType.intron].tolist()[:4], [7, 0, 0, 4])
        self.assertTrue(math.isnan(result[fp.featureType.intron].tolist()[4]))

    def test_create_data_by_gene_extras_only(self):

        """ Test creation of dataframes by gene when only extras are supplied. """

        read_files = [self.test_data_bam]
        a = GffParser(self.test_data_gff)
        fp = FeatureCountParser(a, "featureCounts", "-p -s 2", read_files)

        exon_lengths = a.get_ftr_lengths(a.featureType.exon, bygene=True)
        gene_names = a.get_gene_names()

        total_mapped_reads = 1000

        result = fp.create_data_by_gene([], total_mapped_reads, [exon_lengths, gene_names])

        self.assertListEqual(result["length"].tolist(), exon_lengths["length"].tolist())
        self.assertListEqual(result["Name"].tolist(), gene_names["Name"].tolist())

    def test_create_data_by_gene_no_data(self):

        """ Test creation of dataframes by gene when only extras are supplied. """

        read_files = [self.test_data_bam]
        a = GffParser(self.test_data_gff)
        fp = FeatureCountParser(a, "featureCounts", "-p -s 2", read_files)

        total_mapped_reads = 1000

        with self.assertRaises(ValueError):
            fp.create_data_by_gene([], total_mapped_reads, [])
