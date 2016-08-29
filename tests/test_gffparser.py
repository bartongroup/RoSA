#!/usr/bin/python
"""
Unit tests for gffparser module
"""

import unittest
import tempfile
import os.path
import pandas as pd

from viewseq.gffparser import GffParser
from viewseq.constants import *

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"


class TestGffParser(unittest.TestCase):

    """
    Test functions in parser module.
    """

    def setUp(self):

        """Initialize the framework for testing.

        Define and create a new system temporary directory with the
        `tempfile <http://docs.python.org/2/library/tempfile.html>`_ package
        for use as an existing directory.
        Creates gff files in this directory to test parsing functionality.
        """

        # setup a new temp directory
        self.existing_path = tempfile.mkdtemp()

        # output path for gtf2 file
        self.output_gtf = "{}/gfftogtf.gtf".format(self.existing_path)

        # sequence ontology info
        self.gene_ontology = "Data/gene-SO.txt"
        self.exon_ontology = "Data/exon-SO.txt"

        # make the sequence ontology example gff v3 file
        self.seqont_gff = "{}/seqont_gff_v3.gff".format(self.existing_path)
        temp_file = open(self.seqont_gff, "w")
        temp_file.write("##gff-version 3\n##sequence-region\t\tctg123\t1\t1497228\n")
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

        # make the sequence ontology example gff v3 file with a different gene name
        self.seqont2_gff = "{}/seqont2_gff_v3.gff".format(self.existing_path)
        temp_file = open(self.seqont2_gff, "w")
        temp_file.write("##gff-version 3\n##sequence-region\t\tctg123\t1\t1497228\n")
        temp_file.write("ctg123\t.\tpseudogene\t1000\t9000\t.\t+\t.\tID=gene00001;Name=EDEN\n")
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
        temp_file.write("##gff-version 3\n##sequence-region\t\tctg123\t1\t1497228\n")
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
        temp_file.write("ctg123\t.\tfive_prime_UTR\t1300\t1500\t.\t+\t0\tParent=mRNA00003\n")
        temp_file.write("ctg123\t.\tfive_prime_UTR\t3000\t3300\t.\t+\t0\tParent=mRNA00003\n")
        temp_file.write("ctg123\t.\tfive_prime_UTR\t3000\t3390\t.\t+\t0\tParent=mRNA00003\n")
        temp_file.write("ctg123\t.\tthree_prime_UTR\t7601\t9000\t.\t+\t0\tParent=mRNA00001,mRNA00002,mRNA00003\n")
        temp_file.close()

        # a gff with overlapping +ve and -ve strand exons
        self.overlapping_gff = "{}/overlapping_gff_v3.gff".format(self.existing_path)
        temp_file = open(self.overlapping_gff, "w")
        temp_file.write("##gff-version 3\n##sequence-region\t\tctg123\t1\t1497228\n")
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
        temp_file.write("Chr1\t.\tmRNA\t8900\t12000\t.\t-\t.\tID=mRNA00007;Parent=gene00002\n")
        temp_file.write("Chr1\t.\texon\t5590\t7000\t.\t-\t.\tID=exon00006;Parent=mRNA00006\n")
        temp_file.write("Chr1\t.\texon\t8900\t12000\t.\t-\t.\tID=exon00007;Parent=mRNA00006,mRNA00007\n")
        temp_file.write("Chr1\t.\tthree_prime_UTR\t5590\t6000\t.\t-\t.\tParent=mRNA00006\n")
        temp_file.write("Chr1\t.\tthree_prime_UTR\t8900\t9200\t.\t-\t.\tParent=mRNA00007\n")
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

        # a gff with CDSs and -ve strand exons
        self.negstrandcds_gff = "{}/negstrandcds_gff_v3.gff".format(self.existing_path)
        temp_file = open(self.negstrandcds_gff, "w")
        temp_file.write("##gff-version 3\n##sequence-region\t\tctg123\t1\t1497228\n")
        temp_file.write("Chr1\t.\tgene\t5550\t12000\t.\t-\t.\tID=gene00002;Name=OPP\n")
        temp_file.write("Chr1\t.\tmRNA\t5590\t12000\t.\t-\t.\tID=mRNA00006;Parent=gene00002\n")
        temp_file.write("Chr1\t.\tmRNA\t9000\t12000\t.\t-\t.\tID=mRNA00007;Parent=gene00002\n")
        temp_file.write("Chr1\t.\texon\t5590\t7000\t.\t-\t.\tID=exon00006;Parent=mRNA00006\n")
        temp_file.write("Chr1\t.\texon\t9000\t12000\t.\t-\t.\tID=exon00007;Parent=mRNA00006,mRNA00007\n")
        temp_file.write("Chr1\t.\tCDS\t9000\t11800\t.\t-\t0\tID=cds00001;Parent=mRNA00006;Name=protein1\n")
        temp_file.write("Chr1\t.\tCDS\t6001\t7000\t.\t-\t0\tID=cds00001;Parent=mRNA00006;Name=protein1\n")
        temp_file.write("Chr1\t.\tCDS\t9201\t11800\t.\t-\t0\tID=cds00002;Parent=mRNA00007;Name=protein2\n")
        temp_file.close()

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

        # make a gff file with a different exon name
        self.namechange_gff = "{}/namechange.gff".format(self.existing_path)
        temp_file = open(self.namechange_gff, "w")
        temp_file.write("##gff-version 3\n# These are \n# header\n# lines\n")
        temp_file.write("chrI\t.\tchromosome\t1\t10000\t.\t.\t.\t.\n")
        temp_file.write("chrI\t.\tdiffexon\t1000\t2000\t.\t-\t.\t.\n")
        temp_file.write("chrI\t.\tgene\t900\t2200\t.\t-\t.\t.\n")
        temp_file.write("chrI\t.\tmRNA\t1000\t2000\t.\t.-\t.\t.\n")
        temp_file.write("chr2\t.\tchromosome\t1\t20000\t.\t.\t.\t.\n")
        temp_file.write("chr2\t.\tdiffexon\t2000\t4000\t.\t+\t.\t.\n")
        temp_file.write("chr2\t.\tgene\t1800\t4400\t.\t+\t.\t.\n")
        temp_file.write("chr2\t.\tdiffexon\t4100\t4300\t.\t+\t.\t.\n")
        temp_file.write("chr3\t.\tchromosome\t1\t15000\t.\t.\t.\t.\n")
        temp_file.write("chr3\t.\tdiffexon\t5000\t5450\t.\t+\t.\t.\n")
        temp_file.write("chr3\t.\tgene\t5100\t5300\t.\t+\t.\t.\n")
        temp_file.close()

        # gtf2 version of overlapping.gff for checking gtf output
        self.test_gtf = "{}/overlapping.gtf".format(self.existing_path)
        temp_file = open(self.test_gtf, "w")

        temp_file.write('Chr1\t.\tintron\t1501\t2999\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n')
        temp_file.write('Chr1\t.\tintron\t3903\t4999\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n')
        temp_file.write('Chr1\t.\tintron\t5501\t6999\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n')
        temp_file.write('Chr1\t.\tintron\t1501\t4999\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00002\";\n')
        temp_file.write('Chr1\t.\tintron\t5501\t6999\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00002\";\n')
        temp_file.write('Chr1\t.\tintron\t1501\t2999\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n')
        temp_file.write('Chr1\t.\tintron\t3903\t4999\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n')
        temp_file.write('Chr1\t.\tintron\t5501\t6999\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n')
        temp_file.write('Chr2\t.\tintron\t1501\t2999\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00008\";\n')
        temp_file.write('Chr2\t.\tintron\t3903\t4999\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00008\";\n')
        temp_file.write('Chr2\t.\tintron\t5501\t6999\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00008\";\n')
        temp_file.write('Chr2\t.\tintron\t1501\t4999\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00009\";\n')
        temp_file.write('Chr2\t.\tintron\t5501\t6999\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00009\";\n')
        temp_file.write('Chr2\t.\tintron\t1501\t2999\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00010\";\n')
        temp_file.write('Chr2\t.\tintron\t3903\t4999\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00010\";\n')
        temp_file.write('Chr2\t.\tintron\t5501\t6999\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00010\";\n')
        temp_file.write('Chr1\t.\tintron\t7001\t8899\t.\t-\t.\tgene_id \"gene00002\"; transcript_id \"mRNA00006\";\n')
        temp_file.write("Chr1\t.\texon\t1050\t1500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("Chr1\t.\texon\t3000\t3902\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("Chr1\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("Chr1\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("Chr1\t.\texon\t1050\t1500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("Chr1\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("Chr1\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("Chr1\t.\texon\t1300\t1500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("Chr1\t.\texon\t3000\t3902\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("Chr1\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("Chr1\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")

        temp_file.write("Chr2\t.\texon\t1050\t1500\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00008\";\n")
        temp_file.write("Chr2\t.\texon\t3000\t3902\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00008\";\n")
        temp_file.write("Chr2\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00008\";\n")
        temp_file.write("Chr2\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00008\";\n")
        temp_file.write("Chr2\t.\texon\t1050\t1500\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00009\";\n")
        temp_file.write("Chr2\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00009\";\n")
        temp_file.write("Chr2\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00009\";\n")
        temp_file.write("Chr2\t.\texon\t1300\t1500\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00010\";\n")
        temp_file.write("Chr2\t.\texon\t3000\t3902\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00010\";\n")
        temp_file.write("Chr2\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00010\";\n")
        temp_file.write("Chr2\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00010\";\n")

        temp_file.write("Chr1\t.\texon\t5590\t7000\t.\t-\t.\tgene_id \"gene00002\"; transcript_id \"mRNA00006\";\n")
        temp_file.write("Chr1\t.\texon\t8900\t12000\t.\t-\t.\tgene_id \"gene00002\"; transcript_id \"mRNA00006\";\n")
        temp_file.write("Chr1\t.\texon\t8900\t12000\t.\t-\t.\tgene_id \"gene00002\"; transcript_id \"mRNA00007\";\n")

        temp_file.close()

        self.test_fc = "{}/counts.out".format(self.existing_path)
        temp_file = open(self.test_fc, "w")
        temp_file.write("Geneid\tChr\tStart\tEnd\tStrand\tLength\tmini.bam\n")
        temp_file.write("AT1G01010\t1\t3631\t3713\t+\t83\t10\n")
        temp_file.write("AT1G01010\t1\t3631\t3913\t+\t283\t26\n")
        temp_file.write("AT1G01010\t1\t3996\t4276\t+\t281\t44\n")
        temp_file.write("AT1G01010\t1\t4486\t4605\t+\t120\t9\n")
        temp_file.write("AT1G01010\t1\t4706\t5095\t+\t390\t53\n")
        temp_file.write("AT1G01010\t1\t5174\t5326\t+\t153\t13\n")
        temp_file.write("AT1G01010\t1\t5439\t5899\t+\t461\t57\n")
        temp_file.write("AT1G01040\t1\t23146\t24451\t+\t1306\t102\n")
        temp_file.write("AT1G01040\t1\t24542\t24655\t+\t114\t0\n")
        temp_file.write("AT1G01040\t1\t27618\t27713\t+\t 96\t0\n")
        temp_file.write("AT1G01040\t1\t27803\t28431\t+\t629\t0\n")
        temp_file.write("AT1G01040\t1\t27372\t27536\t+\t165\t23\n")
        temp_file.write("AT1G01040\t1\t30410\t30816\t+\t407\t4\n")
        temp_file.write("AT1G01040\t1\t30902\t31120\t+\t219\t77\n")
        temp_file.write("AT1G01020\t1\t5928\t6263\t-\t336\t12\n")
        temp_file.write("AT1G01020\t1\t6437\t7069\t-\t633\t0\n")
        temp_file.write("AT1G01020\t1\t7157\t7232\t-\t76\t1\n")
        temp_file.write("AT1G01020\t1\t7384\t7450\t-\t67\t2\n")
        temp_file.write("AT1G01020\t1\t8571\t8737\t-\t86\t3\n")
        temp_file.write("AT1G01030\t2\t11649\t13173\t-\t1525\t4\n")
        temp_file.write("AT1G01030\t2\t13335\t13714\t-\t380\t5\n")
        temp_file.write("AT1G01030\t2\t14649\t15173\t-\t524\t6\n")
        temp_file.write("AT1G01030\t2\t16335\t17714\t-\t1379\t7\n")
        temp_file.write("AT1G01050\t2\t18649\t19222\t-\t573\t8\n")
        temp_file.write("AT1G01050\t2\t20145\t20604\t-\t459\t10\n")
        temp_file.write("AT1G01050\t2\t21998\t22313\t-\t315\t11\n")
        temp_file.write("AT1G01050\t2\t22467\t22811\t-\t344\t13\n")
        temp_file.close()

        self.test_fc_t_utr = "{}/tcounts.out".format(self.existing_path)
        temp_file = open(self.test_fc_t_utr, "w")
        temp_file.write("Geneid\tChr\tStart\tEnd\tStrand\tLength\tmini.bam\n")
        temp_file.write("AT1G01010\t1\t5631\t5899\t+\t268\t99\n")
        temp_file.write("AT1G01020\t1\t6437\t6914\t+\t477\t88\n")
        temp_file.write("AT1G01020\t1\t5928\t6263\t+\t335\t77\n")

        temp_file.close()

        self.test_fc_gff = "{}/counts.gff".format(self.existing_path)
        temp_file = open(self.test_fc_gff, "w")
        temp_file.write("1\tTAIR10\tchromosome\t1\t30427671\t.\t.\t.\tID=1;Name=1\n")
        temp_file.write("1\tTAIR10\tgene\t3631\t5899\t.\t+\t.\tID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010\n")
        temp_file.write("1\tTAIR10\tmRNA\t3631\t5899\t.\t+\t.\t"
                        "ID=AT1G01010.1;Parent=AT1G01010;Name=AT1G01010.1;Index=1\n")
        temp_file.write("1\tTAIR10\texon\t3631\t3913\t.\t+\t.\tParent=AT1G01010.1\n")
        temp_file.write("1\tTAIR10\tfive_prime_UTR\t3631\t3759\t.\t+\t.\tParent=AT1G01010.1\n")
        temp_file.write("1\tTAIR10\texon\t3996\t4276\t.\t+\t.\tParent=AT1G01010.1\n")
        temp_file.write("1\tTAIR10\texon\t4486\t4605\t.\t+\t.\tParent=AT1G01010.1\n")
        temp_file.write("1\tTAIR10\texon\t4706\t5095\t.\t+\t.\tParent=AT1G01010.1\n")
        temp_file.write("1\tTAIR10\texon\t5174\t5326\t.\t+\t.\tParent=AT1G01010.1\n")
        temp_file.write("1\tTAIR10\texon\t5439\t5899\t.\t+\t.\tParent=AT1G01010.1\n")
        temp_file.write("1\tTAIR10\tthree_prime_UTR\t5631\t5899\t.\t+\t.\tParent=AT1G01010.1\n")
        temp_file.write("1\tTAIR10\tgene\t5928\t8737\t.\t-\t.\tID=AT1G01020;Note=protein_coding_gene;Name=AT1G01020\n")
        temp_file.write("1\tTAIR10\tmRNA\t5928\t8737\t.\t-\t.\t"
                        "ID=AT1G01020.1;Parent=AT1G01020;Name=AT1G01020.1;Index=1\n")
        temp_file.write("1\tTAIR10\tfive_prime_UTR\t8667\t8737\t.\t-\t.\tParent=AT1G01020.1\n")
        temp_file.write("1\tTAIR10\texon\t8571\t8737\t.\t-\t.\tParent=AT1G01020.1\n")
        temp_file.write("1\tTAIR10\texon\t8417\t8464\t.\t-\t.\tParent=AT1G01020.1\n")
        temp_file.write("1\tTAIR10\texon\t8236\t8325\t.\t-\t.\tParent=AT1G01020.1\n")
        temp_file.write("1\tTAIR10\texon\t7942\t7987\t.\t-\t.\tParent=AT1G01020.1\n")
        temp_file.write("1\tTAIR10\texon\t7762\t7835\t.\t-\t.\tParent=AT1G01020.1\n")
        temp_file.write("1\tTAIR10\texon\t7564\t7649\t.\t-\t.\tParent=AT1G01020.1\n")
        temp_file.write("1\tTAIR10\texon\t7384\t7450\t.\t-\t.\tParent=AT1G01020.1\n")
        temp_file.write("1\tTAIR10\texon\t7157\t7232\t.\t-\t.\tParent=AT1G01020.1\n")
        temp_file.write("1\tTAIR10\tthree_prime_UTR\t6437\t6914\t.\t-\t.\tParent=AT1G01020.1\n")
        temp_file.write("1\tTAIR10\texon\t6437\t7069\t.\t-\t.\tParent=AT1G01020.1\n")
        temp_file.write("1\tTAIR10\tthree_prime_UTR\t5928\t6263\t.\t-\t.\tParent=AT1G01020.1\n")
        temp_file.write("1\tTAIR10\texon\t5928\t6263\t.\t-\t.\tParent=AT1G01020.1\n")

        temp_file.write("1\tTAIR10\tgene\t23146\t31120\t.\t+\t.\tID=AT1G01040;"
                        "Note=protein_coding_gene;Name=AT1G01040\n")
        temp_file.write("1\tTAIR10\tmRNA\t23146\t31120\t.\t+\t.\tID=AT1G01040.1;"
                        "Parent=AT1G01040;Name=AT1G01040.1;Index=1\n")
        temp_file.write("1\tTAIR10\texon\t23146\t24451\t.\t+\t.\tParent=AT1G01040.1\n")
        temp_file.write("1\tTAIR10\texon\t24542\t24655\t.\t+\t.\tParent=AT1G01040.1\n")
        temp_file.write("1\tTAIR10\texon\t27618\t27713\t.\t+\t.\tParent=AT1G01040.1\n")
        temp_file.write("1\tTAIR10\texon\t27803\t28431\t.\t+\t.\tParent=AT1G01040.1\n")
        temp_file.write("1\tTAIR10\texon\t27372\t27536\t.\t+\t.\tParent=AT1G01040.1\n")
        temp_file.write("1\tTAIR10\texon\t30410\t30816\t.\t+\t.\tParent=AT1G01040.1\n")
        temp_file.write("1\tTAIR10\texon\t30902\t31120\t.\t+\t.\tParent=AT1G01040.1\n")

        temp_file.write("2\tTAIR10\tgene\t11649\t17714\t.\t+\t.\tID=AT1G01030;"
                        "Note=protein_coding_gene;Name=AT1G01030\n")
        temp_file.write("2\tTAIR10\tmRNA\t11649\t17714\t.\t+\t.\tID=AT1G01030.1;"
                        "Parent=AT1G01030;Name=AT1G01030.1;Index=1\n")
        temp_file.write("2\tTAIR10\tchromosome\t1\t30427671\t.\t.\t.\tID=2;Name=2\n")
        temp_file.write("2\tTAIR10\texon\t11649\t13173\t.\t+\t.\tParent=AT1G01030.1\n")
        temp_file.write("2\tTAIR10\texon\t13335\t13714\t.\t+\t.\tParent=AT1G01030.1\n")
        temp_file.write("2\tTAIR10\texon\t14649\t15173\t.\t+\t.\tParent=AT1G01030.1\n")
        temp_file.write("2\tTAIR10\texon\t16335\t17714\t.\t+\t.\tParent=AT1G01030.1\n")

        temp_file.write("1\tTAIR10\tgene\t18649\t22811\t.\t+\t.\tID=AT1G01050;"
                        "Note=protein_coding_gene;Name=AT1G01050\n")
        temp_file.write("1\tTAIR10\tmRNA\t18649\t22811\t.\t+\t.\tID=AT1G01050.1;"
                        "Parent=AT1G01050;Name=AT1G01050.1;Index=1\n")
        temp_file.write("2\tTAIR10\texon\t18649\t19222\t.\t+\t.\tParent=AT1G01050.1\n")
        temp_file.write("2\tTAIR10\texon\t20145\t20604\t.\t+\t.\tParent=AT1G01050.1\n")
        temp_file.write("2\tTAIR10\texon\t21998\t22313\t.\t+\t.\tParent=AT1G01050.1\n")
        temp_file.write("2\tTAIR10\texon\t22467\t22811\t.\t+\t.\tParent=AT1G01050.1\n")

        temp_file.close()

        # a gff with overlapping +ve and -ve strand exons
        self.antisense_gff = "{}/antisense_gff_v3.gff".format(self.existing_path)
        temp_file = open(self.antisense_gff, "w")
        temp_file.write("##gff-version 3\n##sequence-region\t\tctg123\t1\t1497228\n")
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
        temp_file.write("Chr1\t.\tmRNA\t8900\t12000\t.\t-\t.\tID=mRNA00007;Parent=gene00002\n")
        temp_file.write("Chr1\t.\texon\t5590\t7000\t.\t-\t.\tID=exon00006;Parent=mRNA00006\n")
        temp_file.write("Chr1\t.\texon\t7100\t7400\t.\t-\t.\tID=exon00007;Parent=mRNA00006\n")
        temp_file.write("Chr1\t.\texon\t7500\t8000\t.\t-\t.\tID=exon00008;Parent=mRNA00006\n")
        temp_file.write("Chr1\t.\texon\t8900\t12000\t.\t-\t.\tID=exon00009;Parent=mRNA00006,mRNA00007\n")
        temp_file.write("Chr1\t.\tthree_prime_UTR\t5590\t6000\t.\t-\t.\tParent=mRNA00006\n")
        temp_file.write("Chr1\t.\tthree_prime_UTR\t9000\t9200\t.\t-\t.\tParent=mRNA00007\n")
        temp_file.write("Chr1\t.\tfive_prime_UTR\t11801\t12000\t.\t-\t.\tParent=mRNA00006,mRNA00007\n")

        temp_file.write("Chr1\t.\tgene\t1300\t3800\t.\t-\t.\tID=gene00004;Name=ANTI\n")
        temp_file.write("Chr1\t.\tmRNA\t1350\t3800\t.\t-\t.\tID=mRNA00011;Parent=gene00004\n")
        temp_file.write("Chr1\t.\texon\t1350\t1480\t.\t-\t.\tID=exon00010;Parent=mRNA00011\n")
        temp_file.write("Chr1\t.\texon\t2500\t3800\t.\t-\t.\tID=exon00011;Parent=mRNA00011\n")
        temp_file.write("Chr1\t.\tthree_prime_UTR\t1350\t1420\t.\t-\t.\tParent=mRNA00011\n")
        temp_file.write("Chr1\t.\tfive_prime_UTR\t3400\t3800\t.\t-\t.\tParent=mRNA00011\n")

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

        os.remove(self.utrseqont_gff)
        os.remove(self.seqont_gff)
        os.remove(self.seqont2_gff)
        os.remove(self.simple_gff)
        os.remove(self.overlapping_gff)
        os.remove(self.antisense_gff)
        os.remove(self.negstrandcds_gff)
        os.remove(self.test_gtf)
        os.remove(self.namechange_gff)
        os.remove(self.test_fc_gff)
        os.remove(self.test_fc)
        os.remove(self.test_fc_t_utr)
        os.removedirs(self.existing_path)

    # def test_are_unexpected_names(self):
    #     """ Test that user is warned about missing UTR names. """
    #
    #     # catch user warning that UTR names are missing
    #     with warnings.catch_warnings(record=True) as w:
    #
    #         GffParser(self.simple_gff)
    #         self.assertTrue(len(w) == 1)
    #         self.assertTrue(issubclass(w[-1].category, UserWarning))
    #         self.assertTrue("Some expected type names are missing from the annotation input file"
    #                         in str(w[-1].message))
    #
    #     # Note: there is some strange interaction with this test and test_exons_per_gene_no_geneIDs
    #     # If this test is run after test_exons_per_gene_no_geneIDs then it fails, with w=[]
    #     # The failure seems to be something to do with test_exons_per_gene_no_geneIDs also
    #     # causing this UserWarning to be raised, but as yet I don't see how to fix it -
    #     # so for now have set this test to run before test_exons_per_gene_no_geneIDs

    def test_intron_lengths(self):
        """ Test intron lengths calc. """

        p = GffParser(self.utrseqont_gff)
        lengths = p.get_ftr_lengths(p.featureType.intron)[LENGTH].tolist()
        self.assertItemsEqual(lengths, [1499, 1499, 3499, 1097])

    def test_exon_lengths(self):
        """ Test exon lengths calc. """

        p = GffParser(self.utrseqont_gff)
        lengths = p.get_ftr_lengths(p.featureType.exon)[LENGTH].tolist()
        self.assertItemsEqual(lengths, [451, 201, 903, 903, 903, 501, 2001])

    def test_exons_per_gene_no_geneIDs(self):
        """ Test empty list is returned when there are no gene IDs available. """

        p = GffParser(self.simple_gff)
        exons_by_gene = p.get_exons_per_gene()[EXON_COUNTS].tolist()
        self.assertItemsEqual(exons_by_gene, [])

    def test_exons_by_gene(self):
        """ Test exon counts calc. """

        p = GffParser(self.utrseqont_gff)
        exons_by_gene = p.get_exons_per_gene()[EXON_COUNTS].tolist()
        self.assertItemsEqual(exons_by_gene, [7])

    def test_transcripts_by_gene(self):
        """ Test transcript counts calc. """

        p = GffParser(self.overlapping_gff)
        ts_by_gene = p.get_transcripts_per_gene()[TRANSCRIPT_COUNTS].tolist()
        self.assertItemsEqual(ts_by_gene, [3, 2, 3])

    def test_overlapping_genes_exon_lengths(self):
        """ Test exon lengths when genes in +ve and -ve strands overlap """

        p = GffParser(self.overlapping_gff)
        lengths = p.get_ftr_lengths(p.featureType.exon)[LENGTH].tolist()
        self.assertItemsEqual(lengths, [451, 201, 903, 903, 903, 501, 2001, 451, 201, 903, 903, 903, 501, 2001,
                                        1411, 3101, 3101])

    def test_overlapping_genes_exon_count(self):
        """ Test exon counts when genes in +ve and -ve strands overlap """

        p = GffParser(self.overlapping_gff)
        exons_by_gene = p.get_exons_per_gene()[EXON_COUNTS].tolist()
        self.assertItemsEqual(exons_by_gene, [7, 3, 7])

    def test_overlapping_genes_intron_lengths(self):
        """ Test intron lengths when genes in +ve and -ve strands overlap """

        p = GffParser(self.overlapping_gff)
        lengths = p.get_ftr_lengths(p.featureType.intron)[LENGTH].tolist()
        self.assertItemsEqual(lengths, [1499, 1499, 3499, 1097, 1499, 1499, 3499, 1097, 1899])

    def test_missing_names(self):
        """ Test that user is notified about unexpected feature names. """

        self.assertRaises(ValueError, GffParser, self.namechange_gff)

    def test_output_to_gtf2(self):
        """ Test gff input correctly output to gtf """

        p = GffParser(self.overlapping_gff)
        p.export_to_gtf2(os.path.splitext(self.output_gtf)[0])

        # check file exists
        self.assertTrue(os.path.exists(self.output_gtf), "Output gtf file does not exist")

        # check contents correspond to gtf2 format
        new_gtf_handle = open(self.output_gtf, "r")
        test_gtf_handle = open(self.test_gtf, "r")

        newdata = new_gtf_handle.readlines()
        testdata = test_gtf_handle.readlines()

        self.assertListEqual(newdata, testdata)

        new_gtf_handle.close()
        test_gtf_handle.close()

        os.remove(self.output_gtf)

    def test_split_output_to_gtf2(self):
        """ Test gff input correctly output to into split gtf files """

        p = GffParser(self.overlapping_gff)
        p.export_to_gtf2(os.path.splitext(self.output_gtf)[0], split=True)

        # check file exists
        files = ["{}Chr1+.gtf".format(os.path.splitext(self.output_gtf)[0]),
                 "{}Chr1-.gtf".format(os.path.splitext(self.output_gtf)[0]),
                 "{}Chr2+.gtf".format(os.path.splitext(self.output_gtf)[0]),
                 "{}Chr2-.gtf".format(os.path.splitext(self.output_gtf)[0])]

        c = 1
        s = "+"
        for f in files:
            self.assertTrue(os.path.exists(f),
                            "Output gtf file does not exist")

            if f == "{}Chr2-.gtf".format(os.path.splitext(self.output_gtf)[0]):
                # this file is empty for this data
                pass
            else:

                # check contents are only correct chromosome and strand
                # assume gtf-ness is checked by earlier test
                data = pd.read_csv(f, header=None, sep="\t")
                chrs = pd.unique(data[[0]].values)
                self.assertEqual(len(chrs[0]), 1)
                self.assertEqual(chrs[0][0], "Chr{}".format(c))

                strand = pd.unique(data[[6]].values)
                self.assertEqual(len(strand[0]), 1)
                self.assertEqual(strand[0][0], s)

            if s == "+":
                s = "-"
            elif s == "-":
                s = "+"
                c = 2

            os.remove(f)

    def test_first_feature_filter(self):
        """ Test filtering by first exon works correctly. """

        data = pd.read_csv(self.test_fc, header=0, sep="\t")
        p = GffParser(self.test_fc_gff)

        filtered = p.filter_by_first_features(data, p.featureType.exon, "Start", "End")
        self.assertEqual(len(filtered.index), 5)
        self.assertListEqual(filtered['mini.bam'].tolist(), [26, 102, 3, 7, 13])

    def test_last_feature_filter(self):
        """ Test filtering by first exon works correctly. """

        data = pd.read_csv(self.test_fc, header=0, sep="\t")
        p = GffParser(self.test_fc_gff)
        filtered = p.filter_by_last_features(data, p.featureType.exon, "Start", "End")

        self.assertEqual(len(filtered.index), 5)

        self.assertListEqual(filtered['mini.bam'].tolist(), [12, 4, 8, 57, 77])

    def test_cds_parsing(self):
        """ Test parsing when UTRs are derived from CDSs """

        p = GffParser(self.seqont_gff)
        p2 = GffParser(self.utrseqont_gff)

        self.assertListEqual(p.get_ftr_lengths(p.featureType.exon)[LENGTH].tolist(),
                             p2.get_ftr_lengths(p.featureType.exon)[LENGTH].tolist())
        self.assertListEqual(p.get_ftr_lengths(p.featureType.intron)[LENGTH].tolist(),
                             p2.get_ftr_lengths(p.featureType.intron)[LENGTH].tolist())
        self.assertListEqual(p.get_exons_per_gene()[EXON_COUNTS].tolist(),
                             p2.get_exons_per_gene()[EXON_COUNTS].tolist())

        data = pd.read_csv(self.test_fc, header=0, sep="\t")
        filtered = p.filter_by_first_features(data, p.featureType.exon, "Start", "End")
        filtered2 = p2.filter_by_first_features(data, p2.featureType.exon, "Start", "End")
        self.assertEqual(len(filtered.index), len(filtered2.index))
        self.assertListEqual(filtered['mini.bam'].tolist(), filtered2['mini.bam'].tolist())

        filtered = p.filter_by_last_features(data, p.featureType.exon, "Start", "End")
        filtered2 = p2.filter_by_last_features(data, p2.featureType.exon, "Start", "End")
        self.assertEqual(len(filtered.index), len(filtered2.index))
        self.assertListEqual(filtered['mini.bam'].tolist(), filtered2['mini.bam'].tolist())

    def test_cds_negstrand(self):
        """ Test CDS parsing of gene on negative strand. """

        p = GffParser(self.negstrandcds_gff)

        self.assertListEqual(p.get_ftr_lengths(p.featureType.exon)[LENGTH].tolist(), [1411, 3001, 3001])
        self.assertListEqual(p.get_ftr_lengths(p.featureType.intron)[LENGTH].tolist(), [1999])

    def test_get_data_by_feature(self):
        """ Test data retrieval by feature type. """

        p = GffParser(self.utrseqont_gff)

        results = p.get_data_by_feature()
        for name, df in results:
            if name == p.featureType.intron:
                lengths = df[LENGTH].tolist()
                self.assertItemsEqual(lengths, [1499, 1499, 3499, 1097])
            elif name == p.featureType.exon:
                lengths = df[LENGTH].tolist()
                self.assertItemsEqual(lengths, [451, 201, 903, 903, 903, 501, 2001])

    def test_gene_ontologies(self):
        """ Test reading of ontology data. """
        p = GffParser(self.seqont2_gff, self.gene_ontology)
        self.assertItemsEqual(p.get_ftr_lengths(p.featureType.exon)[LENGTH].tolist(),
                                [451, 201, 903, 903, 903, 501, 2001])

    def test_build_antisense(self):
        """ Test transcript counts calc. """

        antisense_gtf = "antisense"
        p = GffParser(self.antisense_gff)
        p.build_antisense_gtf([(p.featureType.exon, antisense_gtf)])

        # check file exists
        fullfile = "{}.gtf".format(antisense_gtf)
        self.assertTrue(os.path.exists(fullfile), "Output gtf file does not exist")

        result = pd.read_csv(fullfile, header=None, comment='#', sep='\t', engine='python')
        result.columns = [CHROMOSOME, SOURCE, TYPE, START, STOP, SCORE, STRAND, PHASE, ATTRIBUTES]

        # check we have the correct number of -ve and +ve segments
        self.assertEqual(len(result[result[STRAND] == '-']), 4)
        self.assertEqual(len(result[result[STRAND] == '+']), 31)

        os.remove(fullfile)

    def test_build_antisense_gtf_gene_only(self):
        """ Test transcript counts calc. """

        antisense_gtf = "antisense"
        p = GffParser(self.antisense_gff)
        p.build_antisense_gtf_gene_only(antisense_gtf)

        # check file exists
        fullfile = "{}.gtf".format(antisense_gtf)
        self.assertTrue(os.path.exists(fullfile), "Output gtf file does not exist")

        result = pd.read_csv(fullfile, header=None, comment='#', sep='\t', engine='python')
        result.columns = [CHROMOSOME, SOURCE, TYPE, START, STOP, SCORE, STRAND, PHASE, ATTRIBUTES]

        # check we have the correct number of -ve and +ve segments
        self.assertListEqual(result[START].tolist(), [9001, 1050, 3801, 1050])
        self.assertListEqual(result[STOP].tolist(), [12000, 1349, 5589, 9000])

        os.remove(fullfile)

    def test_utr_parsing(self):
        """ Test parsing when UTRs are derived from CDSs """

        p = GffParser(self.utrseqont_gff)

        self.assertEqual(len(p.cds_tr_join), 10)
        self.assertItemsEqual(p.cds_tr_join[START].tolist(), [1201, 3000, 5000, 7000, 1201,
                                                              5000, 7000, 3391, 5000, 7000])
        self.assertItemsEqual(p.cds_tr_join[STOP].tolist(), [1500, 3902, 5500, 7600, 1500,
                                                             5500, 7600, 3902, 5500, 7600])


if __name__ == '__main__':
    unittest.main(verbosity=2)
