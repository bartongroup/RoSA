#!/usr/bin/env python

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
Unit tests for gffparser module
"""

import unittest
import tempfile
import os.path
import pandas as pd

from antisense.gtfparser import GtfParser
from antisense.constants import *

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"


class TestGtfParser(unittest.TestCase):

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

        # sequence ontology  info
        self.gene_ontology = "Data/gene-SO.txt"
        self.exon_ontology = "Data/exon-SO.txt"

        # output path for temp gtf files
        self.temp_gtf = "{}/temp.gtf".format(self.existing_path)
        self.temp2_gtf = "{}/temp2.gtf".format(self.existing_path)

        # make the sequence ontology example gtf file
        self.seqont_gtf = "{}/seqont.gtf".format(self.existing_path)
        temp_file = open(self.seqont_gtf, "w")
        temp_file.write("ctg123\t.\texon\t1300\t1500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("ctg123\t.\texon\t1050\t1500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("ctg123\t.\texon\t1050\t1500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("ctg123\t.\texon\t3000\t3902\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("ctg123\t.\texon\t3000\t3902\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("ctg123\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("ctg123\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("ctg123\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("ctg123\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("ctg123\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("ctg123\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("ctg123\t.\tCDS\t3000\t3902\t.\t+\t0\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("ctg123\t.\tCDS\t5000\t5500\t.\t+\t0\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("ctg123\t.\tCDS\t7000\t7600\t.\t+\t0\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tgene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("ctg123\t.\tCDS\t5000\t5500\t.\t+\t0\tgene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("ctg123\t.\tCDS\t7000\t7600\t.\t+\t0\tgene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("ctg123\t.\tCDS\t3301\t3902\t.\t+\t0\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("ctg123\t.\tCDS\t5000\t5500\t.\t+\t1\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("ctg123\t.\tCDS\t7000\t7600\t.\t+\t1\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("ctg123\t.\tCDS\t3391\t3902\t.\t+\t0\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("ctg123\t.\tCDS\t5000\t5500\t.\t+\t1\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("ctg123\t.\tCDS\t7000\t7600\t.\t+\t1\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.close()

        # a similar gff file with 3' and 5' UTRs
        self.utrseqont_gtf = "{}/utrseqont.gtf".format(self.existing_path)
        temp_file = open(self.utrseqont_gtf, "w")
        temp_file.write("ctg123\t.\texon\t1300\t1500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("ctg123\t.\texon\t1050\t1500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("ctg123\t.\texon\t1050\t1500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("ctg123\t.\texon\t3000\t3902\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("ctg123\t.\texon\t3000\t3902\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("ctg123\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("ctg123\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("ctg123\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("ctg123\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("ctg123\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("ctg123\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("ctg123\t.\tfive_prime_UTR\t1050\t1200\t.\t+\t0\t"
                        "gene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("ctg123\t.\tfive_prime_UTR\t1050\t1200\t.\t+\t0\t"
                        "gene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("ctg123\t.\tfive_prime_UTR\t1300\t1500\t.\t+\t0\t"
                        "gene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("ctg123\t.\tfive_prime_UTR\t3000\t3300\t.\t+\t0\t"
                        "gene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("ctg123\t.\tfive_prime_UTR\t3000\t3390\t.\t+\t0\t"
                        "gene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("ctg123\t.\tthree_prime_UTR\t7601\t9000\t.\t+\t0\t"
                        "gene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("ctg123\t.\tthree_prime_UTR\t7601\t9000\t.\t+\t0\t"
                        "gene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("ctg123\t.\tthree_prime_UTR\t7601\t9000\t.\t+\t0\t"
                        "gene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.close()

        # a gff with overlapping +ve and -ve strand exons
        self.overlapping_gtf = "{}/overlapping.gtf".format(self.existing_path)
        temp_file = open(self.overlapping_gtf, "w")
        temp_file.write("Chr1\t.\texon\t1300\t1500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\n")
        temp_file.write("Chr1\t.\texon\t1050\t1500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("Chr1\t.\texon\t1050\t1500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("Chr1\t.\texon\t3000\t3902\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("Chr1\t.\texon\t3000\t3902\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("Chr1\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("Chr1\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("Chr1\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("Chr1\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("Chr1\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("Chr1\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("Chr1\t.\tfive_prime_UTR\t1050\t1200\t.\t+\t0\t"
                        "gene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("Chr1\t.\tfive_prime_UTR\t1050\t1200\t.\t+\t0\t"
                        "gene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("Chr1\t.\tfive_prime_UTR\t1300\t1499\t.\t+\t0\t"
                        "gene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("Chr1\t.\tfive_prime_UTR\t3000\t3300\t.\t+\t0\t"
                        "gene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("Chr1\t.\tfive_prime_UTR\t3000\t3390\t.\t+\t0\t"
                        "gene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")
        temp_file.write("Chr1\t.\tthree_prime_UTR\t7601\t9000\t.\t+\t0\t"
                        "gene_id \"gene00001\"; transcript_id \"mRNA00001\";\n")
        temp_file.write("Chr1\t.\tthree_prime_UTR\t7601\t9000\t.\t+\t0\t"
                        "gene_id \"gene00001\"; transcript_id \"mRNA00002\";\n")
        temp_file.write("Chr1\t.\tthree_prime_UTR\t7601\t9000\t.\t+\t0\t"
                        "gene_id \"gene00001\"; transcript_id \"mRNA00003\";\n")

        temp_file.write("Chr1\t.\texon\t5590\t7000\t.\t-\t.\tgene_id \"gene00002\"; transcript_id \"mRNA00006\";\n")
        temp_file.write("Chr1\t.\texon\t9000\t12000\t.\t-\t.\tgene_id \"gene00002\"; transcript_id \"mRNA00006\";\n")
        temp_file.write("Chr1\t.\texon\t9000\t12000\t.\t-\t.\tgene_id \"gene00002\"; transcript_id \"mRNA00007\";\n")
        temp_file.write("Chr1\t.\tthree_prime_UTR\t5590\t6000\t.\t-\t.\t"
                        "gene_id \"gene00002\"; transcript_id \"mRNA00006\";\n")
        temp_file.write("Chr1\t.\tthree_prime_UTR\t9000\t9200\t.\t-\t.\t"
                        "gene_id \"gene00002\"; transcript_id \"mRNA00007\";\n")
        temp_file.write("Chr1\t.\tfive_prime_UTR\t11801\t12000\t.\t-\t.\t"
                        "gene_id \"gene00002\"; transcript_id \"mRNA00006\";\n")
        temp_file.write("Chr1\t.\tfive_prime_UTR\t11801\t12000\t.\t-\t.\t"
                        "gene_id \"gene00002\"; transcript_id \"mRNA00007\";\n")
        
        temp_file.write("Chr2\t.\texon\t1300\t1500\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00010\"\n")
        temp_file.write("Chr2\t.\texon\t1050\t1500\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00008\"\n")
        temp_file.write("Chr2\t.\texon\t1050\t1500\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00009\"\n")
        temp_file.write("Chr2\t.\texon\t3000\t3902\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00008\"\n")
        temp_file.write("Chr2\t.\texon\t3000\t3902\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00010\"\n")
        temp_file.write("Chr2\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00008\"\n")
        temp_file.write("Chr2\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00009\"\n")
        temp_file.write("Chr2\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00010\"\n")
        temp_file.write("Chr2\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00008\"\n")
        temp_file.write("Chr2\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00009\"\n")
        temp_file.write("Chr2\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"mRNA00010\"\n")
        temp_file.write("Chr2\t.\tfive_prime_UTR\t1050\t1200\t.\t+\t0\t"
                        "gene_id \"gene00003\"; transcript_id \"mRNA00008\"\n")
        temp_file.write("Chr2\t.\tfive_prime_UTR\t1050\t1200\t.\t+\t0\t"
                        "gene_id \"gene00003\"; transcript_id \"mRNA00009\"\n")
        temp_file.write("Chr2\t.\tfive_prime_UTR\t1300\t1499\t.\t+\t0\t"
                        "gene_id \"gene00003\"; transcript_id \"mRNA00010\"\n")
        temp_file.write("Chr2\t.\tfive_prime_UTR\t3000\t3300\t.\t+\t0\t"
                        "gene_id \"gene00003\"; transcript_id \"mRNA00010\"\n")
        temp_file.write("Chr2\t.\tfive_prime_UTR\t3000\t3390\t.\t+\t0\t"
                        "gene_id \"gene00003\"; transcript_id \"mRNA00010\"\n")
        temp_file.write("Chr2\t.\tthree_prime_UTR\t7601\t9000\t.\t+\t0\t"
                        "gene_id \"gene00003\"; transcript_id \"mRNA00008\"\n")
        temp_file.write("Chr2\t.\tthree_prime_UTR\t7601\t9000\t.\t+\t0\t"
                        "gene_id \"gene00003\"; transcript_id \"mRNA00009\"\n")
        temp_file.write("Chr2\t.\tthree_prime_UTR\t7601\t9000\t.\t+\t0\t"
                        "gene_id \"gene00003\"; transcript_id \"mRNA00010\"\n")
        temp_file.close()

        # gtf2 version of overlapping.gff for checking gtf output
        self.test_gtf = "{}/testoverlapping.gtf".format(self.existing_path)
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
        temp_file.write('Chr1\t.\tintron\t7001\t8999\t.\t-\t.\tgene_id \"gene00002\"; transcript_id \"mRNA00006\";\n') 
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
        temp_file.write("Chr1\t.\texon\t9000\t12000\t.\t-\t.\tgene_id \"gene00002\"; transcript_id \"mRNA00006\";\n")
        temp_file.write("Chr1\t.\texon\t9000\t12000\t.\t-\t.\tgene_id \"gene00002\"; transcript_id \"mRNA00007\";\n")
        
        temp_file.close()

        # gtf2 version of overlapping.gff for checking gtf output
        self.test_gtf_sq = "{}/overlapping_sq.gtf".format(self.existing_path)
        temp_file = open(self.test_gtf_sq, "w")

        temp_file.write('Chr1\t.\tintron\t1501\t2999\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"gene00001\";\n')
        temp_file.write('Chr1\t.\tintron\t3903\t4999\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"gene00001\";\n')
        temp_file.write('Chr1\t.\tintron\t5501\t6999\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"gene00001\";\n')

        temp_file.write('Chr2\t.\tintron\t1501\t2999\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"gene00003\";\n')
        temp_file.write('Chr2\t.\tintron\t3903\t4999\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"gene00003\";\n')
        temp_file.write('Chr2\t.\tintron\t5501\t6999\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"gene00003\";\n')

        temp_file.write("Chr1\t.\tintron\t7001\t8999\t.\t-\t.\tgene_id \"gene00002\"; transcript_id \"gene00002\";\n")

        temp_file.write("Chr1\t.\texon\t1050\t1500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"gene00001\";\n")
        temp_file.write("Chr1\t.\texon\t3000\t3902\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"gene00001\";\n")
        temp_file.write("Chr1\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"gene00001\";\n")
        temp_file.write("Chr1\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00001\"; transcript_id \"gene00001\";\n")

        temp_file.write("Chr2\t.\texon\t1050\t1500\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"gene00003\";\n")
        temp_file.write("Chr2\t.\texon\t3000\t3902\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"gene00003\";\n")
        temp_file.write("Chr2\t.\texon\t5000\t5500\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"gene00003\";\n")
        temp_file.write("Chr2\t.\texon\t7000\t9000\t.\t+\t.\tgene_id \"gene00003\"; transcript_id \"gene00003\";\n")

        temp_file.write("Chr1\t.\texon\t5590\t7000\t.\t-\t.\tgene_id \"gene00002\"; transcript_id \"gene00002\";\n")
        temp_file.write("Chr1\t.\texon\t9000\t12000\t.\t-\t.\tgene_id \"gene00002\"; transcript_id \"gene00002\";\n")

        temp_file.close()

        self.test_fc = "{}/counts.out".format(self.existing_path)
        temp_file = open(self.test_fc, "w")
        temp_file.write("Geneid\tChr\tStart\tEnd\tStrand\tLength\tmini.bam\n")
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

        self.test_fc_gtf = "{}/counts.gtf".format(self.existing_path)
        temp_file = open(self.test_fc_gtf, "w")
        temp_file.write("1\tTAIR10\texon\t3631\t3913\t.\t+\t.\tgene_id \"AT1G01010\" transcript_id \"AT1G01010.1\"\n")
        temp_file.write("1\tTAIR10\tfive_prime_UTR\t3631\t3759\t.\t+\t.\t"
                        "gene_id \"AT1G01010\" transcript_id \"AT1G01010.1\"\n")
        temp_file.write("1\tTAIR10\texon\t3996\t4276\t.\t+\t.\tgene_id \"AT1G01010\" transcript_id \"AT1G01010.1\"\n")
        temp_file.write("1\tTAIR10\texon\t4486\t4605\t.\t+\t.\tgene_id \"AT1G01010\" transcript_id \"AT1G01010.1\"\n")
        temp_file.write("1\tTAIR10\texon\t4706\t5095\t.\t+\t.\tgene_id \"AT1G01010\" transcript_id \"AT1G01010.1\"\n")
        temp_file.write("1\tTAIR10\texon\t5174\t5326\t.\t+\t.\tgene_id \"AT1G01010\" transcript_id \"AT1G01010.1\"\n")
        temp_file.write("1\tTAIR10\texon\t5439\t5899\t.\t+\t.\tgene_id \"AT1G01010\" transcript_id \"AT1G01010.1\"\n")
        temp_file.write("1\tTAIR10\tthree_prime_UTR\t5631\t5899\t.\t+\t.\t"
                        "gene_id \"AT1G01010\" transcript_id \"AT1G01010.1\"\n")
        
        temp_file.write("1\tTAIR10\tfive_prime_UTR\t8667\t8737\t.\t-\t.\t"
                        "gene_id \"AT1G01020\" transcript_id \"AT1G01020.1\"\n")
        temp_file.write("1\tTAIR10\texon\t8571\t8737\t.\t-\t.\tgene_id \"AT1G01020\" transcript_id \"AT1G01020.1\"\n")
        temp_file.write("1\tTAIR10\texon\t8417\t8464\t.\t-\t.\tgene_id \"AT1G01020\" transcript_id \"AT1G01020.1\"\n")
        temp_file.write("1\tTAIR10\texon\t8236\t8325\t.\t-\t.\tgene_id \"AT1G01020\" transcript_id \"AT1G01020.1\"\n")
        temp_file.write("1\tTAIR10\texon\t7942\t7987\t.\t-\t.\tgene_id \"AT1G01020\" transcript_id \"AT1G01020.1\"\n")
        temp_file.write("1\tTAIR10\texon\t7762\t7835\t.\t-\t.\tgene_id \"AT1G01020\" transcript_id \"AT1G01020.1\"\n")
        temp_file.write("1\tTAIR10\texon\t7564\t7649\t.\t-\t.\tgene_id \"AT1G01020\" transcript_id \"AT1G01020.1\"\n")
        temp_file.write("1\tTAIR10\texon\t7384\t7450\t.\t-\t.\tgene_id \"AT1G01020\" transcript_id \"AT1G01020.1\"\n")
        temp_file.write("1\tTAIR10\texon\t7157\t7232\t.\t-\t.\tgene_id \"AT1G01020\" transcript_id \"AT1G01020.1\"\n")
        temp_file.write("1\tTAIR10\tthree_prime_UTR\t6437\t6914\t.\t-\t.\t"
                        "gene_id \"AT1G01020\" transcript_id \"AT1G01020.1\"\n")
        temp_file.write("1\tTAIR10\texon\t6437\t7069\t.\t-\t.\tgene_id \"AT1G01020\" transcript_id \"AT1G01020.1\"\n")
        temp_file.write("1\tTAIR10\tthree_prime_UTR\t5928\t6263\t.\t-\t.\t"
                        "gene_id \"AT1G01020\" transcript_id \"AT1G01020.1\"\n")
        temp_file.write("1\tTAIR10\texon\t5928\t6263\t.\t-\t.\tgene_id \"AT1G01020\" transcript_id \"AT1G01020.1\"\n")

        temp_file.write("1\tTAIR10\texon\t23146\t24451\t.\t+\t.\tgene_id \"AT1G01040\" transcript_id \"AT1G01040.1\"\n")
        temp_file.write("1\tTAIR10\texon\t24542\t24655\t.\t+\t.\tgene_id \"AT1G01040\" transcript_id \"AT1G01040.1\"\n")
        temp_file.write("1\tTAIR10\texon\t27618\t27713\t.\t+\t.\tgene_id \"AT1G01040\" transcript_id \"AT1G01040.1\"\n")
        temp_file.write("1\tTAIR10\texon\t27803\t28431\t.\t+\t.\tgene_id \"AT1G01040\" transcript_id \"AT1G01040.1\"\n")
        temp_file.write("1\tTAIR10\texon\t27372\t27536\t.\t+\t.\tgene_id \"AT1G01040\" transcript_id \"AT1G01040.1\"\n")
        temp_file.write("1\tTAIR10\texon\t30410\t30816\t.\t+\t.\tgene_id \"AT1G01040\" transcript_id \"AT1G01040.1\"\n")
        temp_file.write("1\tTAIR10\texon\t30902\t31120\t.\t+\t.\tgene_id \"AT1G01040\" transcript_id \"AT1G01040.1\"\n")

        temp_file.write("2\tTAIR10\texon\t11649\t13173\t.\t+\t.\tgene_id \"AT1G01030\" transcript_id \"AT1G01030.1\"\n")
        temp_file.write("2\tTAIR10\texon\t13335\t13714\t.\t+\t.\tgene_id \"AT1G01030\" transcript_id \"AT1G01030.1\"\n")
        temp_file.write("2\tTAIR10\texon\t14649\t15173\t.\t+\t.\tgene_id \"AT1G01030\" transcript_id \"AT1G01030.1\"\n")
        temp_file.write("2\tTAIR10\texon\t16335\t17714\t.\t+\t.\tgene_id \"AT1G01030\" transcript_id \"AT1G01030.1\"\n")

        temp_file.write("2\tTAIR10\texon\t18649\t19222\t.\t+\t.\tgene_id \"AT1G01050\" transcript_id \"AT1G01050.1\"\n")
        temp_file.write("2\tTAIR10\texon\t20145\t20604\t.\t+\t.\tgene_id \"AT1G01050\" transcript_id \"AT1G01050.1\"\n")
        temp_file.write("2\tTAIR10\texon\t21998\t22313\t.\t+\t.\tgene_id \"AT1G01050\" transcript_id \"AT1G01050.1\"\n")
        temp_file.write("2\tTAIR10\texon\t22467\t22811\t.\t+\t.\tgene_id \"AT1G01050\" transcript_id \"AT1G01050.1\"\n")

        temp_file.close()

    def tearDown(self):
        """Tidy up files created for testing."""

        os.remove(self.utrseqont_gtf)
        os.remove(self.seqont_gtf)
        os.remove(self.overlapping_gtf)
        os.remove(self.test_gtf)
        os.remove(self.test_gtf_sq)
        os.remove(self.test_fc_gtf)
        os.remove(self.test_fc)
        os.remove(self.test_fc_t_utr)
        os.removedirs(self.existing_path)

    def test_intron_lengths(self):
        """ Test intron lengths calc. """

        p = GtfParser(self.utrseqont_gtf)
        lengths = p.get_ftr_lengths(p.featureType.intron, usesq=False)[LENGTH].tolist()
        self.assertItemsEqual(lengths, [1499, 1499, 3499, 1097])

        lengths = p.get_ftr_lengths(p.featureType.intron)[LENGTH].tolist()
        self.assertItemsEqual(lengths, [1499, 1097, 1499])

    def test_exon_lengths(self):
        """ Test exon lengths calc. """

        p = GtfParser(self.utrseqont_gtf)
        lengths = p.get_ftr_lengths(p.featureType.exon, usesq=False)[LENGTH].tolist()
        self.assertItemsEqual(lengths, [451, 201, 903, 903, 903, 501, 2001])

        lengths = p.get_ftr_lengths(p.featureType.exon)[LENGTH].tolist()
        self.assertItemsEqual(lengths, [451, 903, 501, 2001])

    def test_exons_by_gene(self):
        """ Test exon counts calc. """

        p = GtfParser(self.utrseqont_gtf)
        exons_by_gene = p.get_exons_per_gene()[EXON_COUNTS].tolist()
        self.assertItemsEqual(exons_by_gene, [11.0/3])

    def test_transcripts_by_gene(self):
        """ Test transcript counts calc. """

        p = GtfParser(self.overlapping_gtf)
        ts_by_gene = p.get_transcripts_per_gene()[TRANSCRIPT_COUNTS].tolist()
        self.assertItemsEqual(ts_by_gene, [3, 2, 3])

    def test_overlapping_genes_exon_lengths(self):
        """ Test exon lengths when genes in +ve and -ve strands overlap """

        p = GtfParser(self.overlapping_gtf)
        lengths = p.get_ftr_lengths(p.featureType.exon, usesq=False)[LENGTH].tolist()
        self.assertItemsEqual(lengths, [451, 201, 903, 903, 903, 501, 2001, 451, 201, 903, 903, 903, 501, 2001,
                                        1411, 3001, 3001])

        lengths = p.get_ftr_lengths(p.featureType.exon)[LENGTH].tolist()
        self.assertItemsEqual(lengths, [451, 903, 501, 2001, 451, 903, 501, 2001,
                                        1411, 3001])

    def test_overlapping_genes_exon_count(self):
        """ Test exon counts when genes in +ve and -ve strands overlap """

        p = GtfParser(self.overlapping_gtf)
        exons_by_gene = p.get_exons_per_gene()[EXON_COUNTS].tolist()
        self.assertItemsEqual(exons_by_gene, [11.0/3, 3.0/2, 11.0/3])

    def test_overlapping_genes_intron_lengths(self):
        """ Test intron lengths when genes in +ve and -ve strands overlap """

        p = GtfParser(self.overlapping_gtf)
        lengths = p.get_ftr_lengths(p.featureType.intron,usesq=False)[LENGTH].tolist()
        self.assertItemsEqual(lengths, [1499, 1499, 3499, 1097, 1499, 1499, 3499, 1097, 1999])

        lengths = p.get_ftr_lengths(p.featureType.intron)[LENGTH].tolist()
        self.assertItemsEqual(lengths, [1499, 1097, 1499, 1499, 1097, 1499, 1999])

    def test_output_to_gtf2_with_sq(self):
        """ Test gff input correctly output to gtf """

        p = GtfParser(self.overlapping_gtf)
        p.export_to_gtf2(os.path.splitext(self.output_gtf)[0])

        # check file exists
        self.assertTrue(os.path.exists(self.output_gtf), "Output gtf file does not exist")

        # check contents correspond to gtf2 format
        new_gtf_handle = open(self.output_gtf, "r")
        test_gtf_handle = open(self.test_gtf_sq, "r")

        newdata = new_gtf_handle.readlines()
        testdata = test_gtf_handle.readlines()

        self.assertListEqual(newdata, testdata)

        new_gtf_handle.close()
        test_gtf_handle.close()

        os.remove(self.output_gtf)

    def test_output_to_gtf2_without_sq(self):
        """ Test gff input correctly output to gtf """

        p = GtfParser(self.overlapping_gtf)
        p.export_to_gtf2(os.path.splitext(self.output_gtf)[0], usesq=False)

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

    def test_first_feature_filter(self):
        """ Test filtering by first exon works correctly. """

        data = pd.read_csv(self.test_fc, header=0, sep="\t")
        p = GtfParser(self.test_fc_gtf)

        filtered = p.filter_by_first_features(data, p.featureType.exon, "Start", "End")
        self.assertEqual(len(filtered.index), 5)
        self.assertListEqual(filtered['mini.bam'].tolist(), [26, 102, 3, 7, 13])

    def test_last_feature_filter(self):
        """ Test filtering by first exon works correctly. """

        data = pd.read_csv(self.test_fc, header=0, sep="\t")
        p = GtfParser(self.test_fc_gtf)
        filtered = p.filter_by_last_features(data, p.featureType.exon, "Start", "End")

        self.assertEqual(len(filtered.index), 5)
        self.assertListEqual(filtered['mini.bam'].tolist(), [12, 4, 8, 57, 77])

    def test_cds_parsing(self):
        """ Test parsing when UTRs are derived from CDSs """

        p = GtfParser(self.seqont_gtf)
        p2 = GtfParser(self.utrseqont_gtf)

        p.export_to_gtf2(os.path.splitext(self.temp_gtf)[0], p.featureType.five_prime_UTR)
        p2.export_to_gtf2(os.path.splitext(self.temp2_gtf)[0], p.featureType.five_prime_UTR)
        data_p = pd.read_csv(self.temp_gtf, header=None, sep="\t")
        data_p2 = pd.read_csv(self.temp2_gtf, header=None, sep="\t")

        data_p.sort_values(by=[3, 4], inplace=True, ascending=True)
        data_p2.sort_values(by=[3, 4], inplace=True, ascending=True)

        self.assertListEqual(data_p[3].tolist(), data_p2[3].tolist())
        self.assertListEqual(data_p[4].tolist(), data_p2[4].tolist())

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

        p = GtfParser("Data/cds.gtf")
        p.export_to_gtf2(os.path.splitext(self.temp_gtf)[0], p.featureType.five_prime_UTR)
        data_p = pd.read_csv(self.temp_gtf, header=None, sep="\t")

        self.assertEquals(len(data_p), 1)   # only one 5' UTR in this data
        self.assertEquals(data_p[3][0], 17494774)
        self.assertEquals(data_p[4][0], 17495033)

        p.export_to_gtf2("temp-p", p.featureType.three_prime_UTR)
        data_p = pd.read_csv("temp-p.gtf", header=None, sep="\t")

        self.assertEquals(len(data_p), 1)   # only one 3' UTR in this data
        self.assertEquals(data_p[3][0], 17493372)
        self.assertEquals(data_p[4][0], 17493624)

        os.remove(self.temp_gtf)
        os.remove(self.temp2_gtf)

        p = GtfParser("Data/cds2.gtf")
        p.export_to_gtf2(os.path.splitext(self.temp_gtf)[0], p.featureType.three_prime_UTR)
        data_p = pd.read_csv(self.temp_gtf, header=None, sep="\t")
        self.assertEquals(len(data_p), 1)

        p = GtfParser("Data/cds3.gtf")
        p.export_to_gtf2(os.path.splitext(self.temp_gtf)[0], [p.featureType.exon, p.featureType.intron,
                                                              p.featureType.three_prime_UTR,
                                                              p.featureType.five_prime_UTR])
        data_p = pd.read_csv(self.temp_gtf, header=None, sep="\t")
        data_p.columns = [CHROMOSOME, SOURCE, TYPE, START, STOP, SCORE, STRAND, PHASE, ATTRIBUTES]
        self.assertEquals(len(data_p[data_p[TYPE] == p.featureType.exon]), 11)
        self.assertEquals(len(data_p[data_p[TYPE] == p.featureType.intron]), 8)
        self.assertEquals(len(data_p[data_p[TYPE] == p.featureType.three_prime_UTR]), 1)
        self.assertEquals(len(data_p[data_p[TYPE] == p.featureType.five_prime_UTR]), 1)

        os.remove(self.temp_gtf)

    def test_build_antisense_excludes_miRNAs(self):
        """ Test transcript counts calc. """

        antisense_gtf = "antisense.gtf"
        p = GtfParser("Data/mini_genes_fixed.gtf")
        p.build_antisense_gtf_gene_only(os.path.splitext(antisense_gtf)[0])

        # check file exists
        self.assertTrue(os.path.exists(antisense_gtf), "Output gtf file does not exist")

        result = pd.read_csv(antisense_gtf, header=None, comment='#', sep='\t', engine='python')
        result.columns = [CHROMOSOME, SOURCE, TYPE, START, STOP, SCORE, STRAND, PHASE, ATTRIBUTES]

        # check we have the correct number of -ve and +ve segments
        self.assertListEqual(result[START].tolist(), [31228, 33379, 23146])
        self.assertListEqual(result[STOP].tolist(), [33153, 37871, 31169])

        os.remove(antisense_gtf)

    def test_utr_parsing(self):
        """ Test parsing when UTRs are derived from CDSs """

        p = GtfParser(self.utrseqont_gtf)

        self.assertEqual(len(p.cds_tr_join), 10)
        self.assertItemsEqual(p.cds_tr_join[START].tolist(), [1201, 3000, 5000, 7000, 1201,
                                                              5000, 7000, 3391, 5000, 7000])
        self.assertItemsEqual(p.cds_tr_join[STOP].tolist(), [1500, 3902, 5500, 7600, 1500,
                                                             5500, 7600, 3902, 5500, 7600])


if __name__ == '__main__':
    unittest.main(verbosity=2)
