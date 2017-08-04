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

import pandas as pd
import logging
from parser import Parser
from constants import *

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"


class GtfParser(Parser):
    def __init__(self, filename, gene_ontology="", exon_ontology="", gene_name="gene", exon_name="exon",
                 transcript_name="mRNA", three_utr_name="three_prime_UTR", five_utr_name="five_prime_UTR",
                 cds_name="CDS"):
        """ Parse gtf v2 file into a pandas dataframe.
        :param filename The gtf v2 file to parse
        :param gene_name The string used to identify genes in the gtf type column
        :param exon_name The string used to identify exons in the gtf type column
        :param transcript_name The string used to identify transcripts in the gtf type column
        :param three_utr_name The string used to identify 3' UTRs in the gtf type column
        :param five_utr_name The string used to identify 5'UTRs in the gtf type column
        :param cds_name The string used to identify CDSs in the gtf type column
        :raises ValueError if exon names cannot be found in the file
        :raises IndexError if there are no transcripts in the file
        """

        self.logger = logging.getLogger(__name__)

        # get gene ontology data for alternative type names
        self._gene_ontology = self._read_ontology(gene_ontology, gene_name)
        self._exon_ontology = self._read_ontology(exon_ontology, exon_name)
        ess_types = [x for x in self._exon_ontology]
        ess_types.append(cds_name)

        self.essential_types = [ess_types]
        super(GtfParser, self).__init__(filename, gene_ontology, exon_ontology, gene_name, exon_name, transcript_name,
                                        three_utr_name, five_utr_name, cds_name)

    def _preprocess_data(self):
        """ Extract the gtf attributes into columns """

        self.df.columns = [CHROMOSOME, SOURCE, TYPE, START, STOP, SCORE, STRAND, PHASE, ATTRIBUTES]

        self.logger.info("Extracting gtf attributes to columns")

        # extract transcript_id and gene_id attributes
        self.df[GENE_ID] = self.df[ATTRIBUTES].str.extract('[Gg]ene_id \"([^\"]*)', expand=False)
        self.df[RNA_ID] = self.df[ATTRIBUTES].str.extract('[Tt]ranscript_id \"([^\"]*)', expand=False)
        self.df[NAME] = self.df[ATTRIBUTES].str.extract('[Gg]ene_name \"([^\"]*)', expand=False)

    def _build_region_data(self, gene_name, exon_name, transcript_name, three_utr_name, five_utr_name, cds_name):

        # get exons and UTRs
        self.logger.info("Extracting regions for different features")

        self.exons = self._extract_regions(self._exon_ontology, EXON_START, EXON_STOP)
        self.exons.drop([ATTRIBUTES], axis=1, inplace=True)

        self.genes = self.exons.groupby(GENE_ID).first()
        self.genes.drop([RNA_ID, TYPE, EXON_START, EXON_STOP, SCORE, STRAND, PHASE], axis=1, inplace=True)

        # make the transcripts
        firstexons = self._get_ftrs_by_position(self.exons, +1)
        lastexons = self._get_ftrs_by_position(self.exons, -1)
        self.transcripts = pd.merge(firstexons[[CHROMOSOME, TYPE, RNA_ID, GENE_ID, STRAND, EXON_START]],
                                    lastexons[[CHROMOSOME, TYPE, RNA_ID, GENE_ID, STRAND, EXON_STOP]],
                                    on=[CHROMOSOME, RNA_ID, GENE_ID, STRAND])
        self.transcripts = self.transcripts.rename(columns={EXON_START: RNA_START, EXON_STOP: RNA_STOP})
        self.transcripts[TYPE] = mRNA

        three_utrs = self._extract_regions([three_utr_name], THREE_UTR_START, THREE_UTR_STOP)
        five_utrs = self._extract_regions([five_utr_name], FIVE_UTR_START, FIVE_UTR_STOP)

        three_utrs.drop([GENE_ID], axis=1, inplace=True)
        five_utrs.drop([GENE_ID], axis=1, inplace=True)
        three_utrs.drop_duplicates([CHROMOSOME, RNA_ID, STRAND, THREE_UTR_START, THREE_UTR_STOP], inplace=True)
        five_utrs.drop_duplicates([CHROMOSOME, RNA_ID, STRAND, FIVE_UTR_START, FIVE_UTR_STOP], inplace=True)

        # match exons to transcripts
        self.ex_tr_join = self.exons
        self.ex_tr_join = self.ex_tr_join.rename(columns={EXON_START: START, EXON_STOP: STOP})

        self.cds_tr_join = pd.DataFrame()
        derive_cds_from_utrs = three_utrs.empty and five_utrs.empty
        if derive_cds_from_utrs:
            cdss = self._extract_regions([cds_name], CDS_START, CDS_STOP)

            self.cds_tr_join = pd.merge(cdss, self.transcripts[[RNA_ID, RNA_START, RNA_STOP]],
                                        left_on=[RNA_ID], right_on=[RNA_ID])
            self.cds_tr_join = self.cds_tr_join.rename(columns={CDS_START: START, CDS_STOP: STOP})
            self.cds_tr_join[ID] = ""

            three_utrs, five_utrs = self._process_cds_to_utrs()

        # initialise unique exons (up to UTR differences)
        self.unique_exons = self.ex_tr_join

        # update unique exons to include UTR differences
        self._build_unique_exons(self.ex_tr_join, three_utrs, five_utrs)
        self.unique_exons.drop_duplicates(
            [START, STOP, GENE_ID, THREE_UTR_START, THREE_UTR_STOP, FIVE_UTR_START, FIVE_UTR_STOP],
            inplace=True)

        # get unique introns
        self.introns = self._build_introns(self.ex_tr_join)
        self.in_tr_join = self.introns
        self.in_tr_join = pd.merge(self.in_tr_join, self.ex_tr_join[[RNA_ID]],
                                   left_on=[PARENT], right_on=[RNA_ID])
        self.in_tr_join.drop([PARENT], axis=1, inplace=True)
        self.in_tr_join.drop_duplicates(subset=[CHROMOSOME, RNA_ID, START, STOP, STRAND], inplace=True)

        # now get rid of duplicates (some duplicates are on different transcripts, so must be done after in_tr_join)
        self.introns.drop_duplicates(subset=[CHROMOSOME, START, STOP, STRAND], inplace=True)

        # match utrs to transcripts
        self.t_utr_tr_join = three_utrs
        self.t_utr_tr_join = self.t_utr_tr_join.rename(columns={THREE_UTR_START: START, THREE_UTR_STOP: STOP})
        self.t_utr_tr_join = pd.merge(self.t_utr_tr_join, self.ex_tr_join[[GENE_ID, RNA_ID]],
                                      left_on=[RNA_ID], right_on=[RNA_ID])
        self.t_utr_tr_join.drop_duplicates(subset=[CHROMOSOME, RNA_ID, START, STOP, STRAND], inplace=True)
        self.f_utr_tr_join = five_utrs
        self.f_utr_tr_join = self.f_utr_tr_join.rename(columns={FIVE_UTR_START: START, FIVE_UTR_STOP: STOP})
        self.f_utr_tr_join = pd.merge(self.f_utr_tr_join, self.ex_tr_join[[GENE_ID, RNA_ID]],
                                      left_on=[RNA_ID], right_on=[RNA_ID])
        self.f_utr_tr_join.drop_duplicates(subset=[CHROMOSOME, RNA_ID, START, STOP, STRAND], inplace=True)

        # tag non-mRNA entries in transcripts
        # Need .ix and to list TYPE after condition in order to get in place update
        self.transcripts.ix[(~self.transcripts[RNA_ID].isin(three_utrs[RNA_ID])), TYPE] = "not{}".format(mRNA)

        # need to do this after creating t_utr_tr_join and f_utr_tr_join as process_cds_to_utrs needs them
        if not derive_cds_from_utrs:
            self._process_utrs_to_cds()
