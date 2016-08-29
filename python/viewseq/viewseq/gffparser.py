#!/usr/bin/env python

import pandas as pd
import logging
from parser import Parser
from constants import *

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"


class GffParser(Parser):
    def __init__(self, filename, gene_ontology="", exon_ontology="", gene_name="gene", exon_name="exon",
                 transcript_name="mRNA", three_utr_name="three_prime_UTR", five_utr_name="five_prime_UTR",
                 cds_name="CDS"):
        """ Parse gff v3 file into a pandas dataframe.
        :param filename The gff v3 file to parse
        :param gene_name The string used to identify genes in the gff type column
        :param exon_name The string used to identify exons in the gff type column
        :param transcript_name The string used to identify transcripts in the gff type column
        :param three_utr_name The string used to identify 3' UTRs in the gff type column
        :param five_utr_name The string used to identify 5'UTRs in the gff type column
        :param cds_name The string used to identify CDSs in the gff type column
        :raises ValueError if gene, exon or transcript names cannot be found in the file
        :raises IndexError if there are no transcripts in the file
        """

        self.logger = logging.getLogger(__name__)

        # get gene ontology data for alternative type names
        self._gene_ontology = self._read_ontology(gene_ontology, gene_name)
        self._exon_ontology = self._read_ontology(exon_ontology, exon_name)

        self.essential_types = [self._gene_ontology, self._exon_ontology, [transcript_name]]
        super(GffParser, self).__init__(filename, gene_ontology, exon_ontology, gene_name, exon_name, transcript_name,
                                        three_utr_name, five_utr_name, cds_name)

    def _build_region_data(self, gene_name, exon_name, transcript_name, three_utr_name, five_utr_name, cds_name):

        # get genes, exons, transcripts and UTRs
        self.logger.info("Extracting regions for different features")
        self.genes = self._extract_regions(self._gene_ontology, GENE_START, GENE_STOP)
        self.genes = self.genes.rename(columns={ID: GENE_ID})
        self.genes.set_index(GENE_ID, inplace=True)

        self.exons = self._extract_regions(self._exon_ontology, EXON_START, EXON_STOP)
        self.exons.drop([ID, ATTRIBUTES], axis=1, inplace=True)

        self.transcripts = self._extract_regions([transcript_name], RNA_START, RNA_STOP)
        self.transcripts.drop([CHROMOSOME, STRAND, SCORE, PHASE, ATTRIBUTES, SOURCE, NAME, INDEX],
                              axis=1, inplace=True)
        self.transcripts = self.transcripts.rename(columns={PARENT: GENE_ID, ID: RNA_ID})

        three_utrs = self._extract_regions([three_utr_name], THREE_UTR_START, THREE_UTR_STOP)
        five_utrs = self._extract_regions([five_utr_name], FIVE_UTR_START, FIVE_UTR_STOP)
        three_utrs = three_utrs.rename(columns={PARENT: RNA_ID})
        five_utrs = five_utrs.rename(columns={PARENT: RNA_ID})
        three_utrs.drop_duplicates([CHROMOSOME, RNA_ID, STRAND, THREE_UTR_START, THREE_UTR_STOP], inplace=True)
        five_utrs.drop_duplicates([CHROMOSOME, RNA_ID, STRAND, FIVE_UTR_START, FIVE_UTR_STOP], inplace=True)

        # match exons to transcripts
        self.logger.info("Matching exons to transcripts")
        self.ex_tr_join = pd.merge(self.exons[[CHROMOSOME, SOURCE, TYPE, EXON_START, EXON_STOP, SCORE,
                                               STRAND, PHASE, PARENT]],
                                   self.transcripts[[RNA_START, RNA_STOP, RNA_ID, GENE_ID]],
                                   left_on=PARENT, right_on=RNA_ID)

        # non-transcripts are things which are parents of exons, but not mRNA
        non_transcript_ids = (self.exons[(~self.exons[PARENT].isin(self.ex_tr_join[PARENT]))])
        non_transcript_ids = non_transcript_ids.rename(columns={PARENT: RNA_ID})
        non_transcripts = pd.merge(self.df[[TYPE, START, STOP, ID, PARENT]],
                                   non_transcript_ids[[RNA_ID]], left_on=ID, right_on=RNA_ID)
        non_transcripts.drop([ID], axis=1, inplace=True)
        non_transcripts = non_transcripts.rename(columns={PARENT: GENE_ID, START: RNA_START, STOP: RNA_STOP})

        self.transcripts = pd.concat([self.transcripts, non_transcripts])

        # make ex_tr_join again with full set of transcripts
        self.ex_tr_join = pd.merge(self.exons[[CHROMOSOME, SOURCE, TYPE, EXON_START, EXON_STOP, SCORE,
                                               STRAND, PHASE, PARENT]],
                                   self.transcripts[[RNA_START, RNA_STOP, RNA_ID, GENE_ID]],
                                   left_on=PARENT, right_on=RNA_ID)

        self.ex_tr_join = self.ex_tr_join.rename(columns={EXON_START: START, EXON_STOP: STOP})

        self.cds_tr_join = pd.DataFrame()
        derive_cds_from_utrs = three_utrs.empty and five_utrs.empty
        if derive_cds_from_utrs:
            cdss = self._extract_regions([cds_name], CDS_START, CDS_STOP)

            self.cds_tr_join = pd.merge(cdss, self.transcripts[[RNA_ID, GENE_ID, RNA_START,
                                                               RNA_STOP]][self.transcripts[TYPE] == mRNA],
                                        left_on=PARENT, right_on=RNA_ID)
            self.cds_tr_join = self.cds_tr_join.rename(columns={CDS_START: START, CDS_STOP: STOP})
            self.cds_tr_join.drop([GENE_ID], axis=1, inplace=True)

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
        self.introns = self.introns.rename(columns={PARENT: RNA_ID})
        self.in_tr_join = pd.merge(self.introns, self.transcripts[[RNA_ID, GENE_ID]], on=[RNA_ID, GENE_ID])
        # now get rid of duplicates (some duplicates are on different transcripts, so must be done after in_tr_join)
        self.introns.drop_duplicates(subset=[CHROMOSOME, START, STOP, STRAND], inplace=True)

        # match utrs to transcripts
        self.t_utr_tr_join = pd.merge(three_utrs, self.transcripts[[RNA_ID, GENE_ID, RNA_START,
                                                                   RNA_STOP]][self.transcripts[TYPE] == mRNA],
                                      left_on=RNA_ID, right_on=RNA_ID)
        self.t_utr_tr_join = self.t_utr_tr_join.rename(columns={THREE_UTR_START: START, THREE_UTR_STOP: STOP})
        self.f_utr_tr_join = pd.merge(five_utrs, self.transcripts[[RNA_ID, GENE_ID, RNA_START,
                                                                  RNA_STOP]][self.transcripts[TYPE] == mRNA],
                                      left_on=RNA_ID, right_on=RNA_ID)
        self.f_utr_tr_join = self.f_utr_tr_join.rename(columns={FIVE_UTR_START: START, FIVE_UTR_STOP: STOP})

        # need to do this after creating t_utr_tr_join and f_utr_tr_join as process_cds_to_utrs needs the
        if not derive_cds_from_utrs:
            self._process_utrs_to_cds()

    def _preprocess_data(self):
        """ Extract the gff attributes into columns """

        self.df.columns = [CHROMOSOME, SOURCE, TYPE, START, STOP, SCORE, STRAND, PHASE, ATTRIBUTES]

        self.logger.info("Extracting gff attributes to columns")

        # convert any specialised feature names to gene/exon/mRNA/etc
        # self.df = pd.merge(self.df, self.config, left_on=TYPE, right_on='feature')
        # del self.df[TYPE]
        # self.df.rename(columns={'new_feature':TYPE}, inplace=True)

        # extract parent and id attributes, can get the rest later if we need them
        self.df[ID] = self.df[ATTRIBUTES].str.extract('ID=([^;]*)', expand=False)
        self.df[PARENT] = self.df[ATTRIBUTES].str.extract('[Pp]arent=([^;]*)', expand=False)
        self.df[NAME] = self.df[ATTRIBUTES].str.extract('Name=([^;]*)', expand=False)  # could also use Alias

        # split out the parents into columns and then pivot into rows
        s = self.df[PARENT].str.split(',', expand=True).stack()

        # replace parent column with the stacked version
        # (so where there are multiple parents, we get a row for each parent)
        s.index = s.index.droplevel(-1)  # to line up with df's index
        s.name = PARENT  # s needs a name in order to join with df
        del self.df[PARENT]  # delete the current parent column
        self.df = self.df.join(s)  # join the new parent column
        self.df.reset_index(inplace=True)
