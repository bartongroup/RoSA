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

import logging
import warnings
import pandas as pd
import numpy as np
import csv
from collections import namedtuple
from abc import ABCMeta, abstractmethod
from scipy.sparse.csgraph import connected_components
from six import iteritems  # for python 2/3 compatibility

from constants import *

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"

FeatureType = namedtuple('FeatureType', 'exon intron three_prime_UTR five_prime_UTR CDS antisense')
FeatureLookup = {'exon': 'Exon', 'intron': 'Intron', 'three_prime_UTR': "3' UTR", 'five_prime_UTR': "5' UTR",
                 'CDS': 'CDS', 'antisense': 'antisense'}


class Parser(object):
    """ Base class for annotation file parsers. Parses the data into a pandas dataframe
        and then provides access to a variety of summary statistics.
    """

    __metaclass__ = ABCMeta

    def __init__(self, filename, gene_ontology="", exon_ontology="", gene_name="gene", exon_name="exon",
                 transcript_name="mRNA", three_utr_name="three_prime_UTR", five_utr_name="five_prime_UTR",
                 cds_name="CDS"):
        """ Parses input file into a pandas dataframe.
        :param filename The file to parse
        :param gene_ontology File containing details of gene ontology,if missing only looks for gene_name labels
        :param exon_ontology File containing details of exon ontology, if missing only looks for exon_name labels
        :param gene_name The string used to identify genes in the type column
        :param exon_name The string used to identify exons in the type column
        :param transcript_name The string used to identify transcripts in the type column
        :param three_utr_name The string used to identify 3' UTRs in the type column
        :param five_utr_name The string used to identify 5'UTRs in the type column
        :param cds_name The string used to identify CDSs in the type column
        :raises ValueError if essential region names cannot be found in the file
        :raises IndexError if there are no transcripts in the file
        """

        self.featureType = FeatureType(exon=exon_name, intron="intron", three_prime_UTR=three_utr_name,
                                       five_prime_UTR=five_utr_name, CDS=cds_name, antisense="antisense")

        self.logger = logging.getLogger(__name__)

        # read in the gff file
        self.logger.info("Reading in annotation file {}".format(filename))

        # separately count '#' comment lines from top of file, to workaround pandas not handling
        # '#' as comment and part of names...
        i = 0   # number of lines to skip
        with open(filename) as fp:
            for i, line in enumerate(fp):
                if not line.startswith('#'):
                    break

        self.df = pd.read_csv(filename, skiprows=i, header=None, sep='\t', engine='python', quotechar='"')
        # Note - don't use comment='#', as some genes might have # in the name (yes, seriously)

        self._preprocess_data()

        # force CHROMOSOME to be string, in case it wasn't (need for later processing)
        self.df[CHROMOSOME] = self.df[CHROMOSOME].astype(str)

        # get the types in this file
        self.alltypes = self.df[TYPE].value_counts().index
        self._check_types(self.essential_types)
        try:
            self._check_types([[three_utr_name, five_utr_name, cds_name]])
        except ValueError as e:
            # if CDS, 3' or 5' UTR names are missing, just give a warning
            warnings.warn(str(e), UserWarning)

        self._build_region_data(gene_name, exon_name, transcript_name, three_utr_name, five_utr_name, cds_name)

        self._make_squashed_extents()

        self.featuredict = {self.featureType.exon: self.ex_tr_join,
                            self.featureType.intron: self.in_tr_join,
                            self.featureType.three_prime_UTR: self.t_utr_tr_join,
                            self.featureType.five_prime_UTR: self.f_utr_tr_join,
                            self.featureType.CDS: self.cds_tr_join}

        self.sqdict = {self.featureType.exon: self.sq_ex,
                       self.featureType.intron: self.sq_in,
                       self.featureType.three_prime_UTR: self.sq_tutr,
                       self.featureType.five_prime_UTR: self.sq_futr}

        # uses self.featuredict in creation
        self.featuredict[self.featureType.antisense] = self.build_antisense_gtf_gene_only()

        # if we used squashed dict then the antisense features will be squashed too
        self.sqdict[self.featureType.antisense] = self.featuredict[self.featureType.antisense]

    @abstractmethod
    def _preprocess_data(self):
        """ Convert the data attributes column into separate columns for parent/gene_id/rna_id
        """
        pass

    @abstractmethod
    def _build_region_data(self, gene_name, exon_name, transcript_name, three_utr_name, five_utr_name, cds_name):
        """ Build data structures for regions in the data
        :param gene_name: name used for genes
        :param exon_name: name used or exons
        :param transcript_name: name used for transcripts
        :param three_utr_name: name used for 3' UTRs
        :param five_utr_name: name used for 5' UTRs
        :param cds_name: name used for CDSs
        """
        pass

    def _make_squashed_extents(self):
        """ Build merged extents of features for every region. """

        self.sq_ex = self._squash_region_data(self.ex_tr_join)
        self.sq_tutr = self._squash_region_data(self.t_utr_tr_join)
        self.sq_futr = self._squash_region_data(self.f_utr_tr_join)

        self.sq_in = self._build_introns(self.sq_ex)

        # add a fake RNAID in to keep later processing happy - some grouping is done by RNAID
        self.sq_in[RNA_ID] = self.sq_in[GENE_ID]

    def _squash_region_data(self, extents):
        """ Squash extents down on a per gene basis so that we have the maximum covered extent
            of the extents in the gene e.g.
            ==============--------------==========
                  ============------====================
            becomes
            ==================------====================
        """

        # don't try to do anything if extents is empty, just return an empty df
        if extents.empty:
            return pd.DataFrame

        #to try to avoid memory issues, break down into chromosomes
        chrs = extents[CHROMOSOME].value_counts().index.tolist()
        strands = ["+", "-"]

        full_sq_set = pd.DataFrame()
        squashset = pd.DataFrame()

        for chro in chrs:
            for strand in strands:

                current_extents = extents.loc[extents[CHROMOSOME] == chro].loc[extents[STRAND]==strand]
                squashset_len = 0 # try to avoid issues measuring length of empty dataframe
                old_len = len(current_extents.index)

                # repeatedly loop trying to squash overlapping exons
                # may not be the fastest approach...
                while old_len != squashset_len:

                    # for extents do an inner join on gene id
                    squashset = pd.merge(current_extents[[GENE_ID, START, STOP, SOURCE, TYPE, SCORE, PHASE]],
                                         current_extents[[GENE_ID, START, STOP]], on=GENE_ID)

                    # select overlapping things
                    squashset = squashset[((squashset["Start_y"] <= squashset["Start_x"]) & (squashset["Stop_y"] >= squashset["Start_x"])) |
                                          ((squashset["Stop_y"] >= squashset["Stop_x"]) & (squashset["Start_y"] <= squashset["Stop_x"]))]

                    # calculate min start and max end for each row
                    squashset["minstart"] = squashset.loc[:, ["Start_x", "Start_y"]].min(axis=1)
                    squashset["maxend"] = squashset.loc[:, ["Stop_x", "Stop_y"]].max(axis=1)

                    # calculate length between min start and max end
                    squashset["sq_length"] = squashset["maxend"] - squashset["minstart"]

                    # pull out row for each exon which has max length
                    indexes = squashset.groupby([GENE_ID, "Start_x", "Stop_x"])["sq_length"].idxmax()

                    squashset = squashset.loc[indexes.reset_index()["sq_length"]]
                    squashset.drop_duplicates([GENE_ID, "minstart", "maxend"], inplace=True)
                    squashset[CHROMOSOME] = chro
                    squashset[STRAND] = strand

                    old_len = len(current_extents.index)
                    squashset_len = len(squashset.index)
                    current_extents = squashset[[GENE_ID, CHROMOSOME, STRAND, SOURCE, TYPE, SCORE, PHASE, "minstart", "maxend"]]
                    current_extents.columns = [GENE_ID, CHROMOSOME, STRAND, SOURCE, TYPE, SCORE, PHASE, START, STOP]


                full_sq_set = pd.concat([full_sq_set, squashset], axis=0)

        squashset = full_sq_set[[GENE_ID, CHROMOSOME, STRAND, SOURCE, TYPE, SCORE, PHASE, "minstart", "maxend"]].copy()
        squashset.columns = [GENE_ID, CHROMOSOME, STRAND, SOURCE, TYPE, SCORE, PHASE, START, STOP]
        squashset.drop_duplicates([GENE_ID, START, STOP], inplace=True)
        squashset[RNA_ID] = squashset[GENE_ID]  # set to geneid so that intron calcs can differentiate

        return squashset

    def _read_ontology(self, ontology, basename):
        """ Read in data from ontology files
        :param ontology: File containing ontology
        :param basename: Name of base type in case we have no ontology file
        :return ontology list
        """

        if ontology != "":

            self.logger.info("Reading in data from ontology file: {}".format(ontology))

            df = pd.read_csv(ontology, header=0, sep="\t", engine='python')

            self_ontology = pd.concat([df["Name"], df["Accession"]])
            self_ontology = self_ontology.reset_index(drop=True).tolist()
        else:

            self.logger.debug("Using base ontology")

            self_ontology = [basename]

        return self_ontology

    def get_data_by_feature(self):
        """ Return list of tuples of (ftr name, dataframe) where dataframe corresponds to the ftr
            :return: list of tuples of (ftr name, dataframe)
        """

        result_list = []
        self.logger.debug("Calculating data by feature type")
        for name, feature_df in iteritems(self.featuredict):

            if name == self.featureType.exon:
                # treat exons differently
                # include exons which differ by UTR
                feature_df = self.unique_exons

                # give ftr a length attribute
                feature_df[LENGTH] = feature_df[STOP] - feature_df[START] + 1
                feature_df.drop([SCORE, SOURCE, STRAND, TYPE, PHASE], axis=1, inplace=True)

                result_list.append((name, feature_df[[GENE_ID, LENGTH]]))

            elif (name == self.featureType.five_prime_UTR) or (name == self.featureType.three_prime_UTR):

                # give ftr a length attribute
                feature_df[LENGTH] = feature_df[STOP] - feature_df[START] + 1

                # remove duplicates at transcript level
                feature_df.drop_duplicates([CHROMOSOME, GENE_ID, STRAND, START, STOP], inplace=True)
                feature_df.drop([SCORE, SOURCE, STRAND, TYPE, PHASE], axis=1, inplace=True)

                # UTRs should have lengths summed across transcripts
                feature_df = feature_df.groupby(RNA_ID).aggregate({LENGTH: 'sum',
                                                                   GENE_ID: 'first'}).reset_index()

                result_list.append((name, feature_df[[GENE_ID, LENGTH]]))

            elif not feature_df.empty:

                # give ftr a length attribute
                feature_df[LENGTH] = feature_df[STOP] - feature_df[START] + 1

                # remove duplicates at transcript level
                feature_df.drop_duplicates([CHROMOSOME, GENE_ID, STRAND, START, STOP], inplace=True)
                feature_df.drop([SCORE, SOURCE, STRAND, TYPE, PHASE], axis=1, inplace=True)

                result_list.append((name, feature_df[[GENE_ID, LENGTH]]))

        return result_list

# seem unused
   # def get_ftrdict(self, name):
   #     return self.featuredict[name]

    def get_ftr_lengths(self, ftr, bygene=False, usesq=True):
        """ Return lengths of features of type ftr, individually or by gene. Uses squashed representation.
        :param ftr: the type of feature to get feature lengths for
        :param bygene: True if lengths to be summed across the gene
        :return: pandas Series of feature lengths
        """

        self._check_feature_type(ftr)
        self.logger.debug("Calculating feature lengths: feature is {}".format(ftr))
        if usesq:
            return self._get_lengths(self.sqdict[ftr], bygene)
        else:
            if ftr == self.featureType.exon:
                return self._get_lengths(self.unique_exons, bygene)
            elif ftr == self.featureType.intron:
                return self._get_lengths(self.introns, bygene)
            else:
                self.logger.debug("Calculating feature lengths: feature is {}".format(ftr))
                return self._get_lengths(self.featuredict[ftr], bygene)

    def _get_lengths(self, ftr_df, bygene):

        self.logger.debug("Getting lengths")

        lengths = ftr_df[STOP] - ftr_df[START] + 1
        lengths.name = LENGTH

        result = pd.concat([ftr_df[GENE_ID], lengths], axis=1)

        if not bygene:
            return result
        else:
            return result.groupby(GENE_ID).sum()

    def get_exons_per_gene(self):
        """ Get number of exons per gene
        :return: average number of exons per transcript, per gene,
        or an empty dataframe if no gene can be identified for any exon
        (can happen if gene IDs are missing)
        """
        counts = self._get_ftrcount_per_gene("exon", self.ex_tr_join)
        transcripts = self._get_ftrcount_per_gene("transcript", self.transcripts)

        # translate exon count to average number of exons per transcript per gene
        counts = pd.merge(counts, transcripts, left_index=True, right_index=True)
        counts["num_exons"] = counts["exon_counts"] / counts["transcript_counts"]
        return pd.DataFrame(counts["num_exons"])

    def get_transcripts_per_gene(self):
        """ Get number of transcripts per gene
        :return: number of unique transcripts, per gene, or an empty dataframe if no gene can be identified
        (can happen if gene IDs are missing)
        """
        return self._get_ftrcount_per_gene("transcript", self.transcripts)

    def get_gene_names(self):
        """ Get names of genes
        :return: gene names by gene id
        """
        return self.genes[[NAME]]

    def _get_ftrcount_per_gene(self, name, frame):
        """
        Count number of given feature per gene
        :param name: the name of the ftr to count per gene
        :return: number of a ftr per gene
        """

        self.logger.debug("Calculating {} per gene".format(name))
        counts_name = "{}_counts".format(name.lower())
        try:
            counts = frame.groupby(GENE_ID).size()
            counts.name = counts_name
            return counts.to_frame()
        except ValueError:
            # all geneIDs are NaNs, no exon counts per gene possible
            df = pd.DataFrame(columns=[GENE_ID, counts_name])
            return df

    def get_metadata(self):
        """
        Get totals of genes, transcripts etc into a single dataframe
        :return: dataframe containing totals
        """

        l = {}

        for name, feature_df in iteritems(self.featuredict):

            if name == self.featureType.exon:
                # treat exons differently
                # include exons which differ by UTR
                feature_df = self.unique_exons

            elif name == self.featureType.antisense:
                # don't output antisense 'ftrs'
                continue

            l["Total number of {}s".format(FeatureLookup[name])] = "{:,}".format(len(feature_df.index))

        # do transcripts as well
        l["Total number of transcripts"] = "{:,}".format(len(self.transcripts[self.transcripts[TYPE] == mRNA].index))

        df = pd.DataFrame(data=list(l.items()), columns=["name", "value"])
        return df

    def filter_by_first_features(self, data, feature, start_name, stop_name, usesq=True):
        """ Filter data so that it only contains data for the first feature in each transcript.

        :param data: dataframe to be filtered by first features
        :param feature: the feature type to filter on (e.g. exon, intron). Must be of type parser.FeatureType.
        :param start_name: the name of the column containing feature start positions
        :param stop_name: the name of the column containing feature stop positions
        :return: data limited to first features
        :raise ValueError if feature is not a parser.FeatureType
        """
        self._check_feature_type(feature)
        return self._filter_by_feature_pos(data, feature, start_name, stop_name, 1, usesq)

    def filter_by_last_features(self, data, feature, start_name, stop_name, usesq=True):
        """ Filter data so that it only contains data for the last feature in each transcript.

        :param data: dataframe to be filtered by last features
        :param feature: the feature type to filter on (e.g. exon, intron)
        :param start_name: the name of the column containing feature start positions
        :param stop_name: the name of the column containing feature stop positions
        :return: data limited to last features
        :raise ValueError if feature is not a parser.FeatureType
        """
        self._check_feature_type(feature)
        return self._filter_by_feature_pos(data, feature, start_name, stop_name, -1, usesq)

    def _filter_by_feature_pos(self, data, feature, start_name, stop_name, pos, usesq):
        """ Filter data so that only features in first or last position (set by pos) in each transcript are returned.
        :param data: dataframe to be filtered by first/last features
        :param feature: the feature type to filter on (e.g. exon, intron)
        :param start_name: the name of the column containing feature start positions
        :param stop_name: the name of the column containing feature stop positions
        :param pos: 1 for first feature, -1 for last feature
        :param usesq: True if using squashed features, false if using full dataset
        :return: data limited to first/last features
        :raise ValueError if feature is not a parser.FeatureType
        """
        self._check_feature_type(feature)

        self.logger.debug("Filtering {} data by {position} entry in gene".
                         format(feature, position="first" if pos == 1 else "last"))

        if pos not in [-1, 1]:
            raise ValueError("Feature position must be -1 or 1")

        # rename start and stop columns in data to match ours
        data = data.rename(columns={start_name: START, stop_name: STOP})

        if usesq:
            join = self.sqdict[feature]
        else:
            join = self.featuredict[feature]
        # NB both start and stop needed in case of features in different transcripts starting or
        # ending at same position - we consider features with the same start AND stop points to be the same feature.
        full_filter = pd.merge(join[[START, STOP, RNA_ID]], data, on=[START, STOP])

        return self._get_ftrs_by_position(full_filter, pos)

    def _get_ftrs_by_position(self, full_filter, pos):
        """ Filter features from dataframe full_filter by first/last position, depending on value of pos.
        :param full_filter dataframe to filter
        :param pos: first (1) or last (-1) feature
        :return filtered dataframe
        """

        self.logger.debug("Filtering features by position")

        # default to first feature
        read_first_from_strand = '+'
        read_last_from_strand = '-'
        if pos == -1:  # last feature
            read_first_from_strand = '-'
            read_last_from_strand = '+'

        pos_filter = full_filter.loc[full_filter[STRAND] == read_first_from_strand]
        neg_filter = full_filter.loc[full_filter[STRAND] == read_last_from_strand]
        pos_filter = pos_filter.groupby(RNA_ID).first()
        neg_filter = neg_filter.groupby(RNA_ID).last()
        full_filter = pd.concat([pos_filter, neg_filter])

        full_filter = full_filter.reset_index()
        return full_filter

    def export_to_gtf2(self, filepath, exportlist=None, feature_name=None, split=False, out_to_file=True, usesq=True):
        """ Export the annotation to gtf2 (e.g. for using with featureCounts).
        :param filepath: path to export file(s) to, without gtf extension
        :param exportlist: list of features to be exported
        :param feature_name: name to apply to *all* features exported (unset, all features retain their existing name)
        :param split: True if to split gtf files out by chromosome and strand; default False
        :param out_to_file: True if gtf2 data should be output to file (otherwise returned as dataframe)
        :param usesq: True if squashed dataset to be used
        :return list of gtf files created
        """

        # Assume that:
        # - gtf2 has the same column format as gff3
        # - only need to output exportlist defined rows, with extra attribution
        # - extra attribution is gene_id (parent of feature) and transcript_id (joined mRNA id)

        if exportlist is None:
            exportlist = ['exon', 'intron']

        self.logger.debug("Exporting to gtf file: {}.gtf".format(filepath))

        # set up output columns
        alldata = []

        # select correct feature dict to use
        if usesq:
            dict = self.sqdict
        else:
            dict = self.featuredict

        for name, feature in iteritems(dict):
            if name in exportlist:
                try:
                    data = feature[[CHROMOSOME, SOURCE, TYPE, START, STOP, SCORE, STRAND, PHASE, GENE_ID, RNA_ID]]
                except KeyError:
                    # KeyError can be if RNA_ID is missing because we are using squashed representation, so try again
                    feature[RNA_ID] = ""
                    data = feature[[CHROMOSOME, SOURCE, TYPE, START, STOP, SCORE, STRAND, PHASE, GENE_ID, RNA_ID]]
                alldata.append(data)

        if len(exportlist) != len(alldata):
            self.logger.warn("When exporting gtfs at least one entry in the export list was not found. "
                             "Export list = [{}]".format(", ".join(exportlist)))

        try:
            features = pd.concat(alldata)
        except ValueError as e:
            self.logger.error("Error while concatenating gtf output data\n")
            self.logger.error(e.message)
            self.logger.error("Selected output features are: {}".format(exportlist))
            raise

        if feature_name is not None:
            features[TYPE] = feature_name

        # set START and STOP column to int - if the columns had/have NaNs they will have been set to float
        # which produces undesirable gtf output
        #features[[START, STOP]] = features[[START, STOP]].astype(int)
        self._fix_gtf2_col_types(features)

        if out_to_file:
            files = []
            if not split:
                f = self._output_gtf2(features, filepath)
                files.append(f)
            else:
                chrs = self.get_chromosomes()
                for c in chrs:
                    for strand in ["+", "-"]:
                        strand_slice = features[(features[CHROMOSOME] == c) & (features[STRAND] == strand)]
                        slicepath = "{}{}{}".format(filepath, c, strand)
                        f = self._output_gtf2(strand_slice, slicepath)
                        files.append(f)

            return files

        else:
            return features

    def _fix_gtf2_col_types(self, df):
        """ Pre-empt any output format issues by changing type of columns """

        # Would be better in _output_gtf2 but working with slices in export_to_gtf2
        # means we get a SettingWithCopy warning if not done in advance

        self.logger.debug("Setting gtf2 column types to fixed types")

        # convert start/stop to int - can be set to float if NaNs present
        df[[START, STOP]] = df[[START, STOP]].astype(int)

        # build new attribute column for gtf2
        df[GTF_ATTRIBUTES] = 'gene_id "' + df[GENE_ID].astype(str) + '"; ' + \
                             'transcript_id "' + df[RNA_ID].astype(str) + '";'

    def _output_gtf2(self, features, filepath):
        """ Actual csv output code for output to gtf2
        :param features: features to output
        :param filepath: basename of file to output to
        :return: full output filepath
        """

        filepath = "{}.gtf".format(filepath)
        self.logger.debug("Outputting gtf2 to file {}".format(filepath))

        # Output to gtf
        features.to_csv(path_or_buf=filepath,
                        sep='\t',
                        header=False,
                        index=False,
                        quoting=csv.QUOTE_NONE,
                        columns=[CHROMOSOME, SOURCE, TYPE, START, STOP, SCORE, STRAND, PHASE, GTF_ATTRIBUTES])

        return filepath

    def get_chromosomes(self):
        """ Get a list of names in the chromosome column of the annotation.
        :return: list of names of chromosomes in the annotation
        """

        chromosomes = self.df[CHROMOSOME].value_counts().index.tolist()
        chromosomes = [str(i) for i in chromosomes]  # make sure they are strings
        return chromosomes

    def _check_types(self, type_names):
        """ Check that each name in type_names is present in the values of the types column
        :param type_names: list of lists of type_names we need: at least one entry in each list is required
        :raise: ValueError if any type name is missing
        """
        self.logger.debug("Checking types exist in gff file: {}".format(
            ", ".join([item for sublist in type_names for item in sublist])))

        for typelist in type_names:
            setdiff = set(typelist) - set(self.alltypes)
            if len(setdiff) == len(typelist):
                msg = "Some expected type names are missing from the annotation input file " +\
                      "(at least one of these should be present): {}."\
                      .format(" ".join(setdiff))
                raise ValueError(msg)

    def _extract_regions(self, region_types, startname, stopname):
        """ Extract all regions corresponding to region_type in the data
        :param region_types: list of gff types to extract
        :param startname: new name for start position
        :param stopname: new name for stop position
        :return: dataframe for this type
        """
        regions = self.df.loc[self.df[TYPE].isin(region_types)]
        regions = regions.rename(columns={
            START: startname,
            STOP: stopname
        })
        return regions

    def _build_unique_exons(self, ex_tr_join, three_utrs, five_utrs):
        """ Identify and separate out exons with different UTRs. """

        # Identify unique exons at the level of different UTRs by joining the 3' and 5' UTRs to the existing
        # exon-transcript join. Then join the result to the existing list of unique exons to get an exhaustive list.
        #
        # For joining the UTRs to the exon-transcript join:
        # On the +ve strand 5' UTRs should have same start as their corresponding exon;
        # on the -ve strand 5' UTRs should have the same stop as their corresponding exon.
        # On the +ve strand 3' UTRs should have same stop as their corresponding exon;
        # on the -ve strand 3' UTRs should have the same start as their corresponding exon.
        #
        # Account for directionality by joining negative strands separately and then concatenating.

        self.logger.debug("Identifying unique exons")

        # identify UTR-unique exons
        ex_3utr_join = pd.merge(ex_tr_join[[CHROMOSOME, START, STOP, RNA_ID]], three_utrs,
                                left_on=[RNA_ID, STOP], right_on=[RNA_ID, THREE_UTR_STOP])
        ex_3utr_join_neg = pd.merge(ex_tr_join[[CHROMOSOME, START, STOP, RNA_ID]], three_utrs,
                                    left_on=[RNA_ID, START], right_on=[RNA_ID, THREE_UTR_START])
        ex_3utr_join = pd.concat([ex_3utr_join, ex_3utr_join_neg])
        ex_5utr_join = pd.merge(ex_tr_join[[CHROMOSOME, START, STOP, RNA_ID]], five_utrs,
                                left_on=[RNA_ID, START], right_on=[RNA_ID, FIVE_UTR_START])
        ex_5utr_join_neg = pd.merge(ex_tr_join[[CHROMOSOME, START, STOP, RNA_ID]], five_utrs,
                                    left_on=[RNA_ID, STOP], right_on=[RNA_ID, FIVE_UTR_STOP])
        ex_5utr_join = pd.concat([ex_5utr_join, ex_5utr_join_neg])

        # join UTR-unique exons to the unique exons list
        self.unique_exons = pd.merge(self.unique_exons[[CHROMOSOME, START, STOP, RNA_ID, GENE_ID]], ex_5utr_join,
                                     on=[START, STOP, RNA_ID], how='left')
        self.unique_exons = pd.merge(
            self.unique_exons[[CHROMOSOME, START, STOP, RNA_ID, GENE_ID, FIVE_UTR_START, FIVE_UTR_STOP]],
            ex_3utr_join, on=[START, STOP, RNA_ID], how='left')

    def _build_introns(self, ex_tr_join):
        """ Create table of introns from exon-transcript join ex_tr_join. """

        # Idea is to build a next_exon_start column which holds the start position (in bp) of the exon after the
        # exon on the current line. Then the introns can be calculated as the region between exon_stop and
        # next_exon_start.
        #
        # next_exon_start is calculated by first sorting the exon start points, accounting for strand, chromosome
        # and transcripts, to avoid exons from different transcripts becoming muddled up. Then next_exon_start is
        # populated with exon_start and the whole column shifted up 1 position. We do the same for RNA_ID because
        # we can then identify when the next exon no longer belongs to the same transcript, by comparing RNA ids.
        # When this happens next_exon_start is set to NaN which eliminates it from processing as an intron end point.
        #
        # Introns are now produced by removing rows with NaNs, and duplicates, then adding 1 to the previous exon
        # stop position, and subtracting 1 from the next exon start position.

        self.logger.debug("Building introns")

        ex_tr_join.sort_values(by=[STRAND, CHROMOSOME, RNA_ID, START], inplace=True, ascending=True)

        # copy RNA_ID and exon_start into next_RNA and next_exon_start and shift column up 1
        ex_tr_join[[NEXT_RNA, NEXT_EXON_START]] = ex_tr_join[[RNA_ID, START]].shift(-1)

        # set next_exon_start to NaN when the RNA_ID <> next_RNA
        ex_tr_join[[NEXT_EXON_START]] = \
            ex_tr_join[[NEXT_EXON_START]].where(ex_tr_join[RNA_ID] == ex_tr_join[NEXT_RNA], np.nan)

        # create intron columns
        intron_start = ex_tr_join[STOP]
        intron_stop = ex_tr_join[NEXT_EXON_START]

        # combine introns with details from ex_tr_join
        introns = pd.concat([ex_tr_join[CHROMOSOME], ex_tr_join[RNA_ID], ex_tr_join[SOURCE],
                             ex_tr_join[SCORE], ex_tr_join[STRAND], ex_tr_join[PHASE], ex_tr_join[GENE_ID],
                             intron_start, intron_stop],
                            keys=[CHROMOSOME, PARENT, SOURCE, SCORE, STRAND, PHASE, GENE_ID, START, STOP], axis=1)
        introns[TYPE] = self.featureType.intron

        introns.dropna(inplace=True)  # get rid of any rows with NaNs
        introns[START] += 1  # intron starts 1bp after exon stop
        introns[STOP] -= 1  # intron ends 1bp before exon start
        introns[STOP] = introns[STOP].astype(int)  # without this sometimes these are output to gtf as float

        return introns

    def _check_feature_type(self, feature):
        """
        :param feature: a feature to check
        :raise ValueError is the feature is not a FeatureType
        """

        if feature not in self.featureType._fields:
            raise ValueError("Feature type {} is not valid".format(feature))

    def _process_utrs_to_cds(self):
        """ Convert UTR data to CDS data. """

        # first aggregate the 3' and 5' UTRs for each transcript into a single extent per transcript
        tutrs = self.t_utr_tr_join.groupby(RNA_ID).aggregate({START: 'min',
                                                              STOP: 'max',
                                                              CHROMOSOME: 'first',
                                                              STRAND: 'first'}).reset_index()

        futrs = self.f_utr_tr_join.groupby(RNA_ID).aggregate({START: 'min',
                                                              STOP: 'max',
                                                              CHROMOSOME: 'first',
                                                              STRAND: 'first'}).reset_index()

        # join 3'/5' UTR to each exon, by RNA_ID
        self.cds_tr_join = self.ex_tr_join.merge(tutrs[[START, STOP, RNA_ID]], on=[RNA_ID])
        self.cds_tr_join = self.cds_tr_join.rename(columns={START+"_y": "t_Start", STOP+"_y": "t_Stop",
                                                            START+"_x": START, STOP+"_x": STOP})
        self.cds_tr_join = self.cds_tr_join.merge(futrs[[START, STOP, RNA_ID]], on=[RNA_ID])
        self.cds_tr_join = self.cds_tr_join.rename(columns={START + "_y": "f_Start", STOP + "_y": "f_Stop",
                                                            START + "_x": START, STOP + "_x": STOP})

        # keep rows where the start base is before the 3'UTR start or where the stop base is after the 5' UTR stop
        # this gives all the exons which are in the CDS, but some might require trimming
        # 555555555555555555555555555                                     333333333333333
        # eeeeeee eeeeeeee   eeeeeeeeee   eeeeeeeeee   eeeeeeeeeeee    eeeeeeee   eeeeeee
        self.cds_tr_join = self.cds_tr_join[((self.cds_tr_join[STRAND] == "+") &
                                            (self.cds_tr_join[START] < self.cds_tr_join["t_Start"]) &
                                            (self.cds_tr_join[STOP] > self.cds_tr_join["f_Stop"])) |
                                            ((self.cds_tr_join[STRAND] == "-") &
                                             (self.cds_tr_join[START] < self.cds_tr_join["f_Start"]) &
                                             (self.cds_tr_join[STOP] > self.cds_tr_join["t_Stop"]))]

        # trim any exons where the start is before the 5'UTR stop, and the stop is more than the 5' UTR stop
        # 555555555555
        #    eeeeeeeeeeeeee <- can't end at or before the 5, as we already excluded that
        self.cds_tr_join[START] = np.where((self.cds_tr_join[STRAND] == "+") &
                                           (self.cds_tr_join[STOP] > self.cds_tr_join["f_Stop"]) &
                                           (self.cds_tr_join[START] < self.cds_tr_join["f_Stop"]),
                                            self.cds_tr_join["f_Stop"] + 1, self.cds_tr_join[START])

        # 33333333333
        #      eeeeeeeeeeee
        self.cds_tr_join[START] = np.where((self.cds_tr_join[STRAND] == "-") &
                                           (self.cds_tr_join[STOP] > self.cds_tr_join["t_Stop"]) &
                                           (self.cds_tr_join[START] < self.cds_tr_join["t_Stop"]),
                                            self.cds_tr_join["t_Stop"] + 1, self.cds_tr_join[START])

        # trim any exons where the stop is after the 3'UTR start, and the start is less than the 3'UTR start
        #        3333333333333
        #     eeeeeeee <- can't start after the 3' UTR start, as we already excluded that
        self.cds_tr_join[STOP] = np.where((self.cds_tr_join[STRAND] == "+") &
                                          (self.cds_tr_join[START] < self.cds_tr_join["t_Start"]) &
                                          (self.cds_tr_join[STOP] > self.cds_tr_join["t_Start"]),
                                           self.cds_tr_join["t_Start"] - 1, self.cds_tr_join[STOP])

        #      555555555555
        #  eeeeeeeeee
        self.cds_tr_join[STOP] = np.where((self.cds_tr_join[STRAND] == "-") &
                                          (self.cds_tr_join[START] < self.cds_tr_join["f_Start"]) &
                                          (self.cds_tr_join[STOP] > self.cds_tr_join["f_Start"]),
                                          self.cds_tr_join["f_Start"] - 1, self.cds_tr_join[STOP])

        self.cds_tr_join.drop_duplicates(subset=[CHROMOSOME, START, STOP, STRAND, RNA_ID], inplace=True)
        self.cds_tr_join = self.cds_tr_join[[CHROMOSOME, SOURCE, TYPE, START, STOP, SCORE,
                                             STRAND, PHASE, RNA_ID, GENE_ID]]

    def _process_cds_to_utrs(self):
        """ Convert CDS data to UTR data. """

        # merge CDS and RNA_IDs so that multiple CDS ids for single RNA are all retained
        TEMP_RNA_ID = "temp_{}".format(RNA_ID)
        temp = self.cds_tr_join
        temp[TEMP_RNA_ID] = temp[RNA_ID]
        temp[RNA_ID] = self.cds_tr_join[RNA_ID] + self.cds_tr_join[ID]

        # Get first and last CDS
        first_cds = self._get_ftrs_by_position(temp, 1)
        last_cds = self._get_ftrs_by_position(temp, -1)

        # Also add in any CDSs which match on RNA_ID and overlap
        # (seems to be the only way to identify multiple CDSs as in the sequence ontology example)
        first_cds = first_cds.rename(columns={STOP: CDS_STOP, START: CDS_START})
        first_cds = self._get_overlaps_by_transcript(self.cds_tr_join, first_cds, START, CDS_START, CDS_STOP)

        # repeat for last CDSs
        last_cds = last_cds.rename(columns={STOP: CDS_STOP, START: CDS_START})
        last_cds = self._get_overlaps_by_transcript(self.cds_tr_join, last_cds, STOP, CDS_START, CDS_STOP)

        five_utrs_pos = self._get_UTR_from_CDS(first_cds[first_cds[STRAND] == '+'], True, TEMP_RNA_ID)
        five_utrs_neg = self._get_UTR_from_CDS(first_cds[first_cds[STRAND] == '-'], False, TEMP_RNA_ID)
        three_utrs_pos = self._get_UTR_from_CDS(last_cds[last_cds[STRAND] == '+'], False, TEMP_RNA_ID)
        three_utrs_neg = self._get_UTR_from_CDS(last_cds[last_cds[STRAND] == '-'], True, TEMP_RNA_ID)

        three_utrs = pd.concat([three_utrs_pos, three_utrs_neg])
        three_utrs = three_utrs.rename(columns={UTR_STOP: THREE_UTR_STOP, UTR_START: THREE_UTR_START})
        three_utrs[TYPE] = self.featureType.three_prime_UTR
        three_utrs[PHASE] = "."

        five_utrs = pd.concat([five_utrs_pos, five_utrs_neg])
        five_utrs = five_utrs.rename(columns={UTR_STOP: FIVE_UTR_STOP, UTR_START: FIVE_UTR_START})
        five_utrs[TYPE] = self.featureType.five_prime_UTR
        five_utrs[PHASE] = "."

        return three_utrs, five_utrs

    def _get_UTR_from_CDS(self, cds, is_left_cds, temp_name):
        """ Build a UTR from a CDS
        :param cds: coding sequence dataframe
        :param is_left_cds: whether this is a CDS at left or right side of gene
        :param temp_name: temporary name for RNA_ID
        :return dataframe of UTRs
        """

        if is_left_cds:
            selection_range = cds[[temp_name, CDS_START]]
            selection_range = selection_range.rename(columns={temp_name: RNA_ID})

            # for left UTR, exon start should be before CDS start, so get exons which are are candidate UTRs
            # by finding all exons which start before or at CDS start point

            utrs = pd.merge(self.ex_tr_join[[CHROMOSOME, SOURCE, SCORE, STRAND, RNA_ID, START, STOP]],
                                 selection_range, on=RNA_ID)
            utrs = utrs.loc[(utrs[START] <= utrs[CDS_START])]

            utrs[UTR_STOP] = np.where(
                (utrs[CDS_START] > utrs[START]) & (utrs[CDS_START] < utrs[STOP]), utrs[CDS_START] - 1, utrs[STOP])
            utrs[UTR_START] = utrs[START]

            # clean up unwanted cols
            utrs.drop([CDS_START, START, STOP], axis=1, inplace=True)
            utrs.drop_duplicates([CHROMOSOME, RNA_ID, STRAND, "UTR_START", "UTR_STOP"], inplace=True)

        else:

            selection_range = cds[[temp_name, CDS_STOP]]
            selection_range = selection_range.rename(columns={temp_name: RNA_ID})

            # get exons which fall in selection range
            # for right hand UTR, exon stop should be after CDS stop, so get exons which are candidate UTRs
            # by finding all exons which stop after or at the CDS stop point
            utrs = pd.merge(self.ex_tr_join[[CHROMOSOME, SOURCE, SCORE, STRAND, RNA_ID, START, STOP]],
                                  selection_range, on=RNA_ID)
            utrs = utrs.loc[(utrs[STOP] >= utrs[CDS_STOP])]

            # if the CDS_STOP is between the exon start and stop, set UTR_START to CDS_STOP+1, otherwise exon start
            utrs[UTR_START] = np.where(
                (utrs[CDS_STOP] > utrs[START]) & (utrs[CDS_STOP] < utrs[STOP]), utrs[CDS_STOP] + 1, utrs[START])
            utrs[UTR_STOP] = utrs[STOP]

            # clean up unwanted cols
            utrs.drop([CDS_STOP, START, STOP], axis=1, inplace=True)
            utrs.drop_duplicates([CHROMOSOME, RNA_ID, STRAND, UTR_START, UTR_STOP], inplace=True)

        return utrs

    def _get_overlaps_by_transcript(self, tr_join, to_overlap, testcol, start, stop):
        """ Get elements in tr_join whose colname entry is in [start,stop] of to_overlap, on the same transcript.
        :param tr_join: a join of transcripts and the elements of interest
        :param to_overlap: dataframe of elements we are seeking overlaps for
        :param testcol: name of column (START or STOP) in tr_join to test
        :param start: name of start column in to_overlap
        :param stop: name of stop column in to_overlap
        """

        self.logger.debug("Getting overlaps by transcript")

        # merge each region with with start/stop details of it's test region
        overlaps = pd.merge(tr_join, to_overlap[[start, stop, RNA_ID]], on=RNA_ID)

        # pull out any rows where the region overlaps with it's test region
        overlaps = overlaps.loc[(overlaps[testcol] >= overlaps[start]) & (overlaps[testcol] <= overlaps[stop]) &
                                (~((overlaps[START] == overlaps[start]) & (overlaps[STOP] == overlaps[stop])))]

        # tidy up columns and add any overlapping region to the set of test regions
        overlaps.drop([start, stop], axis=1, inplace=True)
        overlaps = overlaps.rename(columns={START: start, STOP: stop})
        return pd.concat([to_overlap, overlaps])

    def build_antisense_gtf_gene_only(self, output=None):
        """ Build a gtf with complements of sense/antisense genes.
        :param output: name of output file
        """

        self.logger.info("Building antisense features: this will take several minutes")

        # get extents of entire genes: we will only look for antisense opposite these
        #sense_df = self.featuredict[self.featureType.exon][[CHROMOSOME, SOURCE, TYPE, START, STOP,
        #                                                    SCORE, STRAND, PHASE, GENE_ID, RNA_ID]]
        sense_df = self.sqdict[self.featureType.exon][[CHROMOSOME, SOURCE, TYPE, START, STOP,
                                                            SCORE, STRAND, PHASE, GENE_ID, RNA_ID]]

        gene_extents = sense_df.groupby([GENE_ID, CHROMOSOME, STRAND]).agg({START: min,
                                                                            STOP: max,
                                                                            SOURCE: lambda x: x.iloc[0],
                                                                            TYPE: lambda x: "antisense",
                                                                            SCORE: lambda x: x.iloc[0],
                                                                            PHASE: lambda x: x.iloc[0],
                                                                            RNA_ID: lambda x: ""}).reset_index()
        if gene_extents.empty:
            self.logger.info("Gene extents empty when building antisense features (are gene ids set?)")
            return gene_extents

        gene_extents[CHROMOSOME] = gene_extents[CHROMOSOME].astype(str)

        # account for extents within other extents e.g. miRNA within gene
        # this call builds connected components by chromosome and strand
        # thereby removing extents contained within others
        contigs = gene_extents.groupby([CHROMOSOME, STRAND]).apply(reductionFunction)
        gene_extents = gene_extents.merge(contigs, on=[CHROMOSOME, STRAND, STOP, START])

        # get complement of -ve anti_df and positive sense_df and vice versa
        negs = self._calc_complement(gene_extents, gene_extents, '-', '+')
        poss = self._calc_complement(gene_extents, gene_extents, '+', '-')

        antisense = pd.concat([negs, poss])

        # swap the strand over
        antisense.replace({"+": "temp",  "-": "+", "temp": "-"}, inplace=True)

        # convert start/stop to int - can be set to float if NaNs present
        # antisense[[START, STOP]] = antisense[[START, STOP]].astype(int)
        self._fix_gtf2_col_types(antisense)

        if not (output is None):
            self._output_gtf2(antisense, output)
        else:
            return antisense

    def build_antisense_gtf(self, gtf_data, out_to_file=True):
        """ Build a gtf with complements of sense/antisense genes, considering exons individually
        :param gtf_data: list of tuples: features which will have antisense counts, paired with output gtf names
        e.g [('exon','exon-antisense.gtf'),('intron','intron-antisense.gtf')...]
        :param out_to_file: output result to file
        """

        self.logger.info("Building antisense features")

        merged = None
        sense_df = None
        for feature, output in gtf_data:

            if merged is None:

                # if a read is aligned to an exon/UTR then it can't be counted as antisense
                # so we exclude those regions here (but UTRs are included in exons, so just use exons)
                sense_df = self.featuredict[self.featureType.exon][[CHROMOSOME, SOURCE, TYPE, START, STOP,
                                                                    SCORE, STRAND, PHASE, GENE_ID, RNA_ID]]

                # note pandas runs apply twice on the first group - this is by design to
                # determine whether the shape of the result changes... sigh
                merged = sense_df.groupby([STRAND, CHROMOSOME]).apply(reductionFunction)
                merged.reset_index(inplace=True, drop=True)

            elif feature == self.featureType.intron:
                # if we are looking for antisense counts in introns,
                # also exclude reads aligned to introns on sense strand
                sense_df = \
                    pd.concat([sense_df, self.featuredict[self.featureType.intron]
                    [[CHROMOSOME, SOURCE, TYPE, START, STOP, SCORE, STRAND, PHASE, GENE_ID, RNA_ID]]],
                              axis=0)

                merged = sense_df.groupby([STRAND, CHROMOSOME]).apply(reductionFunction)
                merged.reset_index(inplace=True, drop=True)

            # make sure CHROMOSOME column is string (required for later processing)
            # use featuredict instead of modifying anti_df to avoid SettingWithCopy warning
            merged[CHROMOSOME] = merged[CHROMOSOME].astype(str)
            self.featuredict[feature][CHROMOSOME] = self.featuredict[feature][CHROMOSOME].astype(str)

            # get the regions we want to consider for antisense counts
            anti_df = self.featuredict[feature][[CHROMOSOME, SOURCE, TYPE, START, STOP, SCORE, STRAND, PHASE,
                                                 GENE_ID, RNA_ID]]

            # merge connected components on each individual strand so that when we merge with
            # the opposite strand we have fewer regions to contend with
            self.logger.debug("Merge connected components on each strand")

            # get complement of -ve anti_df and positive sense_df and vice versa
            negs = self._calc_complement(anti_df, merged, '-', '+')
            poss = self._calc_complement(anti_df, merged, '+', '-')

            antisense = pd.concat([negs, poss])

            # convert start/stop to int - can be set to float if NaNs present
            #antisense[[START,STOP]] = antisense[[START, STOP]].astype(int)
            self._fix_gtf2_col_types(antisense)

            if out_to_file:
                self._output_gtf2(antisense, output)
            else:
                return antisense

    def _calc_complement(self, anti_df, merged, strand, opp_strand):
        """
        :param anti_df:
        :param merged: contiguous regions on each strand
        :param strand: string corresponding to strand ('+' or '-')
        :param opp_strand: string corresponding to opposite strand ('+' or '-')
        :return: complement of regions
        """

        self.logger.debug("Calculating complement of strand: {}".format(strand))

        # get contiguous regions of opposite strand only
        merged = merged[merged[STRAND] == opp_strand]

        # get all regions of current strand
        anti_df_strand = anti_df[anti_df[STRAND] == strand]

        # merge current strand with contiguous regions of opposite strand
        combined = pd.concat([merged, anti_df_strand])

        # group regions (both strands) into connected components - same group has same index
        # calculate per chromosome, as scipy's connected components runs out of memory with the full size graph
        indices = combined.groupby([CHROMOSOME]).apply(connected_comp_indices)

        # loop by chromosome, and cut out the regions on the current strand which are masked by regions
        # on the opposite strand
        # slight hacky way to handle chromosome-based indices returned above
        # there should be a more pandas way to do it?
        complement = combined[[]]

        if not indices.empty:

            cs = self.get_chromosomes()
            for c in cs:
                try:
                    self.logger.info("Building antisense for chromosome: " + c + " " + strand + "ve strand")
                    result = combined[combined[CHROMOSOME] == c].\
                        groupby(indices[c]).apply(mask, opp_strand).reset_index()
                except KeyError as e:
                    self.logger.debug("Could not find index while calculating complement of strand: {}"
                                     .format(e.message))
                    continue
                if complement.empty:
                    complement = result
                elif not result.empty:
                    complement = pd.concat([complement, result], axis=0)

        return complement


def mask(x, strand):

    # we want to return -ve strand values if removing strand='+', and +ve strand value if removing strand='-'
    result = x[x[STRAND] != strand]
    candidates = x[x[STRAND] == strand][[START, STOP]]

    if result.empty | candidates.empty:
        return result

    C_START = "Candidate_start"
    C_STOP = "Candidate_stop"
    KEY = "key"
    DELETE = "delete"

    OLD_START = "Old_start"
    OLD_STOP = "Old_stop"
    
    candidates = candidates.rename(columns={
            START: C_START,
            STOP: C_STOP
        })

    result = result.rename(columns={
        START: OLD_START,
        STOP: OLD_STOP
    })

    # possible scenarios:
    # result ------[==============]--------
    # candidate1 ------[================]--   starts after, finishes after: result[STOP] = candidate[START]
    # candidate2 ------[======]------------   starts after, finishes before: result[STOP] = candidate[START] AND
    #                                         result'[START] = candidate[STOP], result'[STOP] = result[STOP]
    # candidate3 [=====================]---   starts before, ends after: delete result
    # candidate4 [============]------------   starts before, ends before: result[START] = candidate[STOP]

    # mark results with index so we can group afterwards
    result.reset_index(level=0, inplace=True)

    # form cartesian product
    result[KEY] = 1
    candidates[KEY] = 1
    result = pd.merge(result, candidates, on=KEY)
    del result[KEY]

    # delete complete overlaps
    result[DELETE] = (result[C_START] <= result[OLD_START]) & (result[C_STOP] >= result[OLD_STOP])
    result = result[result[DELETE] == False]

    # compare each primary strand region with each region on other strand
    # primary region starts before other, set new stop value

    result.loc[result[OLD_START] < result[C_START], STOP] = result[C_START] - 1
    result.loc[result[OLD_STOP] > result[C_STOP], START] = result[C_STOP] + 1
    result.sort_values(by=["index", C_START], inplace=True, ascending="True")

    if not result.empty:
        result = result.groupby("index").apply(adjust_start_stop, OLD_STOP, OLD_START)

    result.reset_index(drop=True, inplace=True)
    return result[[CHROMOSOME, SOURCE, TYPE, START, STOP, SCORE, STRAND, PHASE, GENE_ID, RNA_ID]]


def adjust_start_stop(data, OLD_STOP, OLD_START):
    """ Adjust intervals to account for different possible overlap configurations
    :param data: a  group of data representing a set of intervals
    :param OLD_STOP:
    :param OLD_START:
    :return:
    """

    # Could have: primary strand starts first or second, and overlapping strand starts first or second.
    # This code adjusts for these possibilities so we end up with a set of ordered rows.
    # A value in the START column of the last row indicates that we need to supply the STOP value from
    # primary strand: the primary strand ends the set of regions
    # A value in the STOP column of the first row indicates we need to supply the START value from the primary
    # strand: the primary strand starts the set of regions

    data.reset_index(inplace=True)
    last = len(data)
    added = False

    # test for NaN in START column of last row, and copy the last row to the bottom of the frame if not NaN
    if not np.isnan(data.loc[last-1, START]):
        data.loc[last] = data.loc[last-1]  # add a row at the end
        data.loc[last, START] = np.nan
        data.loc[last, STOP] = data.loc[last, OLD_STOP]
        added = True

    # test for NaN in STOP column of first row, and copy first row to top of the frame if not NaN
    if not np.isnan(data.loc[0, STOP]):
        data.loc[-1] = data.loc[0]  # add a row at the start
        data.loc[-1, START] = data.loc[0, OLD_START]
        data.loc[-1, STOP] = np.nan
        added = True

    if added:
        data.index = data.index + 1   # recalc index to account for -1 we just added
        data = data.sort_index()      # and sort so we can rely on order when we shift STOP column upwards next

    data[STOP] = data[STOP].shift(-1)
    data.dropna(subset=[START, STOP], inplace=True)

    return data


def reductionFunction(data):

    indices = connected_comp_indices(data)

    # group the results by these connected components
    return data.groupby(indices).aggregate({START: 'min', STOP: 'max', CHROMOSOME: 'first', STRAND: 'first'})


def connected_comp_indices(data):

    # create a 2D graph of connectivity between start/stop ranges
    start = data[START].values
    end = data[STOP].values
    chro = data[CHROMOSOME].values

    logger = logging.getLogger()
    logger.debug("Building 2D graph. Data length: {}".format(len(data)))

    graph = (start <= end[:, None]) & (end >= start[:, None]) & (chro == chro[:, None])

    logger.debug("Graph built")

    # find connected components in this graph
    _, indices = connected_components(graph)

    return indices
