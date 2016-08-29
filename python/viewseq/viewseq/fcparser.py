#!/usr/bin/env python

import pandas as pd
import pysam
import logging
import os
import tempfile
import sys
import numpy as np
from collections import namedtuple

from constants import *
from runcommand import RunCommand
from normaliser import Normaliser

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"


# Used to store command line data
class CommandStruct:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


class FeatureCountParser:
    def __init__(self, annotation, fc_exe, fc_options, reads_filenames, old_output_dir=''):
        """ Get feature counts for this annotation and the input bam files.
        :param annotation: the annotation to match reads to
        :param fc_exe: path to featureCounts executable
        :param fc_options: featureCounts options
        :return parsed featureCounts object
        :raise IOError: if there were problem reading the reads files
        """
        self.__logger = logging.getLogger(__name__)
        self.__annotation = annotation

        FeatureType = namedtuple('FeatureType', 'exon intron three_prime_UTR five_prime_UTR antisense')
        self.featureType = FeatureType(exon=annotation.featureType.exon,
                                       intron=annotation.featureType.intron,
                                       three_prime_UTR=annotation.featureType.three_prime_UTR,
                                       five_prime_UTR=annotation.featureType.five_prime_UTR,
                                       antisense=annotation.featureType.antisense)

        # constants defined by featureCounts
        self.__start_col_name = "Start"  # name of Start column in featureCounts output
        self.__stop_col_name = "End"    # name of Stop column in featureCounts output
        self.__geneid_name = "Geneid"   # name of Gene id column in featureCounts output
        self.__chr_name = "Chr"         # name of Chromosome column in featureCounts output
        self.__strand_name = "Strand"   # name of strand column in featureCounts output

        self.__anti_prefix = "antisense"  # prefix for antisense features

        self.__logger.info("Processing alignments: {}".format(" ".join(reads_filenames)))

        # just return if we don't have all bam files
        for readfile in reads_filenames:
            if not os.path.isfile(readfile):
                raise IOError("Could not find reads file: {}".format(readfile))

        # check the chromosome ids match
        if "-A" not in fc_options:
            try:
                self._check_chromosomes(reads_filenames)
            except IOError as e:
                # output this error to the command line as it will mean featureCounts won't be run
                sys.stderr.write(e.message)
                raise

        self.dfs = {}

        ftr_to_cmds = self._build_command_lines(fc_exe, fc_options, reads_filenames)

        self.__logger.info("Running featurecounts")
        r = RunCommand()

        try:
            dfs_by_group = {}

            # make a placeholder for each dataframe for the groups (all or UTRs)
            for c in ftr_to_cmds:
                dfs_by_group[c.ftr] = []

            if len(old_output_dir) == 0:
                r.run_commands([o.cmd for o in ftr_to_cmds])

                self.__logger.info("Reading featurecounts output")
                for c in ftr_to_cmds:
                    # read in the result file
                    self.__logger.debug("Reading featureCounts file: {}".format(c.to_output))
                    dfs_by_group[c.ftr].append(pd.read_csv(c.to_output, header=1, sep='\t', engine='python'))
            else:
                for c in ftr_to_cmds:
                    outfile = "{}/{}-counts.out".format(old_output_dir, c.ftr)
                    dfs_by_group[c.ftr].append(pd.read_csv(outfile, header=1, sep='\t', engine='python'))

            # concat dataframes from same group together
            for group in dfs_by_group.keys():
                self.dfs[group] = pd.concat(dfs_by_group[group])
                self.dfs[group] = self.dfs[group].rename(columns={self.__geneid_name: GENE_ID,
                                                                  self.__chr_name: CHROMOSOME,
                                     self.__start_col_name: START, self.__stop_col_name: STOP,
                                     self.__strand_name: STRAND})

            # now convert the counts features back to the features we are expecting to see
            # join extents with each feature type and copy type into type column
            features = [("all", [self.__annotation.featureType.exon,
                                 self.__annotation.featureType.intron,
                                 self.__annotation.featureType.antisense]),
                        ("UTR", [self.__annotation.featureType.three_prime_UTR,
                                 self.__annotation.featureType.five_prime_UTR])]
            for d, flist in features:
                for ftr in flist:
                    self.dfs[ftr] = self.dfs[d].merge(
                        annotation.featuredict[ftr][[STRAND, START, STOP, GENE_ID]],
                        on=[STRAND, GENE_ID, START, STOP])
                    self.dfs[ftr].drop_duplicates([STRAND, GENE_ID, START, STOP], inplace=True)

        except Exception as e:
            self.__logger.info(e.message)
        finally:
            # tidy up
            self._clean_up(ftr_to_cmds)

    def _clean_up(self, ftr_to_cmds):
        """ Delete temporary files and directories used in running commands in ftr_to_cmds
        :param ftr_to_cmds: list of command line calls and associated output and temporary files
        """

        last_file = None
        for o in ftr_to_cmds:
            for file_to_delete in o.to_delete:
                try:
                    os.remove(file_to_delete)
                    last_file = file_to_delete
                except OSError as e:
                    if e.errno == 2:
                        # no such file, assume already deleted
                        pass

        # we assume all the temp files are in the same temp directory
        try:
            if last_file is not None:
                os.removedirs(os.path.dirname(last_file))
        except OSError as e:
            self.__logger.info(e.message)
            raise

    def _build_command_lines(self, fc_exe, fc_options, reads_filenames):
        """ Build featureCounts specific command lines
        :param fc_exe: featureCounts executable
        :param fc_options: command line options
        :param reads_filenames: bam files containing reads
        :return: list of featureCounts command lines, each as a CommandStruct
        """

        results = []

        # setup a new temp directory
        # put the temp directory in HOME so that if we're running on the cluster
        # the files created by slave processes can be found by the master afterwards
        tempdir = "{}/tempjobs".format(os.environ['HOME'])
        try:
            os.mkdir(tempdir)
        except OSError:
            pass
        temppath = tempfile.mkdtemp(dir=tempdir)

        bam_file_list = " ".join("{}".format(i) for v, i in enumerate(reads_filenames))
        fc_file_path = "{}/{{}}counts.out".format(temppath)
        fc_string = "{} -T 8 -f -t {{}} -a {{}} -o {} {}". \
            format(fc_options, fc_file_path, bam_file_list)

        # export features from the annotation to gtf2, for input to featureCounts
        features = [self.__annotation.featureType.exon,
                    self.__annotation.featureType.intron,
                    self.__annotation.featureType.antisense]

        # first create gtf for exon/intron/antisense - this will be used for exon/intron sense and antisense counts
        # output path for gtf2 file
        results = results + self._build_count_group_command("all", features, temppath, fc_file_path, fc_exe, fc_string)

        # now do the UTRs - need to be separate as they overlap with exons
        features = [self.__annotation.featureType.three_prime_UTR,
                    self.__annotation.featureType.five_prime_UTR,
                    self.__annotation.featureType.CDS,
                    self.__annotation.featureType.antisense]

        results = results + self._build_count_group_command("UTR", features, temppath, fc_file_path, fc_exe, fc_string)

        return results

    def _build_count_group_command(self, group, features, temppath, fc_file_path, fc_exe, fc_string):

        results = []

        gtf_path = "{}/{}".format(temppath, group)
        gtf_paths = self.__annotation.export_to_gtf2(gtf_path, features, "feature")

        for gtf_path in gtf_paths:
            subgroup = os.path.basename(gtf_path)
            subgroup = os.path.splitext(subgroup)[0].split(group, 1)[1]

            output = fc_file_path.format("{}-{}".format(group, subgroup))
            c = CommandStruct(ftr=group,
                              cmd=("{} {}".format(fc_exe, fc_string.format("feature", gtf_path,
                                                                           "{}-{}".format(group, subgroup)))).split(),
                              to_delete=[gtf_path, output, "{}.summary".format(output)],
                              to_output=output)
            results.append(c)
            self.__logger.debug("Appended command {}".format(c.cmd))

        return results

    def _check_chromosomes(self, reads_files):
        """ Check chromosome IDs in current data match those in bam_files
        :param reads_files: list of bam_files whose chromosome IDs we want to match
        :raise IOError: if a reads file is incorrect in some way
        """
        self.__logger.info("Checking chromosome names in annotation and reads files match")

        # run samtools on the bam files
        for reads_file in reads_files:

            try:
                samfile = pysam.AlignmentFile(reads_file, "rb")

                # extract only the @SQ lines from the result
                bam_names = [d["SN"] for d in samfile.header["SQ"]]

            except ValueError:
                raise IOError("Samtools found no header in the reads file {}".format(reads_file))

            gff_names = set(self.__annotation.get_chromosomes())

            # check the names at least intersect
            if not (set(bam_names)).intersection(set(gff_names)):
                raise IOError("The chromosomes in the reads file do not match those in the annotation.\n"
                              "Reads file chromosomes: {}\nAnnotation file chromosomes: {}\n"
                              "Please specify the '-A' option with an aliases file for featureCounts."
                              .format(tuple(bam_names), tuple(gff_names)))

    def counts_per_feature(self, feature, get_antisense=False):
        """ Calculate the raw read counts per feature.
        :param feature: the feature (exon, intron, etc) to use
        :param get_antisense: get antisense counts rather than sense counts for this feature
        :return: list of counts per feature or empty list if feature is not found
        :raise ValueError if feature is not a FeatureType or dataframe not found
        """
        name_of_frame = self._check_feature_type(feature, get_antisense)

        # sum across each column
        cols = [col for col in self.dfs[name_of_frame].columns if ".bam" in col]
        counts = self.dfs[name_of_frame][cols].sum(axis=1)
        return counts.tolist()

    def reads_in_first_feature(self, feature, get_antisense=False):
        """ Calculate counts of reads in first features.
        :param feature: the feature (exon, intron, etc) to use
        :param get_antisense: get antisense counts rather than sense counts for this feature
        :return List of read counts in first exons
        :raise ValueError if feature is not a FeatureType
        """

        name_of_frame = self._check_feature_type(feature, get_antisense)

        result = self.__annotation.filter_by_first_features(self.dfs[name_of_frame], feature,
                                                            self.__start_col_name, self.__stop_col_name)


        result.drop([START, STOP, STRAND], axis=1, inplace=True)

        cols = [col for col in self.dfs[name_of_frame].columns if ".bam" in col]
        result["counts"] = result[cols].sum(axis=1)
        return result[[GENE_ID, "counts"]]

    def reads_in_last_feature(self, feature, get_antisense=False):
        """ Calculate counts of reads in last features.
        :param feature: the feature (exon, intron, etc) to use
        :param get_antisense: get antisense counts rather than sense counts for this feature
        :return List of read counts in last exons
        :raise ValueError if feature is not a FeatureType
        """

        name_of_frame = self._check_feature_type(feature, get_antisense)

        try:
            result = self.__annotation.filter_by_last_features(self.dfs[name_of_frame], feature,
                                                               self.__start_col_name, self.__stop_col_name)
        except KeyError:
            raise ValueError("No matching data was found for feature: {}".format(name_of_frame))

        result.drop([START, STOP, STRAND], axis=1, inplace=True)

        cols = [col for col in self.dfs[name_of_frame].columns if ".bam" in col]
        result["counts"] = result[cols].sum(axis=1)
        return result[[GENE_ID, "counts"]]

    def _check_feature_type(self, feature, get_antisense):
        """
        :param feature: a feature to check
        :return: name of dataframe for feature
        :raise ValueError if the feature is not a FeatureType
        """

        # originally set up to allow for multiple antisense dataframes, e.g. antisense-exon, antisense-intron etc.
        # but currently only one is supported

        if get_antisense:
            name_of_frame = "{}".format(self.__anti_prefix)
        else:
            name_of_frame = feature

            if feature not in self.featureType._fields:
                raise ValueError("Feature type {} is not valid".format(feature))

        if name_of_frame not in self.dfs.keys():
            raise ValueError("Dataframe {} does not exist".format(name_of_frame))

        return name_of_frame

    def _create_feature_data_by_gene(self, name, ftr):
        """ Create summary of reads for ftr by gene
        :param name: name of column
        :param ftr: feature type of column
        :return: dataframe containing reads summed by gene
        """

        # get reads columns
        cols = [col for col in self.dfs[ftr].columns if ".bam" in col]
        # sum reads counts over each gene
        sumreads = self.dfs[ftr].groupby(GENE_ID).sum()
        # sum across all reads columns
        sumreads[name] = sumreads[cols].sum(axis=1)

        return sumreads

    def create_data_by_gene(self, features, total_mapped_reads, extras):
        """ Build a dataframe containing the supplied features read counts, by gene, for sense and antisense
        :param features: list of features to add to dataframe
        :param total_mapped_reads: total number of assigned reads (for FPKM calc)
        :param extras: list of extra dataframes to add, must be by gene already
        :return: dataframe
        """

        cols = []
        emptycols = []

        normer = Normaliser(self.__annotation, total_mapped_reads)
        for ftr in features:

            try:
                col = self._create_feature_data_by_gene(ftr, ftr)
                cols.append(col[[ftr]])

                fpkm_colname = "{}-FPKM".format(ftr)
                fpkm_col = normer.normalise_fpkm(col[ftr], ftr, fpkm_colname)

                cols.append(fpkm_col[[fpkm_colname]])
            except KeyError as e:
                self.__logger.debug("A dataframe could not be found for one of the supplied feature names: {}"
                                   .format(e.message))
                emptycols.append(ftr)
                emptycols.append("{}-FPKM".format(ftr))  # store names of columns we have no data for
                continue

        # merge the whole lot together
        if len(cols) > 0:
            result = cols[0]
            for i, frame in enumerate(cols[1:], 2):
                result = pd.concat([frame, result], axis=1)

            for frame in extras:
                result = pd.concat([frame, result], axis=1)
        else:

            if len(extras) > 0:
                result = extras[0]
                for frame in extras[1:]:
                    result = pd.concat([frame, result], axis=1)

            else:
                raise ValueError("Could not create gene-level data,"
                                 " as no matching data was found for features: ".join(features))

        for name in emptycols:
            result[name] = np.nan

        return result
