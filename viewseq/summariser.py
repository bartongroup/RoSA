#!/usr/bin/env python

"""Summarise annotation and read count data."""


import logging.config
import datetime
import os.path
import argparse
import json
import locale
import yaml
import shutil as su
import inspect

from parser_factory import ParserFactory
from fcparser import FeatureCountParser
from samstats import SamStats
from summariseSTARAlignments import STARLogScraper


__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"


class InputData:
    """Holds and controls access to the input data."""

    __annotation_exts = [".gff", ".gtf"]  # known extensions of annotation files
    __reads_exts = [".bam"]  # known extensions of reads files
    __hist_exts = [".bigwig"]  # known extensions of histogram files

    def __init__(self):
        """Return an empty InputData object"""

        self.__logger = logging.getLogger(__name__)

        self.annotation_filename = ""
        self.reads_filenames = []
        self.hist_filename = ""

    def add_input_file(self, input_file):
        """Add an input file to the input data.
        :param input_file: an input file in the format gff3, bam or wig
        :raises IOError if input_file does not exist
        :raises UnknownFileExtError if it does not know how to handle files with the input_data file extension
        """

        # check the file exists - unpythonic but don't really want to get part
        # way through summarisation and then find we're missing a file
        if not os.path.exists(input_file):
            raise IOError("{} not found".format(input_file))

        # determine the file type from the file extension
        # annotation, reads or histogram
        _, extension = os.path.splitext(input_file)
        if extension in self.__annotation_exts:
            self.annotation_filename = input_file
        elif extension in self.__reads_exts:
            self.reads_filenames.append(input_file)
        elif extension in self.__hist_exts:
            self.hist_filename = input_file
        else:
            raise self.UnknownFileExtError("Unknown file extension: {}".format(extension))

    def get_annotation(self):
        """Get the annotation corresponding to the annotation data file."""

        if self.annotation_filename == "":
            return None
        else:
            return ParserFactory.create_parser(self.annotation_filename)

    def get_feature_counts(self, annotation, fc_exe, fc_options, fc_outputs=None):
        """ Get feature counts for this annotation and the input bam files.
        :param annotation: the annotation to match reads to
        :param fc_exe: path to featureCounts executable
        :param fc_options: featureCounts options string
        :param fc_outputs: list of old outputs, if any
        """

        if fc_outputs is None:
            fc_outputs = []

        # pass through to fcparser to create featurecounts file etc
        if len(self.reads_filenames) == 0:
            return None
        try:
            fcinput = FeatureCountParser(annotation, fc_exe, fc_options, self.reads_filenames, fc_outputs)
        except IOError as e:
            self.__logger.info(e.message)
            fcinput = None

        return fcinput

    def get_reads_data(self):
        """ Get reads counts by strand and by mismatch.
        :return: sam mismatch counts
        """

        if len(self.reads_filenames) == 0:
            return None

        saminput = SamStats(self.reads_filenames)

        return saminput

    class UnknownFileExtError(Exception):
        def __init__(self, value):
            self.value = value


class Stats:
    total_reads = 0
    reads_per_feature = []
    first_feature_reads = []
    last_feature_reads = []


class Summary:
    """Builds summary statistics of the input data."""

    __summarisation_string = "Data summary for the following files:"

    def __init__(self, inputs, fc_exe, fc_options, alignment_log_params=None, old_files=None):
        """Return a Summary object containing summary statistics of the inputs."""

        self.__logger = logging.getLogger(__name__)
        self.__logger.info("Building Summary")

        locale.setlocale(locale.LC_ALL, 'en_GB.UTF-8')

        if alignment_log_params is None:
            alignment_log_params = {}

        assert isinstance(inputs, InputData)
        self.inputdata = inputs
        self.a = None
        self.f = None
        self.s = None

        # set featureCounts executables
        self._fc_exe = fc_exe
        self._fc_options = fc_options

        self.alignment_log_params = alignment_log_params

        # set up placeholders for stats
        self.exon_stats = Stats()
        self.intron_stats = Stats()
        self.three_utr_stats = Stats()
        self.five_utr_stats = Stats()
        self.antisense_stats = Stats()

        self.reads_stats = []

        self.exon_lengths = []
        self.intron_lengths = []
        self.exon_count_per_gene = []

        # for debugging only really
        # re-use old temp files generated by previous run

        self.old_fc_outputs = []
        self.old_sam_outputs = []
        self.old_samcounts = []

        if old_files is not None:
            with open(old_files, 'r') as f:
                l = f.readline().strip()
                if l:
                    self.old_fc_outputs = l
                l = f.readline().strip()
                if l:
                    self.old_sam_outputs = l.split(",")
                l = f.readline().strip()
                if l:
                    self.old_samcounts = l.split(",")

        # run the stats calculations
        self.__calc_annotation_stats()

    def __calc_annotation_stats(self):
        """Calculate statistics for the annotation and reads data."""

        self.__logger.info("Calculating statistics...")

        self.a = self.inputdata.get_annotation()
        self.f = self.inputdata.get_feature_counts(self.a, self._fc_exe, self._fc_options, self.old_fc_outputs)
        self.s = self.inputdata.get_reads_data()

        if self.a is not None:
            # get exon lengths
            self.exon_lengths = self.a.get_ftr_lengths(self.a.featureType.exon)

            # get exon count/gene
            self.exon_count_per_gene = self.a.get_exons_per_gene()

            self.transcript_count_per_gene = self.a.get_transcripts_per_gene()

            # get intron lengths
            try:
                self.intron_lengths = self.a.get_ftr_lengths(self.a.featureType.intron)
            except ValueError:
                # if the introns calc failed there will be no intron features
                # continue with empty intron lengths list
                pass

        if self.s is not None:
            self.reads_stats = self.s.get_stats()

    def output_for_visualisation(self, output_dir):
        """
        Output raw data for visualisation as scatter plots etc.
        :param output_dir: Output directory name
        """

        self.__logger.info("Outputting visualisation data")

        total_mapped_reads = locale.atoi(self.reads_stats["forward strand"][0][1] +
                                         self.reads_stats["reverse strand"][0][1])
        exon_lengths = self.a.get_ftr_lengths(self.a.featureType.exon, bygene=True)
        transcript_counts = self.transcript_count_per_gene
        ave_exon_lengths = exon_lengths.div(transcript_counts.transcript_counts, axis='index')
        gene_names = self.a.get_gene_names()

        # clean up the output directory and recreate
        try:
            su.rmtree(output_dir)
            os.makedirs(output_dir)
        except OSError as e:
            if e.errno == 2:
                pass   # output_dir did not exist

        # put the summary html file, css and scripts into the output directory
        dirname = os.path.dirname(os.path.abspath(inspect.stack()[0][1]))
        parent_dir = os.path.dirname(dirname)

        su.copytree(os.path.join(parent_dir, "viewseq", "websummary", "scripts"),
                    os.path.join(output_dir, "scripts"))
        su.copytree(os.path.join(parent_dir, "viewseq", "websummary", "css"),
                    os.path.join(output_dir, "css"))
        su.copy(os.path.join(parent_dir, "viewseq", "websummary", "summary.html"), output_dir)

        try:
            result = self.f.create_data_by_gene([self.f.featureType.exon, self.f.featureType.intron,
                                                 self.f.featureType.three_prime_UTR, self.f.featureType.five_prime_UTR,
                                                 self.f.featureType.antisense],
                                                total_mapped_reads,
                                                [ave_exon_lengths, self.exon_count_per_gene,
                                                 transcript_counts, gene_names])

            self.__logger.info("Outputting genes JSON file")
            outfile = os.path.join(output_dir, "genes.json")
            result.fillna(0)
            result.reset_index().to_json(path_or_buf=outfile, orient="split", double_precision=4)
        except ValueError as e:
            self.__logger.debug("ValueError when creating genes.json: {}".format(e.message))
            pass

        try:
            self.__logger.info("Outputting first exons JSON file")
            fresult = self.f.reads_in_first_feature(self.f.featureType.exon)
            outfile = os.path.join(output_dir, "firstexons.json")
            fresult.fillna(0)
            fresult.to_json(path_or_buf=outfile, orient="split", double_precision=4)
        except ValueError:
            pass

        try:
            self.__logger.info("Outputting last exons JSON file")
            lresult = self.f.reads_in_last_feature(self.f.featureType.exon)
            outfile = os.path.join(output_dir, "lastexons.json")
            lresult.fillna(0)
            lresult.to_json(path_or_buf=outfile, orient="split", double_precision=4)
        except ValueError:
            pass

        self.__logger.info("Outputting annotation features JSON files")
        results = self.a.get_data_by_feature()
        for ftr, result in results:
            outfile = os.path.join(output_dir, "{}-data.json".format(ftr))
            result.to_json(path_or_buf=outfile, orient="split", double_precision=4)

        self.__logger.info("Outputting totals JSON file")
        result = self.a.get_metadata()
        outfile = os.path.join(output_dir, "totals.json")
        result.fillna(0)
        result.to_json(path_or_buf=outfile, orient="records")

        # input file names to json
        names = {"annotation": self.inputdata.annotation_filename,
                 "reads_files": self.inputdata.reads_filenames}

        self.__logger.info("Outputting metadata JSON file")
        with open(os.path.join(output_dir, "metadata.json"), "w") as f:
            json.dump(names, f)

        self.__logger.info("Outputting matches JSON file")
        with open(os.path.join(output_dir, "matches.json"), "w") as f:
            json.dump(self.reads_stats, f)

        # get the STAR logs data
        if self.alignment_log_params:
            s = STARLogScraper()
            s.scrape(self.alignment_log_params["rootpath"],
                     os.path.join(output_dir, "logsummaries.json"),
                     self.alignment_log_params["logfilename"])


def build_summary(annotation_data, read_data, hist_data, output="summary", featurecounts_exe="featureCounts",
                  fc_options="", alignment_log_params=None, old_files=None):
    """Read in the data files and output summary information.
    :param annotation_data: the annotation data file (gff)
    :param read_data: a list of reads data files (bam)
    :param hist_data: the histograms data file (wig)
    :param output: summary will be written to directory output
    :param featurecounts_exe: optional path to, and name of, featureCounts executable
    :param fc_options: options string for featureCounts

    :param alignment_log_params: option set of params for script to analyse STAR alignment logs
    :param old_files: list of old temp files to use instead of regenerating
    """

    # set up logging
    logger = logging.getLogger(__name__)
    logger.info("Initialising summary")

    # initialise alignment log params
    if alignment_log_params is None:
        alignment_log_params = {}

    # get the data files
    d = InputData()
    logger.info("Adding annotation file: {}".format(annotation_data))
    d.add_input_file(annotation_data)
    for read_file in read_data:
        logger.info("Adding read file: {}".format(read_file))
        d.add_input_file(read_file)
    #d.add_input_file(hist_data)

    # run summarisation
    s = Summary(d, featurecounts_exe, fc_options, alignment_log_params)
    s.output_for_visualisation(output)

    logger.info("Summary complete. Results can be browsed in {}.".format(os.path.join(output, "summary.html")))


def run_main(args):

    dirname = os.path.dirname(os.path.abspath(inspect.stack()[0][1]))
    parent_dir = os.path.dirname(dirname)

    logconfig_path = os.path.join(parent_dir, "viewseq", "configdata", "logging.yaml")

    with open(logconfig_path) as f:
        logconfig = yaml.load(f)

        # Append the date stamp to the file name
        log_filename = logconfig['handlers']['file']['filename']
        base, extension = os.path.splitext(log_filename)
        justnow = datetime.datetime.now().isoformat()
        log_filename = '{}{}{}'.format(base, justnow, extension)
        logconfig['handlers']['file']['filename'] = log_filename

        logging.config.dictConfig(logconfig)

    # parse arguments, which should just be a config file
    argparser = argparse.ArgumentParser(description="A program to summarise annotation and read counts data")
    argparser.add_argument('config_file', help="The configuration file.")

    args = argparser.parse_args(args)

    with open(args.config_file, 'r') as yaml_file:
        config = yaml.safe_load(yaml_file)

        # expected fields
        annotation = "annotation"
        readsfiles = "readsfiles"
        outputdir = "outputdir"

        fc_exe = "featureCounts_exe"
        fc_cmd_line = "featureCounts_options"

        star_rootpath = "STAR_rootpath"
        star_logfilename = "STAR_logfilename"

        build_summary(config[annotation], config[readsfiles], "",
                      config[outputdir], config[fc_exe],
                      config[fc_cmd_line],
                      {"rootpath": config[star_rootpath], "logfilename": config[star_logfilename]}, None)


def main():
    import sys
    run_main(sys.argv[1:])

if __name__ == "__main__":
    import sys
    run_main(sys.argv[1:])
