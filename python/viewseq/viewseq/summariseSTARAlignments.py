#!/usr/bin/python

import re
import json
import os
import logging

from six import itervalues  # for python 2/3 compatibility

__author__ = "Nick Schurch, Kira Mourao"
__email__ = "n.schurch@dundee.ac.uk, k.mourao@dundee.ac.uk"

# Constant string values used as keys
AXIS_LABEL = "axis_label"
SEARCH_STRING = "search_string"
LIMIT = "limit"
SCALING = "scaling"


class STARLogScraper:

    def __init__(self):
        """
        Constructor.
        """

        self.__logger = logging.getLogger(__name__)

        # Set up STAR log fields to scrape
        # axis label: label for (ultimate) plot
        # search string: string to search for in STAR log
        # limit: default axis limit
        # scaling: amount to scale values by, None if no scaling
        self.__fields = {
              1: {
                  AXIS_LABEL: "No. Mreads",
                  SEARCH_STRING: "Number of input reads",
                  LIMIT: 100,
                  SCALING: 1E6
                  },
              2: {
                  AXIS_LABEL: "Uniquely-mapping (%)",
                  SEARCH_STRING: "Uniquely mapped reads %",
                  LIMIT: 100,
                  SCALING: None
                  },
              3: {
                  AXIS_LABEL: "Multi-mapping (%)",
                  SEARCH_STRING: "% of reads mapped to multiple loci",
                  LIMIT: 100,
                  SCALING: None
                  },
              4: {
                  AXIS_LABEL: "Unmapped - multi-mapping (%)",
                  SEARCH_STRING: "% of reads mapped to too many loci",
                  LIMIT: 10,
                  SCALING: None
                  },
              5: {
                  AXIS_LABEL: "Unmapped - mismatches (%)",
                  SEARCH_STRING: "% of reads unmapped: too many mismatches",
                  LIMIT: 10,
                  SCALING: None
                  },
              6: {
                  AXIS_LABEL: "Unmapped - too short (%)",
                  SEARCH_STRING: "% of reads unmapped: too short",
                  LIMIT: 10,
                  SCALING: None
                  },
              7: {
                  AXIS_LABEL: "Unmapped - other (%)",
                  SEARCH_STRING: "% of reads unmapped: other",
                  LIMIT: 10,
                  SCALING: None
                  }
              }

    def scrape(self, rootpath, outfile, logfilename='Log.final.out'):
        """ Scrape STAR logs named as logfilename, found in rootpath, and output to outfile
        :param rootpath: the root directory to start the file search from. All subdirectories will be scanned.
        :param outfile: the name of the output json file
        :param logfilename: the STAR log filename to look for. Defaults to the standard 'Log.final.out'.
        """

        self.__logger.info("Scraping STAR logs")

        # find STAR log files
        matched_files = self._regex_find_files(rootpath, logfilename)

        # set limits to defaults predefined in fields
        limits = {}
        for field in itervalues(self.__fields):
            limits[field[AXIS_LABEL]] = field[LIMIT]

        # pull out data from each matched file, storing in filedata[filename]
        filedata = {}       # data corresponding to each axis label
        for filename in matched_files:
            self.__logger.debug("Processing file: {}...".format(filename))

            # keep rep and condition in filename in case reps have same names in different conditions
            json_filename = os.path.relpath(filename, rootpath)
            filedata[json_filename], limits = self._scrape_star_log(filename, limits)

        # remove axes where all values are 0
        axes = self._find_zero_axes(filedata)
        for axis in axes:
            for filename in filedata:
                del filedata[filename][axis]
            del limits[axis]

        # make the json data structure and output
        outdata = {"data": filedata,
                   "limits": limits,
                   "conditions": filedata.keys()}

        self.__logger.debug("Writing data to output file: {}...".format(outfile))
        with open(outfile, "w") as f:
            json.dump(outdata, f, sort_keys=True, indent=4)

        self.__logger.info("Finished scraping STAR logs")

    def _find_zero_axes(self, filedata):
        """ Find any axis with all values exactly zero
        :param filedata: data corresponding to each field for each log
        :return: list of axes to be removed if any
        """

        axes = []
        for field in itervalues(self.__fields):
            allzero = True
            for filename in filedata:
                if filedata[filename][field[AXIS_LABEL]] != 0:
                    allzero = False
                    break
            if allzero:
                axes.append(field[AXIS_LABEL])

        return axes

    def _scrape_star_log(self, filename, limits):
    
        """ This subroutine reads and scrapes a STAR log for the fields given.
        :param filename: name of log to scrape
        :param limits: current limits on each axis
        :returns: dictionary of data corresponding to each field, updated limits
        """

        self.__logger.debug("Scraping log file {}".format(filename))

        with open(filename, "r") as starlogdata:

            fielddata = {}
            for line in starlogdata:
                linedata = line.strip().split("|")

                # check if this line of the file matches any of the search strings
                for field in itervalues(self.__fields):
                    if linedata[0].strip() == field[SEARCH_STRING]:
                        thisdata = float(linedata[1].strip().split("%")[0])

                        # apply scaling
                        if field[SCALING] is not None:
                            thisdata = thisdata / field[SCALING]

                        # check limits
                        if thisdata > field[LIMIT]:
                            limits[field[AXIS_LABEL]] = int(thisdata * 1.05)

                        fielddata[field[AXIS_LABEL]] = thisdata
    
        return fielddata, limits

    def _regex_find_files(self, path, regex):

        """ Find all files matching the regex string in the path
        :param path: path to search
        :param regex: regular expression to match
        :return set of files which match regex
        """

        self.__logger.debug("Locating files...")

        matcher = re.compile(regex+"$")

        matching_files = set()

        for root, dirs, files in os.walk(path):
            for filename in files:
                if matcher.match(filename):
                    matching_files.add(os.path.join(root, filename))

        self.__logger.debug("Found {} files.".format(len(matching_files)))

        return matching_files
