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

""" Script to make antisense annotation """

import argparse
import datetime
import logging

from parser_factory import ParserFactory

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"


def main():

    logname = "makeannotation-log-{}.txt".format(datetime.datetime.now().isoformat())
    logging.basicConfig(filename=logname, filemode='w', level=logging.DEBUG,
                        format='%(asctime)s %(message)s')  # include timestamp

    # parse arguments, which should be an annotation file
    argparser = argparse.ArgumentParser(description="A script to build an antisense annotation from a sense annotation")
    required_named = argparser.add_argument_group('required named arguments')
    required_named.add_argument('-a', action='store', dest='annot_file',
                                help="Specify the location of the annotation file.", required=True)

    required_named.add_argument('-o', action='store', dest='output_file',
                                help="Specify the name of the output file.", required=True)

    args = argparser.parse_args()

    logger = logging.getLogger(__name__)
    logger.info("Creating annotation file parser")
    parser = ParserFactory.create_parser(args.annot_file)

    logger.info("Building antisense annotation")
    parser.build_antisense_gtf_gene_only(args.output_file)


if __name__ == "__main__":
    main()

