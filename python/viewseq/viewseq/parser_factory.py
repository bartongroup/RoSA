#!/usr/bin/env python

import os
from gffparser import GffParser
from gtfparser import GtfParser

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"


class ParserFactory:
    """ Implementation of a simple factory class to create Parser objects. """

    @staticmethod
    def create_parser(input_file):
        """ Create a parser appropriate for the input_file.
        :param input_file: a file for which a parser is required
        :return: a parser matching the file extension
        """
        _, extension = os.path.splitext(input_file)
        parser_objects = {'.gff': GffParser, '.gtf': GtfParser}
        return parser_objects[extension](input_file)
