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
