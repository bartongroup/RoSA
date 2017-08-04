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

from __future__ import print_function
import logging
import subprocess
import tempfile
import os
import json
import datetime
import sys


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

try:
    import drmaa
except RuntimeError:
    eprint("No drmaa library: Running in serial mode.")  # yuk, but logger not initialised yet...

__author__ = "Kira Mourao"
__email__ = "k.mourao@dundee.ac.uk"


class RunCommand:
    def __init__(self):
        self.__logger = logging.getLogger(__name__)

    def run_commands(self, cmds, fnames=None, output_file=False, output_values=False, forceSerial=False):
        """ Run executable and arguments given in cmd
        :param cmds list of commands (exe + command line args) to run
        :param output_file: whether each command outputs to file or not
        :param output_values: whether each command outputs a string or not (mutually excludes output_file)
        :param forceSerial: force the call to use serial mode
        :param fnames: names for output file
        :return result of call as list of output files
        """
        assert(not(output_file & output_values))  # only one of output_file and output_values may be set

        # initialise fnames here,not in param list (default value is mutable)
        if fnames is None:
            fnames = []

        try:
            if forceSerial:
                result = self._run_serial(cmds, fnames, output_file, output_values)
            else:
                result = self._run_parallel(cmds, fnames, output_file, output_values)
        except NameError:
            result = self._run_serial(cmds, fnames, output_file, output_values)

        return result

    def _run_serial(self, cmds, fnames, output_file=False, output_values=False):
        """ Run each command in cmds serially
        :param cmds: list of commands to be run, where each command is a list of strings
        :param output_file: whether each command outputs to file or not
        :param output_values: whether each command outputs a string or not (mutually excludes output_file)
        :return: list of outputs corresponding to each command
        """

        outputs = []

        tempdir = tempfile.gettempdir()

        # loop over each command, use the loop index for naming the output file
        for index, cmd in enumerate(cmds):

            self.__logger.debug("Running {}".format(cmd))
            try:
                if output_file:
                    # create an output file name
                    if len(fnames) == 0:
                        outfile = os.path.join(tempdir,"{}out{}.txt".format(
                            datetime.datetime.now().isoformat().replace(":", "."), index))
                    else:
                        outfile = os.path.join(tempdir, fnames[index])

                    # open the file and call the command
                    f = open(outfile, "w")
                    subprocess.check_call(cmd, stdout=f)
                    f.close()

                    # append the output file to the list of outputs
                    outputs.append(outfile)
                elif output_values:

                    try:
                        output = subprocess.check_output(cmd)
                        outputs.append(json.dumps(json.loads(output)))
                    except ValueError:
                        # output was empty
                        self.__logger.info("Appended empty output as result of command was empty: {}".format(cmd))
                        outputs.append("")

                else:
                    subprocess.check_call(cmd)
            except OSError as e:
                self.__logger.info("Error when attempting to run {}: {}".format(cmd[0], e.message))
                raise e
            except subprocess.CalledProcessError as e:
                self.__logger.info("Error when attempting to run {}: {}".format(cmd[0], e.message))
                raise e

        return outputs

    def _run_parallel(self, cmds, fnames, output_file=False, output_values=False):
        """ Run list of commands in parallel, accounting for high IO behaviour
        :param cmds list of commands (exe + command line args) to run
        :param output_file: whether each command outputs to file or not
        :param output_values: whether each command outputs a string or not (mutually excludes output_file)
        :return outputs list of outputs corresponding to each commmand
        """

        out_files = []
        outputs = []

        hometempdir = "{}/tempjobs".format(os.environ['HOME'])
        try:
            os.mkdir(hometempdir)
        except OSError:
            pass
        hometemppath = tempfile.mkdtemp(dir=hometempdir)

        with drmaa.Session() as s:
            self.__logger.debug("Creating job template")
            jt = s.createJobTemplate()
            jobid = -1

            for index, cmd in enumerate(cmds):

                if output_file | output_values:
                    # create an output file name
                    if len(fnames) == 0:
                        outfile = os.path.join(hometemppath,"{}out{}.txt".format(
                            datetime.datetime.now().isoformat().replace(":", "."), index))
                    else:
                        outfile = os.path.join(hometemppath, fnames[index])
                    # add the file name to the list of outputs
                    out_files.append(outfile)
                    # specify the output file for this cmd
                    jt.nativeSpecification = "-cwd -V -o {}".format(outfile)  # qsub arguments
                else:
                    jt.nativeSpecification = "-cwd -V"  # qsub arguments

                jt.remoteCommand = cmd[0]  # the command to run: note modules required need to be loaded
                jt.args = cmd[1:]  # the command line args
                jt.joinFiles = True  # join qsub output and error files together

                jobid = s.runJob(jt)
                self.__logger.debug("Job submitted with ID {}".format(jobid))

            # wait for every job to complete
            s.synchronize([drmaa.Session.JOB_IDS_SESSION_ALL], drmaa.Session.TIMEOUT_WAIT_FOREVER, True)

            self.__logger.debug("Removing job template {}".format(jobid))
            s.deleteJobTemplate(jt)

        if output_values:
            # open each file and put the contents into a list
            # assumes output is in JSON format
            for filename in out_files:
                self.__logger.debug("Reading output value file: {}".format(filename))
                with open(filename, 'r') as f:
                    try:
                        outputs.append(json.dumps(json.load(f)))
                    except ValueError:
                        # file f was empty
                        self.__logger.debug("Appended empty output as file {} was empty".format(filename))
                        outputs.append("")
        else:
            # copy the output file to its designated location
            for outfile in out_files:
                homeoutfile = os.path.join(hometemppath, os.path.basename(outfile))
                outputs.append(homeoutfile)

        return outputs
