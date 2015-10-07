#!/usr/bin/python
################################################################################
#
# sageWiggle.py was written by J.F.J. Laros at Leiden University Medical Center
# and released under the MIT license:
# https://git.lumc.nl/j.f.j.laros/piletools/raw/master/LICENSE
#
# Then updated by I.F.A.C. Fokkema at Leiden University Medical Center to be
# compatible to Samtools 0.1.19 and up.
#
# Created     : 2011-10-06
# Modified    : 2015-10-06
# Version     : 0.5
#
# Copyright   : 2011-2015 Leiden University Medical Center; http://www.LUMC.nl/
# Programmers : Dr. Jeroen F.J. Laros <J.F.J.Laros@LUMC.nl>
#               Ing. Ivo F.A.C. Fokkema <I.F.A.C.Fokkema@LUMC.nl>
#
# Changelog   : 0.5     2015-10-06
#               Samtools 0.1.19 and higher creates mpileup files that have
#               positions with 0 coverage, leaving out other columns and causing
#               index errors.
#
#
# This work is licensed under the MIT license.
# Copyright (c) 2011-2015 Leiden University Medical Center <humgen@LUMC.nl>
#                         Jeroen F.J. Laros <J.F.J.Laros@LUMC.nl>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
##############

"""
@requires: sys
@requires: re
@requires: argparse
"""

import argparse # ArgumentParser().
import sys      # stdin.
import re       # finditer().
import os       # basename(), splitext().

def findall(string, substring) :
    """
    Find the index of all occurences of a substring in a string.

    @arg string:
    @type string: string
    @arg substring:
    @type substring: string

    @returns: List of indices.
    @rtype: list[int]
    """

    occurences = []
    for i in re.finditer(substring, string) :
        occurences.append(i.start())
    return occurences
#findall

def makeWiggles(pileupHandle, forwardHandle, reverseHandle, name, threePrime) :
    """
    Scan a pileup file for start of reads, determine whether the read is in
    forward or reverse orientation and store this information in either the
    forward wiggle file of the reverse wiggle file, dependant on the
    orientation.

    @arg pileupHandle: Open handle to the pileup file.
    @type pileupHandle: file handle
    @arg forwardHandle: Open handle to the forward wiggle file.
    @type forwardHandle: file handle
    @arg reverseHandle: Open handle to the reverse wiggle file.
    @type reverseHandle: file handle
    @arg name: Short name of the sample.
    @type name: string
    @arg threePrime: Record the 3' end of the read instead of the 5' end.
    @arg threePrime: bool
    """

    # Headers for the start of the wiggle file and for chromosomes.
    wiggleHeader = \
        "track type=wiggle_0 name=%s description=%s visibility=full\n"
    chromosomeHeader = "variableStep chrom=%s\n"
    forwardWiggleName = name + "_forward"
    reverseWiggleName = name + "_reverse"

    # Write the header to each wiggle file.
    forwardHandle.write(wiggleHeader % (forwardWiggleName,
        forwardWiggleName))
    reverseHandle.write(wiggleHeader % (reverseWiggleName,
        reverseWiggleName))

    chromosome = ""
    line = pileupHandle.readline()
    while line :
        data = line.split()
        # Samtools 0.1.19 and higher creates mpileup files that have positions
        # with 0 coverage, leaving out other columns and causing index errors.
        if len(data) < 5:
            line = pileupHandle.readline()
            continue
        #if

        if chromosome != data[0] : # We found a new chromosome.
            chromosome = data[0]
            forwardHandle.write(chromosomeHeader % chromosome)
            reverseHandle.write(chromosomeHeader % chromosome)
        #if
        position = data[1]
        pile = data[4]

        forwardDepth = 0
        reverseDepth = 0
        for i in findall(pile, '\^') :   # Find all starts of reads.
            forwardMapped = False
            if pile[i + 2] in ".ACGTN" : # Is it mapped to the forward strand?
                forwardMapped = True
            if not threePrime and forwardMapped :
                forwardDepth += 1
            if threePrime and not forwardMapped :
                reverseDepth += 1
        #for

        for i in findall(pile, '\$') :   # Find all ends of reads.
            reverseMapped = False
            if pile[i - 1] in ",acgtn" : # Is it mapped to the reverse strand?
                reverseMapped = True
            if not threePrime and reverseMapped :
                reverseDepth += 1
            if threePrime and not reverseMapped :
                forwardDepth += 1
        #for

        # Only write something to the wiggle files if the number of found start
        # sites is larger than 0.
        if forwardDepth :
            forwardHandle.write("%s\t%i\n" % (position, forwardDepth))
        if reverseDepth :
            reverseHandle.write("%s\t%i\n" % (position, reverseDepth))
        line = pileupHandle.readline()
    #while
#makeWiggles

def main() :
    """
    Do argument checking and call the makeWiggles() function.
    """

    parser = argparse.ArgumentParser(
        prog = 'sageWiggle',
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = '',
        epilog = """""")

    parser.add_argument('-i', dest = 'input', type = argparse.FileType('r'),
        default = sys.stdin, help = 'Pileup file.')
    parser.add_argument('-o', dest = 'output', type = argparse.FileType('w'),
        required = True, nargs = 2,
        help = 'The forward and reverse wiggle files.')
    parser.add_argument('-p', dest = 'threePrime', action = "store_true",
        default = False, help = "Record the 3' end of the reads.")

    arguments = parser.parse_args()

    makeWiggles(arguments.input, arguments.output[0], arguments.output[1],
        os.path.splitext(os.path.basename(arguments.input.name))[0],
        arguments.threePrime)
#main

if __name__ == "__main__" :
    main()
