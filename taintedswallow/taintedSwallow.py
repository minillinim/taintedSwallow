#!/usr/bin/env python
###############################################################################
#                                                                             #
#    taintedSwallow.py                                                        #
#                                                                             #
#    Utilities for parsing the nucmer outputs and gettin' the job done        #
#                                                                             #
#    Copyright (C) Michael Imelfort                                           #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2013"
__credits__ = ["Michael Imelfort"]
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Dev"

###############################################################################
###############################################################################
###############################################################################
###############################################################################


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class TSUtilityClass():
    """Utilities wrapper"""
    def __init__(self):
        self.contigs = {}
    
    def marryContigs(self,
                     contigFileName1,
                     contigFileName2,
                     nucmerFileName):
        """Marry two sets of contigs to eachother
        
        High level workflow wrapper"""
        self.parseContigs(contigFileName1)
        self.parseContigs(contigFileName2)
        self.parseCoords(nucmerFileName)

    def parseContigs(self, contigFileName):
        """Parse the contigs
        
        contigFileNames should have exactly two entries
        """
        CP = ContigParser()
        with open(contigFileName, "r") as c_fh:
            for cid, seq in CP.readFasta(c_fh):
                self.contigs[cid] = (seq, len(seq)) 
                
    def parseCoords(self, nucmerFileName):
        """parse the nucmer file"""    
        try:
            nuc_fh = open(nucmerFileName, "r")
        except: 
            print "Error opening file:", nucmerFileName, exc_info()[0]
            raise    
        NP = NucMerParser()
        for record in NP.readNuc(nuc_fh):
            print record, self.contigs[record[7]][1], self.contigs[record[8]][1]
            
        nuc_fh.close()

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class NucMerParser:
    """Wrapper class for parsing nucmer output"""
    def __init__(self):
        self.prepped = False
    
    def readNuc(self, fp):
        """Read through a nucmer coords file
        
        this is a generator function
        """
        line = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not self.prepped:
                # we still need to strip out the header
                    for l in fp: # search for the first record
                        if l[0] == '=': # next line is good
                            self.prepped = True
                            break
            # file should be prepped now
            for l in fp:
                fields = l.split('|')
                yield ([int(i) for i in fields[0].split()] +
                       [int(i) for i in fields[1].split()] +
                       [int(i) for i in fields[2].split()] +
                       [float(i) for i in fields[3].split()] +
                       fields[4].split())
            break # done!
        
###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ContigParser:
    """Main class for reading in and parsing contigs"""
    def __init__(self): pass

    def readFasta(self, fp): # this is a generator function
        header = None
        seq = None
        while True:
            for l in fp:
                if l[0] == '>': # fasta header line
                    if header is not None:
                        # we have reached a new sequence
                        yield header, "".join(seq)
                    header = l.rstrip()[1:].partition(" ")[0] # save the header we just saw
                    seq = []
                else:
                    seq.append(l.rstrip())
            # anything left in the barrel?
            if header is not None:
                yield header, "".join(seq)
            break
                               
###############################################################################
###############################################################################
###############################################################################
###############################################################################
                               