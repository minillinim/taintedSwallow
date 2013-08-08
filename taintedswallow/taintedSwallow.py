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
        self.contigString2Int = {}
        self.contigInt2String = {}
        self.contigIDCounter = 0
    
    def getIntID(self,
                 contigID):
        """Get an integer-based contigID given a string-based contigID"""
        try:
            # seen before, second time try to store
            intID = self.contigString2Int[contigID]
            return intID
        except KeyError:
            # first time seen 
            self.contigIDCounter += 1
            # store in the hashes
            self.contigString2Int[contigID] = self.contigIDCounter
            self.contigInt2String[self.contigIDCounter] = contigID
            return self.contigIDCounter
    
    def getStringID(self,
                    contigID):
        """Get an integer-based contigID given a string-based contigID"""
        try:
            # seen before, second time try to store
            return self.contigInt2String[contigID]
        except KeyError:
            print "ERROR: ID %d not in store" % contigID
            raise 
     
    def marryContigs(self,
                     contigFileName1,
                     contigFileName2,
                     nucmerFileName,
                     identity=0.95):
        """Marry two sets of contigs to eachother
        
        High level workflow wrapper"""
        self.parseContigs(contigFileName1)
        self.parseContigs(contigFileName2)
        self.parseCoords(nucmerFileName)
        self.processLinks()
        self.printLinks()

    def parseContigs(self, contigFileName):
        """Parse the contigs
        
        contigFileNames should have exactly two entries
        """
        CP = ContigParser()
        with open(contigFileName, "r") as c_fh:
            for cid, seq in CP.readFasta(c_fh):
                ID = self.getIntID(cid)
                self.contigs[ID] = (ContigLinker(ID, len(seq)), len(seq))  
                
    def parseCoords(self, nucmerFileName):
        """parse the nucmer file"""    
        try:
            nuc_fh = open(nucmerFileName, "r")
        except: 
            print "Error opening file:", nucmerFileName, exc_info()[0]
            raise    
        NP = NucMerParser()
        for record in NP.readNuc(nuc_fh):
            c1_int_id = self.getIntID(record[NP._ID_1])
            c2_int_id = self.getIntID(record[NP._ID_2])
            self.contigs[c1_int_id][0].addLink(c2_int_id,
                                               record[NP._START_1],
                                               record[NP._END_1],
                                               record[NP._START_2],
                                               record[NP._END_2],
                                               self.contigs[c2_int_id][1],
                                               record[NP._IDENTITY]
                                               )

            self.contigs[c2_int_id][0].addLink(c1_int_id,
                                               record[NP._START_2],
                                               record[NP._END_2],
                                               record[NP._START_1],
                                               record[NP._END_1],
                                               self.contigs[c1_int_id][1],
                                               record[NP._IDENTITY]
                                               )
            
        nuc_fh.close()

    def processLinks(self):
        """process all the links we've found"""
        for intID in self.contigs:
            self.contigs[intID][0].processLinks(nameDict=self.contigInt2String)
        
    def printLinks(self):
        """printing!"""
        for intID in self.contigs:
            self.contigs[intID][0].printLinks(nameDict=self.contigInt2String)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ContigLinker:
    """Store information about how contigs are linked"""
    # constants used for managing info within links
    def __init__(self,
                 ID,
                 length):
        self.ID = ID
        self.contigLength = length
        self.links = {}
        self.stretch=400
        self.minKeep=400
    
    def addLink(self,
                intID,
                thisStart,
                thisEnd,
                thatStart,
                thatEnd,
                thatLength,
                identity):
        """add a link between two contigs"""
        # start before end!
        if thisEnd < thisStart:
            tmp = thisStart
            thisStart = thisEnd
            thisEnd = tmp
            tmp = thatStart
            thatStart = thatEnd
            thatEnd = tmp
        
        try:
            chain = self.links[intID][0]
            self.links[intID][0] = self.appendLink(chain,
                                                thisStart,
                                                thisEnd,
                                                thatStart,
                                                thatEnd,
                                                identity) 
        except KeyError:
            # add the first link
            # (thisStart, thisEnd, thatStart, thatEnd, identity, nextLink)
            self.links[intID] = [contigLink(thisStart, thisEnd, thatStart, thatEnd, identity), thatLength]

    def appendLink(self,
                   chain,
                   thisStart,
                   thisEnd,
                   thatStart,
                   thatEnd,
                   identity):
        """append a link to an already existing chain
        
        takes care of link sorting
        """
        # first we check to see if the new link lies before 
        # the current head of the chain
        current_link = chain
        next_link = current_link.nextLink
        if thisStart <= current_link.thisStart:
            # new link lies before the beginning of the chain
            # make a new link at the start and return it to replace the chain
            return contigLink(thisStart, thisEnd, thatStart, thatEnd, identity, nextLink=chain)
        while True:
            if next_link is not None:
                if thisStart <= next_link.thisStart:
                    # insert between current and next
                    current_link.nextLink = contigLink(thisStart,
                                                       thisEnd,
                                                       thatStart,
                                                       thatEnd,
                                                       identity,
                                                       nextLink=next_link)
                    return chain
                else:
                    # lies after the start of the next. Hop along one...
                    current_link = next_link
                    next_link = current_link.nextLink
            else:
                # append to the end of the list
                current_link.nextLink = contigLink(thisStart,
                                                   thisEnd,
                                                   thatStart,
                                                   thatEnd,
                                                   identity)
                return chain

    def processLinks(self, nameDict=None):
        """Keep only those links which make sense"""
        dead_uns = []
        for intID in self.links:
            # go through this guy and work out
            # if we can compress any links here
            current_link = self.links[intID][0] 
            next_link = current_link.nextLink
            skipped_one = False
            while(next_link is not None):
                merged_this_round = False
                if abs(next_link.thisStart - current_link.thisEnd) <= self.stretch:
                    # perhaps we can overlap and join these two fellas
                    if abs(next_link.thatStart - current_link.thatEnd) <= self.stretch:
                        # we can collapse these two
                        curr_len = abs(current_link.thisEnd - current_link.thisStart)+1
                        next_len = abs(next_link.thisEnd - next_link.thisStart)+1
                        new_identity = (curr_len*current_link.identity + next_len*next_link.identity) / (curr_len + next_len + 2) 
                        current_link.thisEnd = next_link.thisEnd
                        current_link.thatEnd = next_link.thatEnd
                        current_link.identity = new_identity
                        current_link.nextLink = next_link.nextLink
                        merged_this_round = True
                if not merged_this_round:
                    if skipped_one:
                        # failed twice, leave next the same and shuffle up current
                        current_link = current_link.nextLink
                        skipped_one = False
                    else:
                        # keep current and skip one ahead
                        next_link = next_link.nextLink
                        skipped_one = True
                else:
                    # merged, so current stays put and we need to look at who he (now) points to
                    next_link = current_link.nextLink

            # now go through and kill all the short links
            last_link = None
            current_link = self.links[intID][0] 
            while(current_link is not None):
                if abs(current_link.thisStart - current_link.thisEnd) < self.minKeep:
                    # too short
                    if last_link is None:
                        # at the start of the list
                        current_link = current_link.nextLink
                        self.links[intID][0] = current_link
                    else:
                        # skip over this one
                        current_link = current_link.nextLink
                        last_link.nextLink = current_link
                else:
                    # OK
                    last_link = current_link
                    current_link = current_link.nextLink
            
            if self.links[intID][0] is None:
                dead_uns.append(intID)

        # delete any useless links
        for intID in dead_uns:
            del self.links[intID]
        
        # now go through all the links for this guy and work out which ones are
        # the best representative for each region in the host
        current_link = self.links[intID][0] 
        while(current_link is not None):
            current_link = current_link.nextLink
        

    def printLinks(self, nameDict=None):
        """printing!"""
        print "==============================================================================="        
        print "==============================================================================="        
        if nameDict is None:
            print "BASE: %d (%d bp)" % (self.ID, self.contigLength)
        else:
            print "BASE: %s (%d bp)" % (nameDict[self.ID], self.contigLength)
        
        for intID in self.links:
            print "..............................................................................."        
            if nameDict is None:
                print "LINK: %d (%d bp)" % (intID, self.links[intID][1])
            else:
                print "LINK: %s (%d bp)" % (nameDict[intID], self.links[intID][1])
            current_link = self.links[intID][0]
            while(current_link is not None):
                print current_link
                current_link = current_link.nextLink
            
###############################################################################
###############################################################################
###############################################################################
###############################################################################
            
class contigLink:
    """Generic linked list"""
    def __init__(self, thisStart, thisEnd, thatStart, thatEnd, identity, nextLink=None):
        self.thisStart = thisStart
        self.thisEnd = thisEnd
        self.thatStart = thatStart
        self.thatEnd = thatEnd
        self.identity = identity
        self.nextLink = nextLink
        
    def __str__(self):
        """print statement"""
        return "%d -> %d (%d) | %d -> %d (%d) identity: %f" % (self.thisStart,
                                                               self.thisEnd,
                                                               (self.thisEnd-self.thisStart)+1,
                                                               self.thatStart,
                                                               self.thatEnd,
                                                               (self.thatEnd-self.thatStart)+1,
                                                               self.identity)
                   
###############################################################################
###############################################################################
###############################################################################
###############################################################################

class NucMerParser:
    """Wrapper class for parsing nucmer output"""
    # constants to make the code more readable
    _START_1  = 0
    _END_1    = 1
    _START_2  = 2
    _END_2    = 3
    _LEN_1    = 4
    _LEN_2    = 5
    _IDENTITY = 6
    _ID_1     = 7
    _ID_2     = 8

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
                               