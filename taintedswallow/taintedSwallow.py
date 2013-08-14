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

import numpy as np

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class TSUtilityClass():
    """Utilities wrapper"""
    def __init__(self):
        self.contigs = {}
        
        # use token to represent string based IDs
        self.IDString2Int = {"NONE":0}
        self.IDInt2String = {0:"NONE"}
        self.IDCounter = 0
        
        # binz
        self.contigsProcessed = False
        self.bin1Assignments = {}           # contigs assigned to each bin
        self.bin2Assignments = {}
        self.bin2BinLinkCounts = {}         # number of links shared between pairwise bins
        self.perBinLinkCounts = {}          # number of links in total for each bin
    
    def getIntID(self,
                 stringID):
        """Get an integer-based ID given a string-based ID"""
        try:
            # seen before, second time try to store
            intID = self.IDString2Int[stringID]
            return intID
        except KeyError:
            # first time seen 
            self.IDCounter += 1
            # store in the hashes
            self.IDString2Int[stringID] = self.IDCounter
            self.IDInt2String[self.IDCounter] = stringID
            return self.IDCounter
    
    def getStringID(self,
                    intID):
        """Get an string-based ID given a integer-based ID"""
        try:
            # seen before, second time try to store
            return self.IDInt2String[intID]
        except KeyError:
            print "ERROR: ID %d not in store" % intID
            raise 
     
    def marryContigs(self,
                     contigFileName1,
                     contigFileName2,
                     nucmerFileName,
                     identity=0.95):
        """Marry two sets of contigs to eachother
        
        High level workflow wrapper"""
        self.parseContigs(contigFileName1, 1)
        self.parseContigs(contigFileName2, 2)
        self.parseCoords(nucmerFileName)
        self.processLinks()
        self.contigsProcessed = True
        return 0

    def marryBins(self, bins1, bins2):
        """Try to marry bins to eachother based on contig links"""
        if not self.contigsProcessed:
            print "ERROR: You need to marry the contigs first"
            return 1
        
        # parse the bins files!
        self.parsebins(bins1, 1)
        self.parsebins(bins2, 2)
        
        # work out the total amount of links available to each bin
        # and the inter bin links
        # for each bin in group 1
        for int_bid in self.bin1Assignments:
            self.perBinLinkCounts[int_bid] = 0
            # for each contig in this bin
            for int_cid in self.bin1Assignments[int_bid]:
                self.perBinLinkCounts[int_bid] += len(self.contigs[int_cid][0].links)
                linking_contigs = self.contigs[int_cid][0].links.keys()
                # for each linking contig...
                for link_int_cid in linking_contigs:
                    link_int_bid = self.contigs[link_int_cid][3]
                    key = self.makeLinkKey(link_int_bid, int_bid)
                    try:
                        self.bin2BinLinkCounts[key] += 1
                    except KeyError:
                        self.bin2BinLinkCounts[key] = 1
        # and then group 2 
        for int_bid in self.bin2Assignments:
            self.perBinLinkCounts[int_bid] = 0
            # for each contig in this bin
            for int_cid in self.bin2Assignments[int_bid]:
                self.perBinLinkCounts[int_bid] += len(self.contigs[int_cid][0].links)
                linking_contigs = self.contigs[int_cid][0].links.keys()
                for link_int_cid in linking_contigs:
                    link_int_bid = self.contigs[link_int_cid][3]
                    key = self.makeLinkKey(link_int_bid, int_bid)
                    try:
                        self.bin2BinLinkCounts[key] += 1
                    except KeyError:
                        self.bin2BinLinkCounts[key] = 1

        self.buildBinGraph()

    def parsebins(self, binsFileName, groupNumber):
        """parse bin assignments"""
        if groupNumber == 1:
            ba = self.bin1Assignments
        else:
            ba = self.bin2Assignments
            
        with open(binsFileName, "r") as b_fh:
            for line in b_fh:
                [cid, bid] = line.rstrip().split('\t')
                int_cid = self.getIntID(cid)
                int_bid = self.getIntID(bid)
                try:
                    self.contigs[int_cid][3] = int_bid
                    try:
                        ba[int_bid].append(int_cid)
                    except KeyError:
                        ba[int_bid] = [int_cid]
                except KeyError:
                    print "ERROR: contig %s assigned to bin %s but not in contigs file!"
                    return
                
    def parseContigs(self, contigFileName, groupNumber):
        """Parse the contigs
        
        contigFileNames should have exactly two entries
        """
        CP = ContigParser()
        with open(contigFileName, "r") as c_fh:
            for cid, seq in CP.readFasta(c_fh):
                ID = self.getIntID(cid)
                self.contigs[ID] = [ContigLinker(ID, len(seq)), len(seq), groupNumber, 0]
                
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
            self.contigs[intID][0].processLinks(nameDict=self.IDInt2String)
        
    def makeLinkKey(self, c1, c2):
        """make a unique number from two ints"""
        if c1 < c2:
            return c2 * 10000000 + c1
        else:
            return c1 * 10000000 + c2

    def splitLinkKey(self, key):
        k1 = int(key/10000000)
        k2 = key - k1 * 10000000 
        return (k1, k2)

    def buildBinGraph(self):
        """Build a graph showing inter bin relationships"""
        max_perc = 0
        perc_save = {}
        for b_key in self.bin2BinLinkCounts:
            (k1, k2) = self.splitLinkKey(b_key)
            if k1 != 0 and k2 != 0:
                totes = (self.perBinLinkCounts[k1] + self.perBinLinkCounts[k2]) / 2.
                perc = float(self.bin2BinLinkCounts[b_key]) / float(totes)
                perc_save[b_key] = perc
                if perc > max_perc:
                    max_perc = perc

        reducto = 3./ max_perc
        
        print "graph taintedSwallowBinLinks {"
        print "    rankdir=LR;"
        print "    ranksep=10;"
        print "    subgraph bg_1 {"
        print "        rank=\"same\";"
        print "        node [shape=box,fixedsize=true,width=2];",  
        for int_bid in self.bin1Assignments.keys():
            print "        %s;" % self.getStringID(int_bid),
        print   
        print "    }"

        print "    subgraph bg_2 {"
        print "        rank=\"same\";"
        print "        node [shape=box,fixedsize=true,width=2];",  
        for int_bid in self.bin2Assignments.keys():
            # add the prefix "GroopM_" here
            print "        GroopM_%s;" % self.getStringID(int_bid),
        print   
        print "    }"
        
        for b_key in self.bin2BinLinkCounts:
            (k1, k2) = self.splitLinkKey(b_key)
            if k1 != 0 and k2 != 0:
                # add the prefix "GroopM_" here too
                s1 = self.getStringID(k1)
                s2 = self.getStringID(k2)
                if s1[0] == 'C':
                    s2 = "GroopM_" + s2
                else:
                    s1 = "GroopM_" + s1
                print "    %s -- %s [penwidth=%0.4f];" % (s1,
                                                          s2,
                                                          perc_save[b_key] * reducto)
        print "};"
        
        
    def buildLinkGraph(self):
        """Run through the links and build a graphviz-type graph of contig associations"""
        print "graph taintedSwallow {" 
        # first print the node shapes
        group1_nodes = []
        group2_nodes = []
        seen_links = {}
        for intID in self.contigs:
            if self.contigs[intID][2] == 1:
                group1_nodes.append(intID)
            else:
                group2_nodes.append(intID)
        print "node [shape=box,fixedsize=true,width=0.9];",  
        for node in group1_nodes:
            print "%d;" % node,  
        print "node [shape=circle,fixedsize=true,width=0.9];",
        for node in group2_nodes:
            print "%d;" % node,  
        for intID in self.contigs:
            for link in self.contigs[intID][0].links.keys():
                try:
                    seen_links[self.makeLinkKey(intID, link)] += 1
                except KeyError: 
                    seen_links[self.makeLinkKey(intID, link)] = 1
                    print "    %d -- %d;" % (intID, link)
        print "};"
        
    def printLinks(self):
        """printing!"""
        for intID in self.contigs:
            self.contigs[intID][0].printLinks(nameDict=self.IDInt2String)

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

        # force each match to sit at only one position
        # keep the longest one...
        for intID in self.links:
            current_link = self.links[intID][0]
            l_len = 0
            best_link = current_link 
            while(current_link is not None):
                if len(current_link) > l_len:
                    l_len = len(current_link)
                    best_link =  current_link
                current_link = current_link.nextLink
            best_link.nextLink = None
            self.links[intID][0] = best_link
    
        # each link should now only have one contig
        # sort the linking contigs in order of start point
        # go through all the links for this guy and work out which ones are
        # the best representative for each region in the host
        ordered_IDs = np.array(self.links.keys())
        num_links = len(ordered_IDs)
        dead_uns = []
        if len(ordered_IDs) > 1:
            starts = [self.links[x][0].thisStart for x in ordered_IDs]
            sorted_starts = np.argsort(starts)
            # no point sorting single lists
            ignore_next = False 
            for link_count in np.arange(num_links -1):
                if not ignore_next:
                    link_1 = self.links[ordered_IDs[sorted_starts][link_count]][0]
                    link_2 = self.links[ordered_IDs[sorted_starts][link_count+1]][0]
                    overlap = link_1.thisEnd - link_2.thisStart 
                    if overlap > 100:
                        # there is a problem!
                        if len(link_2) < 2000:
                            # link_2 gotta go!
                            dead_uns.append(ordered_IDs[sorted_starts][link_count+1])
                            ignore_next = True
                else:
                    ignore_next = False

        # delete any useless links
        for intID in dead_uns:
            del self.links[intID]

    def printLinks(self, nameDict=None):
        """printing!"""
        print "==============================================================================="        
        print "==============================================================================="        
        if nameDict is None:
            print "BASE: %d (%d bp)" % (self.ID, self.contigLength)
        else:
            print "BASE: %s (%d bp)" % (nameDict[self.ID], self.contigLength)

        ordered_IDs = np.array(self.links.keys())
        num_links = len(ordered_IDs)
        starts = [self.links[x][0].thisStart for x in ordered_IDs]
        sorted_starts = np.argsort(starts)
        for link_count in np.arange(num_links -1):
            intID = ordered_IDs[sorted_starts][link_count]
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

    def __len__(self):
        """length of overlap on this"""
        return self.thisEnd - self.thisStart + 1
        
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
                               