# taintedSwallow

## Overview

Program which tries to marry two sets of contigs together. Like a reciprocal best blast but using nucmer instead.

## Installation

Should be as simple as

    pip install taintedSwallow

## Example usage 1

Still under dev, but if you run nucmer on a set of contigs and get coords using something like this:

    nucmer <fasta1> <fasta2> --mum --coords
    
This will produce a file called out.coords. which should look like this:

    <path to file 1> <path to file 2>
    NUCMER

        [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  | [TAGS]
    =====================================================================================
           1      535  |      109      643  |      535      535  |    99.81  | ACDRDN_C00001(id=20317)  contig_10080_63190_63189_0_l=037867_37867
           1      535  |    73282    72748  |      535      535  |    99.81  | ACDRDN_C00001(id=20317)  contig_10354_64637_64636_0_l=073390_73390
           1      535  |      535        1  |      535      535  |    99.81  | ACDRDN_C00001(id=20317)  contig_10529_65339_65338_0_l=041340_41340
    ...

Then you run TS like this:

    taintedSwallow <fasta1> <fasta2> out.coords

TS munges the coords file and produces output which looks like this:

    ===============================================================================
    BASE: FASTA1_CONTIG_ID_1(id=20491) (19293 bp)
    ...............................................................................
    LINK: FASTA_2_CONTIG_ID_1 (13379 bp)
    1 -> 13877 (13877) | 9 -> 13340 (13332) identity: 99.687824
    ===============================================================================
    ===============================================================================
    BASE: FASTA1_CONTIG_ID_2(id=20492) (105677 bp)
    ...............................................................................
    LINK: FASTA_2_CONTIG_ID_2 (1041 bp)
    1 -> 987 (987) | 1005 -> 19 (-985) identity: 100.000000
    ...............................................................................
    LINK: FASTA_2_CONTIG_ID_3 (18919 bp)
    1056 -> 20077 (19022) | 51 -> 18919 (18869) identity: 99.073147
    ...............................................................................
    LINK: FASTA_2_CONTIG_ID_4 (31085 bp)
    19965 -> 32811 (12847) | 12859 -> 14 (-12844) identity: 99.990000
    ...............................................................................
    LINK: FASTA_2_CONTIG_ID_5 (50197 bp)
    33262 -> 83980 (50719) | 56 -> 50185 (50130) identity: 99.598672
    ===============================================================================
    ===============================================================================
    BASE: FASTA1_CONTIG_ID_3(id=20493) (20737 bp)
    ===============================================================================

In this example contig 1 form fasta 1 matches contig 1 from fasta 2; contig 2 from fasta 1 matches 4 contigs in fasta 2; contigs 3 in fasta 1 has no matches.

But overall, this is a boring way to run this code. I mainly intend to use it as an api. Pydoc the code for more deets peeps.

## Example usage 2

You can make tab separated bin assignments files which have this format:

    FASTA_1_CONTIG_ID_1	BIN_X
    FASTA_1_CONTIG_ID_2	BIN_Y
    FASTA_1_CONTIG_ID_3	BIN_Z
    FASTA_1_CONTIG_ID_4	BIN_X
    ...

You can run TS like this:

    taintedSwallow <fasta1> <fasta2> out.coords --bins1 <bins1> --bins2 <bins2>

And then TS will produce graphviz style dot output which maps links between bins. This graph is useful to see how bins are related to eachother.

## Example usage 3


## Help

If you experience any problems using taintedSwallow, open an [issue](https://github.com/minillinim/taintedSwallow/issues) on GitHub and tell us about it.

## Licence and referencing

Project home page, info on the source tree, documentation, issues and how to contribute, see http://github.com/minillinim/taintedSwallow

This software is currently unpublished

## Copyright

Copyright (c) 2013 Michael Imelfort. See LICENSE.txt for further details.
