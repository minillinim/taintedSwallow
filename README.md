# taintedSwallow

## Overview

Program which tries to marry two sets of contigs together. Like a reciprocal best blast but using nucmer instead.

## Installation

Should be as simple as

    pip install taintedSwallow

## Example usage

Still under dev, but if you run nucmer on a set of contigs and get coords ussing something like this:

    nucmer <fasta1> <fasta2> --mum --coords
    
This will produce a file called out.coords. Then you run TS like this:

    taintedSwallow <fasta1> <fasta2> out.coords

What will happen? This is still a surprise!

## Help

If you experience any problems using taintedSwallow, open an [issue](https://github.com/minillinim/taintedSwallow/issues) on GitHub and tell us about it.

## Licence and referencing

Project home page, info on the source tree, documentation, issues and how to contribute, see http://github.com/minillinim/taintedSwallow

This software is currently unpublished

## Copyright

Copyright (c) 2013 Michael Imelfort. See LICENSE.txt for further details.
