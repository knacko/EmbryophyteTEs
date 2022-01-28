#!/usr/bin/env python

import sys

Usage ="""
FormatBLASTFasta.py
Converts a RM output TEs into an acceptable format for input into RStudio.

Usage:
FormatRMFasta.py sequences.fasta"""

if len(sys.argv) < 2:
   print Usage

else:
   InFileName = sys.argv[1]
   InFile = open(InFileName,'r')

   seq=""

   for line in InFile:
      if line[:1] == '>':
          if len(seq) > 0:  
              print seq
              seq=""
          elements = line.split('(')
	  name = elements[0].strip()
          print name
      else:
          seq+=line.rstrip()

