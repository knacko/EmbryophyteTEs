#!/usr/bin/env python

import sys

if len(sys.argv) < 2:
   print Usage

else:
   InFileName = sys.argv[1]
   InFile = open(InFileName,'r')

   for line in InFile:
      line = line.replace('   ','\t')
      elements = line.split('\t')

      if len(elements) < 3:
         continue;

      else:
         seq = elements[0]
         name = elements[1]
         family = elements[2]
         print ">%s#%s" % (name, family)
         print seq
