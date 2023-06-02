#!/usr/bin/env python3

import numpy as np
import argparse
import re

p = argparse.ArgumentParser(description="(Only for cells with all 90 degree sides) Takes supercell parameters and returns unit cell parameters over time for a .cell file.  If given proper unitcell parameters it instead returns the ratio.")
p.add_argument('cell',nargs=1,help="The .cell file to be analysed.",metavar='FILE.cell')
p.add_argument('-abc',nargs=3,help="The AxBxC of the unit cell into the super cell",required=True,metavar=('A','B','C'))
p.add_argument('-unit',nargs=3,help="The abc side lengths of the original unit cell (for making ratios)",metavar=('a','b','c'))
args = p.parse_args()

with open(args.cell[0]) as fp:
  for line in fp:
    words = line.split()
    if re.match("\d",words[0]):
      a = float(words[2])/int(args.abc[0])
      b = float(words[6])/int(args.abc[1])
      c = float(words[10])/int(args.abc[2])
      if args.unit:
        a /= float(args.unit[0])
        b /= float(args.unit[1])
        c /= float(args.unit[2])
      print("Step {}: a = {}, b = {}, c = {}".format(words[0],a,b,c))
