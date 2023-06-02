#!/usr/bin/env python3

import numpy as np
import sys

file=sys.argv[1]
olist=[]
clist=[]
nlist=[]
hlist=[]

with open(file) as fp:
  for line in fp:
    words = line.split()
    if words[0] == 'C':
      clist.append(line)
    elif words[0] == 'N':
      nlist.append(line)
    elif words[0] == 'H':
      hlist.append(line)
    else:
      olist.append(line)

mlist=olist+clist+nlist+hlist
with open(file[:-4]+'-reordered'+file[-4:],'w') as fp:
  for i in mlist:
    fp.write(i)
