#!/usr/bin/env python3

import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)

def angle(vec1,vec2):
  vec1norm = vec1/np.linalg.norm(vec1)
  vec2norm = vec2/np.linalg.norm(vec2)
  angle = np.arccos(np.clip(np.dot(vec1norm,vec2norm),-1,1))*180/np.pi
  return angle

def distance(atom1,atom2,a,b,c):
  xdist = atom1[0]-atom2[0]
  ydist = atom1[1]-atom2[1]
  zdist = atom1[2]-atom2[2]
  xdist = xdist - int(round(xdist/a))*a
  ydist = ydist - int(round(ydist/b))*b
  zdist = zdist - int(round(zdist/c))*c
  return np.sqrt(xdist**2+ydist**2+zdist**2)

def getkey(dict,val):
  return list(dict.keys())[list(dict.values()).index(val)]

#####
# Start of program

n = 6
a = 18.6625961569
b = 18.7353338603
c = 25.8484418991

file=sys.argv[1]

cdict={}
ndict={}
cndict={}
middict={}
nearmidnum={}
count=1

# Reading file for atoms/coords
with open(file) as fp:
  # Skipping first two lines of xyz file to get to coords
  fp.readline()
  fp.readline()
  # Going through all atoms to make get carbon and nitrogen positions
  for line in fp:
    words=line.split()
    if words[0] == 'C':
      cdict[count] = np.array([float(words[1]),float(words[2]),float(words[3])])
    elif words[0] == 'N':
      ndict[count] = np.array([float(words[1]),float(words[2]),float(words[3])])
    count+=1

# Building C-N dictionary in the form cndict[#N]=#C
for i in ndict.keys():
  dist = 100
  ncoords=ndict.get(i)
  for j in cdict.keys():
    ccoords=cdict.get(j)
    currdist=distance(ncoords,ccoords,a,b,c)
    if currdist < dist:
      dist=currdist
      cnum=j
  cndict[i]=cnum

# Making midpoints dictionary
for i in ndict.keys():
  ncoords=ndict.get(i)
  ccoords=cdict.get(cndict.get(i))
  middict[i]=(ncoords+ccoords)/2

# Now that we have the midpoints, make a list of the n closest ones for each
for i in ndict.keys():
  ddict={}
  point=middict.get(i)
  for j in ndict.keys():
    if i != j:
      ddict[distance(middict.get(i),middict.get(j),a,b,c)] = j
  nearmidnum[i] = [ddict.get(d) for d in sorted(ddict.keys())[0:n]]

# Now for each N we go through its nearest neighbors and get vector angles
for i in sorted(ndict.keys()):
  ncoords=ndict.get(i)
  ccoords=cdict.get(cndict.get(i))
  vec1=ncoords-ccoords
  print('N{}:'.format(i),end=' ')
  for j in nearmidnum.get(i):
    ncoords2=ndict.get(j)
    ccoords2=cdict.get(cndict.get(j))
    vec2=ncoords2-ccoords2
    ang=angle(vec1,vec2)
    print('{:.4f}'.format(ang),end=' ')
  print()
