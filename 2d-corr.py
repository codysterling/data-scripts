#!/usr/bin/env python3

import numpy as np
import sys
import argparse

def distance(atom1,atom2,a,b,c):
  xdist = float(atom1[0])-float(atom2[0])
  ydist = float(atom1[1])-float(atom2[1])
  zdist = float(atom1[2])-float(atom2[2])
  xdist = xdist - int(round(xdist/a))*a
  ydist = ydist - int(round(ydist/b))*b
  zdist = zdist - int(round(zdist/c))*c
  return np.sqrt(xdist**2+ydist**2+zdist**2)

def angle(vec1,vec2):
  vec1norm = vec1/np.linalg.norm(vec1)
  vec2norm = vec2/np.linalg.norm(vec2)
  angle = np.arccos(np.clip(np.dot(vec1norm,vec2norm),-1,1))*180/np.pi
  return angle

#####
# Program introduction

p = argparse.ArgumentParser(description="Reads a trajectory and plots a 2D correlation of two interatomic distance coordinates")
p.add_argument('-xyz',nargs=1,type=str,help="The XYZ file with the trajectory",metavar='',required=True)
p.add_argument('-c',nargs=1,type=str,help="The .cell file with ABC values corresponding to the trajectory",metavar='',required=True)
p.add_argument('-x',nargs=2,default=['H','I'],help="The two atoms defining the coordinate on the X axis")
p.add_argument('-y',nargs=2,default=['N','H'],help="The two atoms defining the coordinate on the Y axis")
p.add_argument('-ang',nargs=1,type=float,help="Only include pairs where the dihedral angle is above cutoff",metavar='D')
args = p.parse_args()

xyzfile = args.xyz[0]
cellfile = args.c[0]
xatoms = args.x
yatoms = args.y

# Checking that there's a common atom between the axes; otherwise we don't know how to order them
catom = list(set(xatoms) & set(yatoms))
if len(catom) == 2:
  print("The lists x: {} and y: {} are the same, the 2D plot would just be a straight line.".format(xatoms,yatoms))
  exit()
elif len(catom) == 0:
  print("The lists x: {} and y: {} have no commonalities, don't know how to order the points.".format(xatoms,yatoms))
  exit()
elif len(catom) != 1:
  print("Something is wrong with the lists x: {} and y: {}, try to fix?".format(xatoms,yatoms))
  exit()

# Getting the atom types to go through
ctype = catom[0]
xtype = xatoms[(xatoms.index(ctype)+1)%2]
ytype = yatoms[(yatoms.index(ctype)+1)%2]

with open(xyzfile) as fp:
  xyz = fp.readlines()

with open(cellfile) as fp:
  cell = fp.readlines()

lnum = 0
if args.ang:
  ang='-ang{:g}'.format(args.ang[0])
else: ang=''
f = open('{}-2dcorr-{}{}-{}{}{}.dat'.format(xyzfile[:-4],args.x[0],args.x[1],args.y[0],args.y[1],ang),'w')
while lnum < len(xyz):
  natoms = int(xyz[lnum].split()[0])
  inum = xyz[lnum+1].split()[2][:-1]
  count = 0
  # Getting ABC from cell file to match the inum
  for line in cell:
    words = line.split()
    if words[0] == inum:
      a = float(words[2])
      b = float(words[6])
      c = float(words[10])
  # Getting the XYZ for this step into array so we can find atomic distances
  tarray = np.empty([natoms,4],dtype=object)
  for i in range(0,natoms):
    tarray[i,:] = np.array(xyz[lnum+2+i].split())
  # Now searching through tarray to get the distances
  for words in tarray:
    if words[0] == ctype: # If we hit the common atom
      cc = np.array(words[1:],dtype=float)
      xdist = 100
      ydist = 100
      for words2 in tarray: # Then look for its partners:
        if words2[0] == xtype: # X partner
          cxdist = distance(words[1:],words2[1:],a,b,c)
          if cxdist < xdist:
            xc = np.array(words2[1:],dtype=float)
            xdist = cxdist
        elif words2[0] == ytype: # Y partner
          cydist = distance(words[1:],words2[1:],a,b,c)
          if cydist < ydist:
            yc = np.array(words2[1:],dtype=float)
            ydist = cydist
      # Now that we have the x and y distances, we print to the file as x, y
      # Only print if angle is above cutoff
      if args.ang:
        if angle(xc-cc,yc-cc) > args.ang[0]:
          if ydist < 1.6:
            f.write("{} {}\n".format(xdist,ydist))
      else:
        if ydist < 1.6:
          f.write("{} {}\n".format(xdist,ydist))
  lnum += natoms+2
