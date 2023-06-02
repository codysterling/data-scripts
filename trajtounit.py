#!/usr/bin/env python3

import numpy as np
import sys
import argparse
from collections import Counter

#####
# Program introduction

p = argparse.ArgumentParser(description="Converts a trajectory .xyz file (using a corresponding cell file) to an averaged unit cell; only for orthogonal cells")
p.add_argument('-xyz',nargs=1,type=str,help="The XYZ file with the trajectory",metavar='',required=True)
p.add_argument('-abc',nargs=1,type=str,help="The cell file with the ABC values",metavar='',required=True)
p.add_argument("-s",nargs=3,help="The size of the supercell in ABC",metavar=('A','B','C'),required=True)
p.add_argument("--start",nargs=1,type=int,default=[0],help="The step number at which to start the averaging",metavar='N')
args = p.parse_args ()

xyzfile = args.xyz[0]
cellfile = args.abc[0]
supers = args.s
for i in range(0,len(supers)):
  supers[i] = float(supers[i])
min = int(args.start[0])

#####
# Defining functions

def combindex(carray,indices):
  isum = np.zeros([3])
  for i in indices:
    isum += carray[i]
  isum /= len(indices)
  return isum

#####
# Start program

with open(xyzfile) as fp:
  xyz = fp.readlines()

with open(cellfile) as fp:
  cell = fp.readlines()

# Getting number of atoms for the final position array
natoms = int(xyz[0].split()[0])
fullpos = np.zeros([natoms,3])

# Making atom names array
atoms = np.array([])
for i in range(2,natoms+2):
  atoms = np.append(atoms,[xyz[i].split()[0]],axis=0)

lnum = 0
count = 0
while lnum < len(xyz):
  natoms = int(xyz[lnum].split()[0])
  inum = int(xyz[lnum+1].split()[2][:-1])
  if inum >= min: # only doing this if we're above the starting step
    # Getting ABC from cell file to match the inum
    for line in cell:
      words = line.split()
      if words[0] == str(inum):
        a = float(words[2])/float(supers[0])
        b = float(words[6])/float(supers[1])
        c = float(words[10])/float(supers[2])
    # Getting the XYZ for this step into array so we can find atomic distances
    for i in range(0,natoms):
      scoords = xyz[lnum+2+i].split()[1:]
      scoords[0] = float(scoords[0])/a
      scoords[1] = float(scoords[1])/b
      scoords[2] = float(scoords[2])/c
      # Updating fullpos with running average
      fullpos[i,:] = (fullpos[i,:]*count + scoords)/(count+1)
    count += 1
  lnum += natoms+2

print("finished the positions, fullpos is {} (size: ({},{})".format(fullpos,fullpos.shape,type(fullpos)))

# First doing some approximate cell wrapping
for i in range(0,natoms):
  for j in range(0,3):
    if fullpos[i,j] > 0: fullpos[i,j] -= np.floor(fullpos[i,j])
    if fullpos[i,j] > 0.85: fullpos[i,j] -= np.round(fullpos[i,j])
print("wrapped the positions, fullpos is {} (size: ({},{})".format(fullpos,fullpos.shape,type(fullpos)))
print("comb 1 and 3 is {}".format(combindex(fullpos,[0,2])))

# Figuring out which atoms are equivalent in the cell
nunits = int(supers[0]*supers[1]*supers[2]) # Number of unit cells in the super cell
nunitatoms = int(len(atoms)/nunits) # Number of atoms in the unit cell
unitatoms = np.array([])
print("unitatoms: ({},{})".format(unitatoms.shape,type(unitatoms)))
unitpos = np.zeros([nunitatoms,3])
cdict = Counter(atoms) # Getting atom counts for each type
atomnums = {}
tsum = 0
for i in set(atoms):
  atomnums[i] = 0
  ntype = int(cdict[i]/nunits)
  print("atom {} has ntype {}, total {}".format(i,ntype,tsum))
  for j in range(0,ntype):
    unitatoms = np.append(unitatoms,[i],axis=0)
  tsum += ntype
print("unitatoms is {}, {}".format(unitatoms,len(unitatoms)))
for i in atomnums.keys():
  print("atom {} has {} atomnums".format(i,atomnums.get(i)))
  print("in unitatoms, atom {} starts at index {}".format(i,np.where(unitatoms == i)[0][0]))
print(nunits,nunitatoms,len(set(atoms)))
print("cdict is {}".format(cdict))
print("num units is {}".format(nunits))

# Going through all atoms and putting them into the unit positions
for index,coords in enumerate(fullpos):
  atype = atoms[index]
  print("index is {}, coords are {}: {}".format(index,atype,coords))

# Wrapping positions into single ABC cell
wpos = fullpos[0,:]
for i in range(0,natoms):
  for j in range(0,3):
    fullpos[i,:] -= wpos
    if fullpos[i,j] > 0: fullpos[i,j] -= np.floor(fullpos[i,j])
    fullpos[i,j] -= np.floor(fullpos[i,j])
fullpos[0,:] = 0
print("wrapped the positions, fullpos is {} (size: ({},{})".format(fullpos,fullpos.shape,type(fullpos)))
print("atoms are {}".format(atoms))
f = open('output.xyz','w')
f.write("{}\n".format(natoms))
f.write("averaged positions for the whole trajectory in abc coordinates\n")
for i in range(0,natoms):
  f.write("{} {} {} {}\n".format(atoms[i],fullpos[i,0],fullpos[i,1],fullpos[i,2]))
f.close()
