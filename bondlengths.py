#!/usr/bin/env python3

import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)

def distance(atom1,atom2,a,b,c):
  xdist = atom1[0]-atom2[0]
  ydist = atom1[1]-atom2[1]
  zdist = atom1[2]-atom2[2]
  xdist = xdist - int(round(xdist/a))*a
  ydist = ydist - int(round(ydist/b))*b
  zdist = zdist - int(round(zdist/c))*c
  return np.sqrt(xdist**2+ydist**2+zdist**2)

#####
# Explaining code

if len(sys.argv) > 2:
  print('Too many arguments, exiting.')
  exit()
elif len(sys.argv) < 2:
  print('This code reads the .xyz file given as the argument, with ABC values requested or hardcoded in.')
  print('For every atom, it returns the atom number and distance of the nearest neighbor for all atom types in the file.')
  print('For C or N, it also returns the closest H-I or H-Pb for each attached hydrogen.')
  exit()

#####
# Start of program

a = 17.71
b = a
c = 25.318

file=sys.argv[1]

with open(file) as fp:
  natoms = int(fp.readline())

# Initializing arrays
atoms = np.array([])
coords = np.empty([0,3])

with open(file) as fp:
  # Skipping first two lines of xyz file to get to coords
  fp.readline()
  fp.readline()
  # Going through all atoms to make atom name and coordinates array
  for line in fp:
    words=line.split()
    atoms = np.append(atoms,[words[0]],axis=0)
    coords = np.append(coords,[[float(words[1]),float(words[2]),float(words[3])]],axis=0)

# Going through each atom and finding its closest neighbor of all atom types
for i in range(0,len(atoms)):
  print('Atom num {}: {}'.format(i+1,atoms[i]))
  # For each atom type:
  for u in np.unique(atoms):
    dist = 100
    for j in range(0,len(atoms)):
      if i != j and atoms[j] == u:
        currdist = distance(coords[i],coords[j],float(a),float(b),float(c))
        if currdist < dist:
          dist = currdist
          atomnum = j+1
    print('  Closest {} is atom {} at {} Å'.format(u,atomnum,dist))
  if atoms[i] == 'N' or atoms[i] == 'C':
    # Getting all N/C-H distances
    hdist = {}
    for k in range(0,len(atoms)):
      if atoms[k] == 'H':
        hdist[distance(coords[i],coords[k],float(a),float(b),float(c))] = k+1
    nearh = sorted(hdist.keys())[0:3]
    nearhnum = [hdist.get(d) for d in nearh]
    for hatom in nearhnum:
      dist = 100
      for l in range(0,len(atoms)):
        if atoms[l] == 'Pb' or atoms[l] == 'I':
           currdist = distance(coords[hatom-1],coords[l],float(a),float(b),float(c))
           if currdist < dist:
             dist = currdist
             atomnum = l+1
      print('    For H{}, closest is {}{} at {} Å'.format(hatom,atoms[atomnum-1],atomnum,dist))
