#!/usr/bin/env python3

import numpy as np
import sys
import glob
import matplotlib.pyplot as plt
import os
import argparse
from datetime import datetime as dt
import random
np.set_printoptions(threshold=sys.maxsize)

#####
# Initializing program

p = argparse.ArgumentParser(description="Correlates some geometric metric in an XYZ file with the associated spectra.\nHere's a new line!\n\nSkipped twice!",formatter_class=argparse.RawTextHelpFormatter)
p.add_argument('-xyz',nargs=1,help="The XYZ file which generated the spectra, read for the geometric analysis",metavar="XYZ.xyz",required=True)
p.add_argument('-abc',nargs=3,type=float,help="The ABC dimensions of the cell",metavar=('A','B','C'),required=True)
p.add_argument('-s','--spec',nargs=1,help="The atom type of the spectra",metavar='I',required=True)
p.add_argument('-a1','--atom1',nargs=1,help="The first atom type (closest to -s atom)",metavar='A1',required=True)
p.add_argument('-n','--natom1',nargs=1,type=int,default=[1],help="The number of A1 (ordered closest to -s atom) to look at (default: 1)",metavar='N')
p.add_argument('-a2','--atom2',nargs=1,help="The second atom type (closest to A1 atom)",metavar='A2',required=True)
p.add_argument('-t','--type',nargs=1,default=['min'],help="When -n > 1, set -t to min or sum to determine how to parse for the spectra ordering (default: min)\nHere's a new line for type")
p.add_argument('-avg','--avg',nargs=1,type=int,help="Enables block averaging of the spectra into N blocks",metavar='N')
p.add_argument('--noscale',action='store_false',help="Turns off scaling of intensities for plotting (enable for comparing separate data)")
p.add_argument('--noprint',action='store_false',help="Turns off printing of GNU and Py files")
p.add_argument('--noplot',action='store_false',help="Turns off plotting of the spectra")
args = p.parse_args()

# Parsing arguments
geom = args.xyz[0]
a = args.abc[0]
b = args.abc[1]
c = args.abc[2]
spec = args.spec[0]
atom1 = args.atom1[0]
natom1 = args.natom1[0]
atom2 = args.atom2[0]
typer = args.type[0]

#####
# Defining functions

def collapse(vec1,a,b,c):
  vec1[0] = vec1[0] - int(round(vec1[0]/a))*a
  vec1[1] = vec1[1] - int(round(vec1[1]/b))*b
  vec1[2] = vec1[2] - int(round(vec1[2]/c))*c
  return vec1

def nlines(file):
  count = 0
  with open(file) as fp:
    for line in fp:
      count+=1
  return count

#####
# Running program

# Catching some mistakes
if natom1 != 1 and atom1 != 'H':
  print("natom1 is {} but atom1 is {}, this might break stuff!".format(natom1,atom1))
if spec == atom1 and natom1 != 1:
  print("spec is {} and atom1 is {} with natom1 {}, resetting natom1 to 1.".format(spec,atom1,natom1))
  natom1 = 1

# Slurping up the coords and atoms
with open(geom) as fp:
  natoms = int(fp.readline())

# Initializing arrays
atoms = np.array([])
coords = np.empty([0,3])

with open(geom) as fp:
  # Skipping first two lines of xyz file to get to coords
  fp.readline()
  fp.readline()
  # Going through all atoms to make atom name and coordinates array
  for line in fp:
    words = line.split()
    atoms = np.append(atoms,[words[0]],axis=0)
    coords = np.append(coords,[[float(words[1]),float(words[2]),float(words[3])]],axis=0)

# Going through atomlist to generate atom numbers for spec atom
speclist = []
for index,atom in enumerate(atoms):
  if atom == spec:
    speclist.append(index+1)

# Making empty dict that will be filled with criteria later
distdict = {}

## This is the main loop, ordering the spectra by some coordinate
for num in speclist:
  speccoords = coords[num-1]
  # Generating the list for atom1
  a1dist = {}
  for index,atom in enumerate(atoms):
    if atom == atom1:
      atom1coords = coords[index]
      a1dist[np.linalg.norm(collapse(speccoords-atom1coords,a,b,c))] = index
  near = sorted(a1dist.keys())[0:natom1]
  nearnum = [a1dist.get(d) for d in near]
  a2dist = []
  # Going through atom1 list to get atom2 distances
  for a1 in nearnum:
#    print("a1 is {}".format(a1))
    atom1coords = coords[a1]
    dist = 100
    for index,atom in enumerate(atoms):
      if atom == atom2:
        if index != a1:
          atom2coords = coords[index]
          cdist = np.linalg.norm(collapse(atom1coords-atom2coords,a,b,c))
          if cdist < dist:
            dist = cdist
            a2 = index+1
#    print("  For a1: {} a2 is {} with dist {}".format(a1+1,a2,dist))
    a2dist.append(dist)
  # Rough check of a2dist
  if len(a2dist) != natom1:
    print("Something is wrong with a2dist: {}, length {} should be {}!".format(a2dist,len(a2dist),natom1))
    exit()
  # Processing a2dist according to -t type for the distdict
  if typer == 'min':
    nadd = min(a2dist)
  elif typer == 'sum':
    nadd = sum(a2dist)
  elif typer == 'diff':
    nadd = max(a2dist)-min(a2dist)
  elif typer == 'diff2':
    nadd = min(a2dist) + np.median(a2dist) - max(a2dist)
  elif typer == 'diff3':
    nadd = max(a2dist) + min(a2dist) - 2*np.median(a2dist)
  elif typer == 'diff4':
    nadd = abs(max(a2dist) + min(a2dist) - 2*np.median(a2dist))
  elif typer == 'avg':
    nadd = sum(a2dist)/len(a2dist)
  else:
    print("Type {} not recognized!".format(typer))
    exit()
  distdict[nadd] = num

# Rough check of distdict
if len(distdict) != len(speclist):
  print("Distdict ({}) {} doesn't match speclist length of {}, code is broken!".format(len(distdict),distdict,len(speclist)))
  print(distdict.values())
  exit()

# Now with the dictionary of criteria, we order this and slurp out the spectra info in order
npts = 5000
carray = np.empty([npts,2,len(speclist)])
sortdist = sorted(distdict.keys())

for index,dist in enumerate(sortdist):
  num = distdict.get(dist)
  flist = glob.glob('*xas*at{}*200*.dat'.format(num))
  # Checking flist args, break if not just one match
  if len(flist) > 1:
    print('There are {} files that match num {}!  Exiting.'.format(len(flist),num))
    exit()
  elif len(flist) < 1:
    print('There are no files that match num {}!  Exiting.'.format(num))
    exit()
  file = flist[0]
  lnum = 0
  with open(file) as fp:
    for line in fp:
      carray[lnum,0,index] = line.split()[0]
      carray[lnum,1,index] = float(line.split()[1])
      lnum+=1

# Cleaning up the spectra so they have the same x values for averaging
for i in range(0,len(speclist)):
  x = carray[:,0,i]
  x[x==0] = np.nan
  y = carray[:,1,i]
  y[y==0] = np.nan
  x[y==0] = np.nan
xmin = min(carray[0,0,:])
xmax = max([max(carray[:,0,i]) for i in range(0,len(speclist))])
xrange = np.linspace(xmin,xmax,npts)
for i in range(0,len(speclist)):
  x = [k for k in carray[:,0,i] if str(k) != 'nan']
  y = [j for j in carray[:,1,i] if str(j) != 'nan']
  carray[:,1,i] = np.interp(xrange,x,y)
  carray[:,0,i] = xrange

# Rescaling/shifting intensities so they plot nicely
# Rescaling
maxsep = sortdist[-1] - sortdist[0]
maxint = np.amax(carray[:,1,:])
if args.noscale:
  print("Enabling intensity scaling!")
  carray[:,1,:] = carray[:,1,:] * maxsep/(2*maxint)
else: print("No scaling!")
maxint = np.amax(carray[:,1,:])
# Shifting
if (args.avg and args.avg[0] != 1) or not args.avg:
  for i in range(0,len(speclist)):
    carray[:,1,i] = carray[:,1,i] + sortdist[i]
maxint = np.amax(carray[:,1,:])

# Now doing block averaging of the spectra
if args.avg:
  nblks = args.avg[0]
  narray = np.zeros([npts,2,nblks])
  ymin = min(carray[0,1,:])
  ymax = max(carray[0,1,:])
  yrng = ymax - ymin
  jcnt = np.zeros([nblks])
  # Setting x values for blocks
  # Adding carray spectra into the blocks and averaging
  for i in range(0,len(speclist)):
    for j in range(0,nblks):
      if carray[0,1,i] <= ymin + yrng*(j+1)/nblks:
        narray[:,1,j] = (narray[:,1,j]*jcnt[j] + carray[:,1,i])/(jcnt[j]+1)
        jcnt[j]+=1
        break

#####
# Initializing plot

titlesize=15
axissize=15
ticksize=10
textsize=10

# General plot parameters
ax = plt.subplot()
plt.xlabel('Energy (eV)',size=axissize)
plt.xlim(carray[0,0,0],carray[-1,0,0])
xlen = carray[-1,0,0] - carray[0,0,0]
plt.xticks(size=ticksize)
ax.set_xticks(minor=True,ticks=np.arange(np.ceil(carray[0,0,0]),np.ceil(carray[-1,0,0]),2))
ax.set_xticks(minor=True,ticks=np.arange(np.ceil(carray[0,0,0]),np.ceil(carray[-1,0,0]),1))
plt.tick_params(axis='x',which='major',direction='in',length=8)
plt.tick_params(axis='x',which='minor',direction='in',length=5)
ylen = maxint - np.amin(carray[0,1,:])
plt.ylim(np.amin(carray[0,1,:])-ylen*0.05,maxint+ylen*0.05)

# Plotting main data
for i in range(0,len(speclist)):
  x = carray[:,0,i]
  x[x==0] = np.nan
  y = carray[:,1,i]
  y[y==0] = np.nan
  atomnum = distdict.get(sortdist[i])
  ax.plot(x,y,label='{}{}'.format(spec,atomnum),linewidth=1)

# Adding block averages
if args.avg:
  for i in range(0,nblks):
    narray[:,0,i] = xrange
    ax.plot(narray[:,0,i],narray[:,1,i],color='black',linewidth=3)
    ax.axhline(narray[0,1,i],color='black',linewidth=2,linestyle='--')
    if args.avg[0] != 1:
      plt.text(narray[0,0,i]+xlen*.01,narray[0,1,i]+ylen*.0075,'{:.2f}'.format(narray[0,1,i]),size=textsize)

# Setting plot parameters for overlay average (1 block)
if args.avg and args.avg[0] == 1:
  plt.ylabel('Intensity (arbitrary units)',size=axissize)
  plt.yticks([])
  handles, labels = ax.get_legend_handles_labels()
  ax.legend(handles[::-1], labels[::-1], loc='center right')
  at = ''
else:
# Setting plot parameters for proper separating/blocking
  plt.ylabel('{}-{} distance (Å)'.format(atom1,atom2),size=axissize)
  plt.yticks(size=ticksize)
  plt.tick_params(axis='y',which='major',direction='inout',length=12)
  handles, labels = ax.get_legend_handles_labels()
  ax.legend(handles[::-1], labels[::-1], loc='center right')

  # Adding low/high values for scale
  tshift = xlen*0.025
  plt.text(carray[0,0,0]-tshift,carray[0,1,0],'{:.2f}'.format(sortdist[0]),size=textsize)
  plt.text(carray[0,0,-1]-tshift,carray[0,1,-1],'{:.2f}'.format(sortdist[-1]),size=textsize)

  at = ', Coord: {}-{}'.format(atom1,atom2)
  if natom1 > 1:
    at = at+'\nNumber of {}: {}, comparison: {}'.format(atom1,natom1,typer)
  if args.avg:
    at = at+' Averaged with {} blocks'.format(nblks)

  ax2 = ax.twinx()
  ax2.set_ylabel('Intensity (arbitrary units)',size=axissize)
  ax2.set_yticks([])

cwds = os.getcwd().split('/')
plt.title('{} K-edge XAS\nTPHH, 6-311++G2d2p\n{}{}\n{}'.format(spec,cwds[-2],at,dt.now().strftime("%d %b %Y")))

# Saving array
# For use combining spectra with the Python script:
with open('spectra-py-{}-{}{}-n{}-t{}.dat'.format(spec,atom1,atom2,natom1,typer),'w') as fp:
  for i in range(0,npts):
    fp.write("{} ".format(carray[i,0,0]))
  fp.write("\n")
  for j in range(0,len(speclist)):
    for i in range(0,npts):
      fp.write("{} ".format(carray[i,1,j]))
    fp.write("\n")

# For plotting in gnuplot:
with open('spectra-gnu-{}-{}{}-n{}-t{}.dat'.format(spec,atom1,atom2,natom1,typer),'w') as fp:
  for i in range(0,npts):
    fp.write(str(carray[i,0,0]))
    for j in range(0,len(speclist)):
      fp.write(" {}".format(carray[i,1,j]))
    fp.write("\n")

if args.noplot:
  plt.show()
