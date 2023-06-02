#!/usr/bin/env python3

import numpy as np
import sys
import glob
import matplotlib.pyplot as plt
import os
import argparse
from datetime import datetime as dt

#####
# Initializing program

p = argparse.ArgumentParser(description="Does block averaging over multiple files of spectral data.  Files must be generated from 1ax-corr.py.")
p.add_argument('files',nargs='+',help="List of files to be combined for averaging",metavar="*.dat")
p.add_argument('-n',nargs=1,type=int,help="Number of blocks to average the data",required=True)
p.add_argument('-o',action='store_true',help="Plots block averages on the same level to compare them directly")
p.add_argument('-c',nargs=1,default=[1],help="Cutoff to not plot block averages without at least C spectra in them")
p.add_argument('-s',action='store_true',help="Saves file to the directory in the script")
p.add_argument('-br',nargs=2,type=float,help="Starting and ending values for the fixing block averaging",metavar=('START','END'))
p.add_argument('-ns',nargs=1,type=int,help="Number of spectra per file",required=True)
args = p.parse_args()

# Parsing arguments
nfiles = len(args.files)
n = int(args.n[0])

npts = 5000
nspec = int(args.ns[0])

#####
# Defining functions

# Taken from: https://realpython.com/python-rounding/
def truncate(n,dec=10):
  x = 10**dec
  return int(n*x)/x

#####
# Running program

farray = np.empty([nspec+1,npts,nfiles])
# Reading the files and getting the info
for index,file in enumerate(args.files):
  lcnt = 0
  with open(file) as fp:
    for line in fp:
      words = line.split()
      farray[lcnt,:,index] = words
      lcnt += 1

# Cleaning up the spectra so they have the same x values for averaging
xmin = np.amin(farray[0,0,:])
xmax = np.amax(farray[0,-1,:])
xrange = np.linspace(xmin,xmax,npts)
yarray = np.empty([nspec*nfiles,npts])
for i in range(0,nfiles):
  for j in range(0,nspec):
    yarray[i*nspec+j,:] = np.interp(xrange,farray[0,:,i],farray[j+1,:,i])
maxint = np.amax(yarray[:,:])

print("Avg: {}, stdev: {}".format(np.mean(yarray[:,0]),np.std(yarray[:,0])))

# Rescaling
if args.br:
  ymin = args.br[0]
  ymax = args.br[1]
  print("Using fixed ymin: {}, ymax: {}".format(ymin,ymax))
else:
  ymin = np.amin(yarray[:,0])
  ymax = np.amax(yarray[:,0])
  print("Using automatic ymin: {}, ymax: {}".format(ymin,ymax))
yrng = ymax - ymin
nmaxint = max([max(yarray[i,:]-yarray[i,0]) for i in range(0,nspec*nfiles)])
for i in range(0,nspec*nfiles):
  yarray[i,:] = yarray[i,0] + (yarray[i,:] - yarray[i,0]) * yrng/(2*nmaxint)
maxint = np.amax(yarray[:,:])

# Doing block averaging
narray=np.zeros([n,npts])
jcnt = np.zeros([n])
for i in range(0,nspec*nfiles):
  for j in range(0,n):
    if truncate(yarray[i,0]) <= truncate(ymin + yrng*(j+1)/n):
      narray[j,:] = (narray[j,:]*jcnt[j] + yarray[i,:])/(jcnt[j]+1)
      jcnt[j]+=1
      break

# If the -o flag is on, we shift all the blocks so they start at zero
carray = np.copy(narray)
if args.o:
  for i in range(0,n):
    narray[i,:] -= narray[i,0]
  ymin = np.amin(narray[:,0])
  maxint = np.amax(narray[:,:])

#####
# Initializing plot

titlesize=15
axissize=35
ticksize=30
textsize=20

# General plot parameters
plt.figure(figsize=(21,12))
ax = plt.subplot()
plt.xlabel('Energy (eV)',size=axissize)
xmin = 403
xmax = 415
plt.xlim(xmin,xmax)
plt.xticks(size=ticksize)
ax.set_xticks(minor=True,ticks=np.arange(np.ceil(xmin),np.ceil(xmax),2))
ax.set_xticks(minor=True,ticks=np.arange(np.ceil(xmin),np.ceil(xmax),1))
plt.tick_params(axis='x',which='major',direction='in',length=8)
plt.tick_params(axis='x',which='minor',direction='in',length=5)
ylen = maxint - ymin
plt.ylim(ymin-ylen*0.05,maxint+ylen*0.05)

# Plotting main data and block averages without the -o flag:
if not args.o:
  for i in range(0,nfiles*nspec):
    ax.plot(xrange,yarray[i,:],linewidth=1)

  for i in range(0,n):
    ax.plot(xrange,narray[i,:],color='black',linewidth=8)
    ax.axhline(narray[i,0],color='black',linewidth=2,linestyle='--')
    if n != 1:
      ax.axhline(ymin+yrng*i/n,color='black')
      ax.axhline(ymin+yrng*(i+1)/n,color='black')

  tshift = (xmax-xmin)*0.025

else: # Otherwise plotting just the block averages
  for i in range(0,n):
    if jcnt[i] >= int(args.c[0]):
      ax.plot(xrange,narray[i,:],label='{:.2f} - {:.2f} Å ({:.0f} spectra)'.format(args.br[0]+i*yrng/n,args.br[0]+(i+1)*yrng/n,jcnt[i]),linewidth=8)
      plt.legend(prop={'size': 25})
  plt.yticks([])

plt.ylabel('Average N-C distance (Å)',size=axissize)
plt.yticks(size=ticksize)
plt.tick_params(axis='y',which='major',direction='inout',length=12)
cwds = os.getcwd().split('/')

ax2 = ax.twinx()
ax2.set_ylabel('Intensity (arbitrary units)',size=axissize)
ax2.set_yticks([])

# Printing the blocks
for i in range(0,n):
  with open('comb-n{}-block{}.blk'.format(n,i),'w') as fp:
    for j in range(0,npts):
      fp.write("{} {}\n".format(xrange[j],narray[i,j]))
  print("Block {} has {:.0f} spectra averaged.".format(i+1,jcnt[i]))
print("Block total: {:.0f} spectra.".format(sum(jcnt[:])))

if args.s:
  if int(args.c[0]) > 1:
    plt.savefig("{}/{}-{}-n{}-c{}.pdf".format(os.getcwd(),cwds[-2],cwds[-1],n,args.c[0]),bbox_inches='tight')
  else:
    plt.savefig("{}/{}-{}-n{}.pdf".format(os.getcwd(),cwds[-2],cwds[-1],n),bbox_inches='tight')
plt.show()
