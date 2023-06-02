#!/usr/bin/env python3

import sys
import numpy as np
from glob import glob
import argparse
from scipy.special import wofz

#####
# Explaining program

# Parsing arguments
p = argparse.ArgumentParser(description="Generates convoluted Gaussian curves for each *xas*.spectrum file in the directory.")
p.add_argument('range',nargs=2,default=[0,1000],type=int,help="Convolutes spectra for atoms in range I to J",metavar='I J')
p.add_argument('-s','--sigma',default=200,type=float,help="value of sigma (in meV) for Gaussian curve",metavar='')
p.add_argument('-l','--lines',default=1000,type=int,help="number of lines from .spectrum file to be read",metavar='')
p.add_argument('-n','--numpts',default=100,type=int,help="number of points per Gaussian curve",metavar='')
group = p.add_mutually_exclusive_group()
group.add_argument('-r',action='store_false',help="print only the individual range spectra, not the summed one")
group.add_argument('-t',action='store_false',help="print only the summed spectrum, not the individual ones")
p.add_argument('-c','--cutoff',nargs=1,help="The number of sigma for cutoff of the Gaussian curve for the convolution",metavar='N')
args = p.parse_args()

# Setting parameters
arange = args.range
print("Range of atoms: {} to {}".format(arange[0],arange[1]))
sigma = args.sigma
nlines = args.lines
npts = args.numpts
print("Running with values: sigma = {} meV, lines = {}, numpts = {}".format(sigma,nlines,npts))

#####
# Defining functions
def convolute(xvals,yvals,sigma,npts):
  sigma = sigma/1000
  lower = min(xvals)-3*sigma
  upper = max(xvals)+3*sigma
  dx = 6*sigma/npts
  tpts = int((upper-lower)/dx)

  garray = np.zeros([tpts,2])
  for i in range(0,tpts):
    x = lower + i*dx
    garray[i,0] = x
    for j in range(0,len(yvals)):
      if args.cutoff:
        if np.abs(x-xvals[j])/sigma <= float(args.cutoff[0]):
          garray[i,1] = garray[i,1] + (1/(sigma*np.sqrt(2*np.pi)))*yvals[j]*np.exp(-0.5*(((x-xvals[j])/sigma)**2))
      else:
        garray[i,1] = garray[i,1] + (1/(sigma*np.sqrt(2*np.pi)))*yvals[j]*np.exp(-0.5*(((x-xvals[j])/sigma)**2)) # Gaussian
  return garray

#####
# Running program

# Getting all the XAS spectra
flist = []
for i in range(arange[0],arange[1]+1):
  f = glob('*xas*at{}*.spectrum'.format(i))
  if len(f) > 1:
    print("Too many files found: {}".format(f))
    exit()
  elif len(f) == 1:
    flist.append(f[0])
nspec = len(flist)

iarray = np.zeros([nlines,2,nspec])
# Going through each spectrum, reading the spectral peak data into array
for i in range(0,nspec):
  f = flist[i]
  with open(f) as fp:
    linecnt = 0
    fp.readline()
    for line in fp:
      if linecnt < nlines:
        x = line.split()[1]
        y = line.split()[5]
        iarray[linecnt,0,i] = x
        iarray[linecnt,1,i] = y
        linecnt+=1

# Now with all the data, going through each spectrum and convoluting the Gaussian
if args.t:
  for i in range(0,nspec):
    conv = convolute(iarray[:,0,i],iarray[:,1,i],sigma,npts)
    print("Convoluting spectrum for {}".format(flist[i]))
    np.savetxt(flist[i].replace('.spectrum','-s{}-l{}-n{}.dat'.format(sigma,nlines,npts)),conv)

# Compiling all the data into a huge list for sum convolution
tarray = np.zeros([nlines*nspec,2])
for i in range(0,nspec):
  tarray[i*nlines:(i+1)*nlines,0] = iarray[:,0,i]
  tarray[i*nlines:(i+1)*nlines,1] = iarray[:,1,i]

if args.r:
  print("Now convoluting the summed spectrum for the range {} to {}".format(arange[0],arange[1]))
  conv = convolute(tarray[:,0],tarray[:,1],sigma,npts)
  conv[:,1] /= nspec # normalising the array to the number of spectra
  np.savetxt('sum-conv-range{}to{}-s{:g}-l{}-n{}.dat'.format(arange[0],arange[1],sigma,nlines,npts),conv)
