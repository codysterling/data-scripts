#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
import os

#####
# Initializing program

p = argparse.ArgumentParser(description="For nicely plotting data, this assumes all files are same structure (# rows and columns)")
p.add_argument('-f','--files',nargs='+',required=True,help="List of files to plot")
p.add_argument('-l','--labels',nargs='+',required=True,help="List of labels for each file in order")
p.add_argument('-x','--xlabel',nargs=1,required=True,help="X-axis label")
p.add_argument('-y','--ylabel',nargs=1,required=True,help="Y-axis label")
p.add_argument('-n','--name',nargs=1,help="Filename to save file")
args = p.parse_args()

files = args.files
labels = args.labels

if len(files) != len(labels):
  print("Number of files ({}) != number of labels ({})!  Exiting...".format(len(files),len(labels)))
  exit()

#####
# Defining functions

def nlines(file):
  count = 0
  with open(file) as fp:
    for line in fp:
      count+=1
  return count

#####
# Running program

# Get number of lines in the file

farray = np.empty([101,3,len(files)])

# Reading the files and getting the info
for index,file in enumerate(files):
  lcnt = 0
  with open(file) as fp:
    for line in fp:
      words = line.split()
      if len(words) == 2:
        words = np.append(words,0)
      farray[lcnt,:,index] = words
      lcnt += 1

titlesize=15
axissize=35
ticksize=30
textsize=20

plt.figure(figsize=(21,12))
ax = plt.subplot()
plt.xlabel("{}".format(args.xlabel[0]),size=axissize)

xmin = 0
xmax = 10
plt.xlim(xmin,xmax)
plt.xticks(size=ticksize)
ax.set_xticks(minor=True,ticks=np.arange(np.ceil(xmin),np.ceil(xmax),2))
ax.set_xticks(minor=True,ticks=np.arange(np.ceil(xmin),np.ceil(xmax),1))
plt.tick_params(axis='x',which='major',direction='in',length=8)
plt.tick_params(axis='x',which='minor',direction='in',length=5)

plt.ylabel("{}".format(args.ylabel[0]),size=axissize)
plt.yticks(size=ticksize)
plt.ylim(-0.025,1.8)

for i in range(farray.shape[2]):
  x = farray[:,0,i]
  y = farray[:,1,i]
  plt.plot(farray[:,0,i],farray[:,1,i],linewidth=8,label="{}".format(labels[i]))

plt.legend(prop={'size': axissize})

# MAPBr lines
for v in [2.8766003514497287, 3.8705422849559663, 5.030390809965464, 5.818883082468649]:
  plt.axvline(v,c='C0',lw=5,ls='--')

# MAPI lines
for v in [3.077000905114847, 4.114970261984184, 5.255137612921149, 6.204658360167027]:
  plt.axvline(v,c='C1',lw=5,ls='--')

if args.name:
  cwds = os.getcwd()
  plt.savefig("{}/{}".format(cwds,args.name[0]),bbox_inches='tight')
plt.show()
