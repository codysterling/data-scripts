#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import argparse
import re
import os

#####
# Initializing program

p = argparse.ArgumentParser(description="Plots X/Y data from files, where later files plot over earlier files.  The format of the files is two columns in X Y order")
p.add_argument('files',nargs='+',help="The list of data files to be plotted.")
p.add_argument('-f',action='store_true',help="Enables lines of best fit for the data sets")
args = p.parse_args()

#####
# Defining functions

def nlines(file):
  count = 0
  with open(file) as fp:
    for line in fp:
      count+=1
  return count

def nround(x, base=.2):
  return round(base * round(float(x)/base),3)

#####
# Running program

lowx = lowy = 100
hix = hiy = 0
plt.figure(figsize=(21,12))
ax = plt.subplot()

# Plotting the data from each file
for index,file in enumerate(args.files):
  xl = []
  yl = []
  print("file is {}, nlines: {}".format(file,nlines(file)))
  with open(file) as fp:
    for line in fp:
      words = line.split()
      if len(words) != 2 and len(words) != 0:
        print("Line doesn't have two columns: {}".format(line))
        exit()
      x = float(words[0])
      y = float(words[1])
      xl.append(x)
      yl.append(y)
  if len(xl) != len(yl):
    print("xl ({}): {} and yl ({}): {} aren't the same!".format(len(xl),xl,len(yl),yl))
  if min(xl) < lowx: lowx = min(xl)
  if max(xl) > hix: hix = max(xl)
  if min(yl) < lowy: lowy = min(yl)
  if max(yl) > hiy: hiy = max(yl)
  plt.scatter(xl,yl,s=100,alpha=1-index*0.5)

# Line of best fit
import scipy.stats
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(xl, yl)
print(slope, intercept, r_value, p_value, std_err)
plt.plot(np.unique(xl),slope*np.unique(xl)+intercept,color='black',lw=8)

# Setting plot parameters
lowx -= 0.2
hix += 0.2
lowy -= 0.05
hiy += 0.05
rlowx = nround(lowx,0.2)
rhix = nround(hix,0.2)
rlowy = nround(lowy,0.1)
rhiy = nround(hiy,0.1)
rlowx, rhix = [1.4,1.65]
rlowy, rhiy = [1,1.12]
lowx, hix = [1.4,1.65]
lowy, hiy = [1,1.12]

titlesize=15
axissize=35
ticksize=30
textsize=20

# best fit text
plt.text(1.5825,1.038,"y = {:4f}x + {:4f}".format(slope,intercept),size=textsize)
plt.text(1.61,1.033,"r² = {:4f}".format(r_value**2),size=textsize)

plt.xlabel("N-C distance (Å)",size=axissize)
plt.xlim(lowx,hix)
plt.xticks(np.append(np.arange(rlowx,rhix,0.05),rhix)[1:],size=ticksize)
ax.set_xticks(minor=True,ticks=np.append(np.arange(rlowx,rhix,0.1),rhix))
plt.tick_params(axis='x',which='major',direction='in',length=8)
plt.tick_params(axis='x',which='minor',direction='in',length=5)

plt.ylabel("N-H distance (Å)",size=axissize)
plt.ylim(lowy,hiy)
plt.yticks(np.append(np.arange(rlowy,rhiy,0.02),rhiy),size=ticksize)
plt.tick_params(axis='y',which='major',direction='in',length=8)
plt.tick_params(axis='y',which='minor',direction='in',length=5)

cwds = os.getcwd()
name = "cnnh.pdf"
plt.savefig("{}/{}".format(cwds,name),bbox_inches='tight')
plt.show()
