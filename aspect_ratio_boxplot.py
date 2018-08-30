#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 18:15:25 2018

@author: madeline
"""
import matplotlib
matplotlib.use('Agg')

import numpy as np
import seaborn as sns
import pandas as pd
import os
import sys
from optparse import OptionParser


## OPTIONS PARSING
usage = '%prog options\n  > one input CSV is required'
parser = OptionParser(usage=usage)
parser.add_option('-d', type='string', dest='datafile', help='path to clustering data CSV file')

(options, args) = parser.parse_args()
error = 0

if not(options.datafile):
    error = 1
    print 'Error: input data file must be provided using -d; use a file of the type produced by aspect_ratio.py (eg. \'clusteringdata_40um_radius_1.csv\')'

def np_read_csv(file):
    da = np.genfromtxt(file, delimiter=',', names=True)  #make an array out of the .csv, called da. max_rows=5 is good for testing!
    return da.view((float, len(da.dtype.names))) #transform everything in the array to a float

try:
    data = np_read_csv(options.datafile)
except KeyError:
    sys.exit('Error: could not parse CSV; please check your file and try again.')
    
max_time = data[-1,0]
data = np.delete(data, np.where((data[:,12] + data[:,13]) <= 2), axis=0) #get rid of all two-cell clusters in dataset, as these skew the aspect ratios
request = "Your input CSV contains data for frames 0 to " + str(int(max_time)) + ". \nEnter the frame numbers you want to include in the box plot figure, separated by commas: "

while True:
    times = input(request)
    if all(i <= max_time for i in times):
        break
    else:
        print ("At least one of your frame numbers is out of range.  Please try again.")
    
times = np.array(times)

aspectratios = np.zeros((times.shape[0],))

def timecollector(time):
    
    timeslist = data[:,0].tolist()  #list format of the frame numbers (as many of each frame # as there are cells in it)
    firsttime = timeslist.index(time)  #index for first instance of frame number
    lasttime = len(timeslist) - timeslist[::-1].index(time) - 1 #index for last instance of frame number
    
    blocklength = lasttime - firsttime + 1
    aspectratio = np.zeros((blocklength,2))
    aspectratio[:,0] = int(time)
    aspectratio[:,1] = data[firsttime:(lasttime + 1),10] ** -1 #plot min_dist / max_dist, not the other way around! (0 to 1 easier to read)
    return aspectratio

violinplotdata = np.concatenate([timecollector(x) for x in times])

#makeviolinplotdata into a pandas dataframe so seaborn can read it
columnlabels = ['Frame', 'Min_Dist / Max_Dist']
df = pd.DataFrame(violinplotdata, columns=columnlabels)

#make violin-and-whisker plot
violinplot = sns.boxplot(x='Frame', y='Min_Dist / Max_Dist', data=df)

#save figure
figname = 'aspect_ratio_violinplot_1.png'

count = 1
while os.path.isfile(figname): #if the csv name already exists, make new files with _1, _2, _3 at the end to differentiate
    count += 1
    figname = 'aspect_ratio_violinplot_' + str(count) + '.png'
    
fig = violinplot.get_figure()
fig.savefig(figname) 

print 'Figure produced: ' + figname


