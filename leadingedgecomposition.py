import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from optparse import OptionParser

## OPTIONS PARSING
usage = '%prog options\n  > both files (-g, -r) are required\n'
parser = OptionParser(usage=usage)
parser.add_option('-g', type='string', dest='greenfile', help='path to green CSV file')
parser.add_option('-r', type='string', dest='redfile', help='path to red CSV file')
parser.add_option('-R', type='string', dest='percent_red', help='percentage of red cells in experiment')
parser.add_option('-G', type='string', dest='percent_green', help='percentage of green cells in experiment')

(options, args) = parser.parse_args()
error = 0
if not(options.greenfile or options.redfile):
    error = 1
    print 'Error: no file input, use -g and/or -r'

colnames = ['Metadata_FrameNumber','ObjectNumber']

try:
    greenfile, redfile = pd.read_csv(options.greenfile, usecols=colnames), pd.read_csv(options.redfile, usecols=colnames)
except KeyError:
    sys.exit('Error: could not parse CSV, use -C if file is for CellProfiler (MATLAB files are default)\n'
             'Use -h for options usage help')

percent_red = float(options.percent_red)
percent_green = float(options.percent_green)

# Build directory to contain outputs
experiment_index = options.greenfile.find('EXP') 
experiment = options.greenfile[experiment_index:experiment_index+5] #experiment name to include in output filenames

underscores = [pos for pos, char in enumerate(options.greenfile) if char == '_']
s = underscores[2]
e = underscores[7]
detail = options.greenfile[s+1:e] #max angle and radius to include in figname

uniquereds = redfile.groupby(['Metadata_FrameNumber'],as_index=False)['ObjectNumber'].nunique().values #this is a vector containing the number of leading cells per frame
uniquegreens = greenfile.groupby(['Metadata_FrameNumber'],as_index=False)['ObjectNumber'].nunique().values #this is a vector containing the number of leading cells per frame

uniquereds = uniquereds.astype(float) * percent_green / percent_red
uniquegreens = uniquegreens.astype(float) * percent_red / percent_green

redproportion = uniquereds / (uniquereds + uniquegreens) 
greenproportion = uniquegreens / (uniquereds + uniquegreens) 
frames = greenfile['Metadata_FrameNumber'].unique()
max_time = int(greenfile['Metadata_FrameNumber'].max())

fig = plt.figure()
plt.scatter(frames, greenproportion, s=10, label='GFP', color='c')
plt.scatter(frames, redproportion, s=10, label='mCherry', color='r')
plt.legend(loc=0, fontsize='x-large', ncol=2)
plt.xlabel('frame', fontsize='x-large')
plt.ylabel('proportion of leading edge cells', fontsize='x-large')
ax = fig.add_subplot(1, 1, 1)
ax.tick_params(axis='both', labelsize='x-large')
plt.xlim(-20,max_time+20)
plt.ylim(0,1.2)

figname = experiment + '_leading_edge_composition_' + detail + '_1.png'
count = 1
while os.path.isfile(figname): #if the csv name already exists, make new files with _1, _2, _3 at the end to differentiate
    count += 1
    figname = experiment + '_leading_edge_composition_'  + detail + '_' + str(count) + '.png'
fig.savefig(figname) 
print 'Figure produced: ' + figname



