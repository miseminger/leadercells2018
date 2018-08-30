import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os
import sys
from optparse import OptionParser


## OPTIONS PARSING
usage = '%prog options\n  > one input CSV is required'
parser = OptionParser(usage=usage)
parser.add_option('-f', type='string', dest='datafile', help='path to clustering data CSV file')

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
request = "Your input CSV contains data for frames 0 to " + str(int(max_time)) + ". \nEnter the frame numbers you want to include in the box plot figure, separated by commas: "

while True:
    times = input(request)
    if all(i <= max_time for i in times):
        break
    else:
        print ("At least one of your frame numbers is out of range.  Please try again.")

df = pd.read_csv(options.datafile, usecols=["Metadata_FrameNumber", "Max_Dist / Min_Dist"])
df['Max_Dist / Min_Dist'] = df['Max_Dist / Min_Dist']**(-1)  #convert max/min to min/max

experiment_index = options.datafile.find('EXP') 
experiment = options.datafile[experiment_index:experiment_index+5] #experiment name to include in output figname

underscores = [pos for pos, char in enumerate(options.datafile) if char == '_']
s = underscores[1]
e = underscores[6]
detail = options.datafile[s+1:e] #max angle and radius to include in figname


#make box-and-whisker plot

#get boxplot data first
def timechunk(x):
	return df[df["Metadata_FrameNumber"]==int(x)]

boxplotdata = pd.concat([timechunk(x) for x in times])
boxplotdata.Metadata_FrameNumber = boxplotdata.Metadata_FrameNumber.astype(int)

fig = plt.figure()
boxplot = sns.boxplot(x='Metadata_FrameNumber', y='Max_Dist / Min_Dist', data=boxplotdata, color='#bc23d3')
plt.xlabel('frame', fontsize='x-large')
plt.ylabel('aspect ratio', fontsize='x-large')
ax = fig.add_subplot(1, 1, 1)
ax.tick_params(axis='both', labelsize='x-large')
figname = experiment + '_aspect_ratio_boxplot_' + detail + '_1.png'
count = 1
while os.path.isfile(figname): #if the csv name already exists, make new files with _1, _2, _3 at the end to differentiate
    count += 1
    figname = experiment + '_aspect_ratio_boxplot_'  + detail + '_' + str(count) + '.png'
fig = boxplot.get_figure()
fig.savefig(figname) 
print 'Boxplot produced: ' + figname

#create errorbar figure
frames = df['Metadata_FrameNumber'].unique()

means = df.groupby(['Metadata_FrameNumber'],as_index=False)[['Max_Dist / Min_Dist']].mean().values[:,1]
stds = df.groupby(['Metadata_FrameNumber'],as_index=False)[['Max_Dist / Min_Dist']].std().values[:,1]

fig = plt.figure()
plt.errorbar(frames, means, yerr=stds, fmt='none', alpha=0.5, ecolor='#5a1a75')
plt.scatter(frames, means, s=10, color='#5a1a75')
plt.xlabel('frame', fontsize='x-large')
plt.ylabel('mean aspect ratio', fontsize='x-large')
ax = fig.add_subplot(1, 1, 1)
ax.tick_params(axis='both', labelsize='x-large')
plt.xlim(-20,max_time+20)

figname2 = experiment + '_meanaspectratios_' + detail + '_1.png'  
count = 1
while os.path.isfile(figname): #if the csv name already exists, make new files with _1, _2, _3 at the end to differentiate
    count += 1
    figname = experiment + '_meanaspectratios_' + detail + '_' + str(count) + '.png'
fig.savefig(figname2, bbox_inches='tight') #bbox_inches prevents the legend from being cut off

print 'Errorbar plot produced: ' + figname2
