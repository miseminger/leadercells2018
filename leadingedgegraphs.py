import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import pandas as pd
from optparse import OptionParser
import seaborn as sns

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

columnnames = ('Metadata_FrameNumber','ObjectNumber','Homotypic Neighbour IDs','Heterotypic Neighbour IDs','Main_Location_Center_X (px)','Main_Location_Center_Y (px)','Homo_Location_Center_X (px)','Homo_Location_Center_Y (px)','Hetero_Location_Center_X (px)','Hetero_Location_Center_Y (px)')

try:
    greendf, reddf = pd.read_csv(options.greenfile, usecols=columnnames), pd.read_csv(options.redfile, usecols=columnnames)
except KeyError:
    sys.exit('Error: could not parse CSV, use -C if file is for CellProfiler (MATLAB files are default)\n'
             'Use -h for options usage help')

frames = greendf['Metadata_FrameNumber'].unique()
max_time = int(greendf['Metadata_FrameNumber'].max())

reddf = reddf.rename(index=str, columns={'Metadata_FrameNumber': "Frame"})
reddf = reddf.rename(index=str, columns={'ObjectNumber': "Leading Cell ObjectNumber"})

greendf = greendf.rename(index=str, columns={'Metadata_FrameNumber': "Frame"})
greendf = greendf.rename(index=str, columns={'ObjectNumber': "Leading Cell ObjectNumber"})

percent_red = float(options.percent_red)
percent_green = float(options.percent_green)

# Build directory to contain outputs
experiment_index = options.greenfile.find('EXP') 
experiment = options.greenfile[experiment_index:experiment_index+5] #experiment name to include in output filenames

underscores = [pos for pos, char in enumerate(options.greenfile) if char == '_']
s = underscores[2]
e = underscores[7]
detail = options.greenfile[s+1:e] #max angle and radius to include in figname

directory = experiment + '_leadingedge_graphs_' + detail
directory_i = 1
while os.path.exists(directory):
    directory_i += 1
    directory = experiment + '_leadingedge_graphs_' + detail + '_' + str(directory_i)
os.makedirs(directory)

print ' '
print 'Directory produced: ' + directory
print ' '


reddf[['Homotypic Neighbour IDs','Heterotypic Neighbour IDs']] = reddf[['Homotypic Neighbour IDs','Heterotypic Neighbour IDs']].astype(bool)
reddf[['Homotypic Neighbour IDs','Heterotypic Neighbour IDs']] = reddf[['Homotypic Neighbour IDs','Heterotypic Neighbour IDs']].astype(int)

rrdf = reddf.groupby(['Frame', 'Leading Cell ObjectNumber'],as_index=False)[['Homotypic Neighbour IDs']].sum()
rrdf = rrdf.drop(['Leading Cell ObjectNumber'], axis=1)
rrdf = rrdf.rename(index=str, columns={'Homotypic Neighbour IDs': "Neighbours"}) #rename 'Heterotypic Neighbour IDs' column as 'Neighbours'

rgdf = reddf.groupby(['Frame', 'Leading Cell ObjectNumber'],as_index=False)[['Heterotypic Neighbour IDs']].sum()
rgdf = rgdf.drop(['Leading Cell ObjectNumber'], axis=1)
rgdf = rgdf.rename(index=str, columns={'Heterotypic Neighbour IDs': "Neighbours"}) #rename 'Heterotypic Neighbour IDs' column as 'Neighbours'

greendf[['Homotypic Neighbour IDs','Heterotypic Neighbour IDs']] = greendf[['Homotypic Neighbour IDs','Heterotypic Neighbour IDs']].astype(bool)
greendf[['Homotypic Neighbour IDs','Heterotypic Neighbour IDs']] = greendf[['Homotypic Neighbour IDs','Heterotypic Neighbour IDs']].astype(int)

ggdf = greendf.groupby(['Frame', 'Leading Cell ObjectNumber'],as_index=False)[['Homotypic Neighbour IDs']].sum()
ggdf = ggdf.drop(['Leading Cell ObjectNumber'], axis=1)
ggdf = ggdf.rename(index=str, columns={'Homotypic Neighbour IDs': "Neighbours"}) #rename 'Heterotypic Neighbour IDs' column as 'Neighbours'

grdf = greendf.groupby(['Frame', 'Leading Cell ObjectNumber'],as_index=False)[['Heterotypic Neighbour IDs']].sum()
grdf = grdf.drop(['Leading Cell ObjectNumber'], axis=1) #finally, remove the 'Leading Cell ObjectNumber' column from each dataframe
grdf = grdf.rename(index=str, columns={'Heterotypic Neighbour IDs': "Neighbours"}) #rename 'Heterotypic Neighbour IDs' column as 'Neighbours'

ggdf2 = ggdf
rrdf2 = rrdf
rgdf2 = rgdf
grdf2 = grdf

#now normalize all the neighbour numbers
rrdf["Neighbours"] = rrdf["Neighbours"] * percent_green / percent_red
rgdf["Neighbours"] = rgdf["Neighbours"] * percent_green / percent_red
ggdf["Neighbours"] = ggdf["Neighbours"] * percent_red / percent_green
grdf["Neighbours"] = grdf["Neighbours"] * percent_red / percent_green



print 'Preparing errorbar plots...'
print ' '

#get means and stdevs per time
rrmeans = rrdf.groupby('Frame').mean().values #mean number of rrneighbours per cell per frame
rrstds = rrdf.groupby('Frame').std().values #stdev of number of rrneighbours per cell per frame

ggmeans = ggdf.groupby('Frame').mean().values #mean number of ggneighbours per cell per frame
ggstds = ggdf.groupby('Frame').std().values #stdev of number of ggneighbours per cell per frame

rgmeans = rgdf.groupby('Frame').mean().values #mean number of rgneighbours per cell per frame
rgstds = rgdf.groupby('Frame').std().values #stdev of number of rgneighbours per cell per frame

grmeans = grdf.groupby('Frame').mean().values #mean number of grneighbours per cell per frame
grstds = grdf.groupby('Frame').std().values #stdev of number of grneighbours per cell per frame

#collect data for boxplots (one for red cells, one for green, for 3 timepoints: start, middle, and end)
def getboxplotdata(dataframe, homo_or_hetero):
	first =  dataframe[dataframe['Frame']==0]
	middle = dataframe[dataframe['Frame']==int(np.median(frames))]
	last = dataframe[dataframe['Frame']==max_time] 
	sections = [first, middle, last]
	combinedsections = pd.concat(sections) #return a pandas dataframe that's the same as the original but with only the selected times in it
	final = combinedsections.assign(Type=homo_or_hetero) #add a column for the category ('homotypic' or heterotypic')
	return final 

ggbox = getboxplotdata(ggdf, 'homotypic')
grbox = getboxplotdata(grdf, 'heterotypic')
both = [ggbox, grbox]
greenboxplotdata = pd.concat(both)

rrbox = getboxplotdata(rrdf, 'homotypic')
rgbox = getboxplotdata(rgdf, 'heterotypic')
both = [rrbox, rgbox]
redboxplotdata = pd.concat(both)

#create two errorbar figures
#first for the red leading edge cells
fig = plt.figure()
plt.errorbar(frames, rrmeans, yerr=rrstds, fmt='none', alpha=0.5, ecolor='r')
plt.scatter(frames, rrmeans, s=10, label='homotypic', color='r')
plt.errorbar(frames, rgmeans, yerr=rgstds, fmt='none', alpha=0.3, ecolor='k')
plt.scatter(frames, rgmeans, s=10, label='heterotypic', color='k')
plt.legend(loc=2, fontsize='x-large', ncol=2)
plt.xlabel('frame', fontsize='x-large')
plt.ylabel('mean neighbours', fontsize='x-large')
ax = fig.add_subplot(1, 1, 1)
ax.tick_params(axis='both', labelsize='x-large')
plt.xlim(-20,max_time+20)
 
#save red figure
redfigname = directory  + '/leadingedge_neighbours_red.png'  
plt.savefig(redfigname, bbox_inches='tight') #bbox_inches prevents the legend from being cut off


#then for the green leading edge cells
fig = plt.figure()
plt.errorbar(frames, ggmeans, yerr=ggstds, fmt='none', alpha=0.5, ecolor='c')
plt.scatter(frames, ggmeans, s=10, label='homotypic', color='c')
plt.errorbar(frames, grmeans, yerr=grstds, fmt='none', alpha=0.3, ecolor='k')
plt.scatter(frames, grmeans, s=10, label='heterotypic', color='k')
ax = fig.add_subplot(1, 1, 1)
ax.tick_params(axis='both', labelsize='x-large')
plt.legend(loc=2, fontsize='x-large', ncol=2)
plt.xlabel('frame', fontsize='x-large')
plt.ylabel('mean neighbours', fontsize='x-large')
plt.xlim(-20,max_time+20)

#save green figure
greenfigname = directory  + '/leadingedge_neighbours_green.png'
plt.savefig(greenfigname, bbox_inches='tight') #bbox_inches prevents the legend from being cut off

print 'Scatterplots of mean neighbours per cell produced: ' 
print greenfigname 
print redfigname
print ' '
print 'Preparing box-and-whisker plots...'
print ' '

#now for the boxplots!

#green cells first
#make box-and-whisker plot
fig = plt.figure()
green_pal = {"homotypic": '#5EECFD', "heterotypic": "k"}
boxplot = sns.boxplot(x='Frame', y='Neighbours', hue='Type', data=greenboxplotdata, palette=green_pal)
plt.legend(loc=1, fontsize='x-large')
plt.xlabel('frame', fontsize='x-large')
plt.ylabel('neighbours', fontsize='x-large')
ax = fig.add_subplot(1, 1, 1)
ax.tick_params(axis='both', labelsize='x-large')
#save figure
greenfigname = directory  + '/leadingedge_neighbours_boxplot_green.png'
boxplot.get_figure()
fig.savefig(greenfigname) 

#red cells next
fig = plt.figure()
red_pal = {"homotypic": '#FF9090', "heterotypic": "k"}
boxplot = sns.boxplot(x='Frame', y='Neighbours', hue='Type', data=redboxplotdata, palette=red_pal)
plt.legend(loc=1, fontsize='x-large')
plt.xlabel('frame', fontsize='x-large')
plt.ylabel('neighbours', fontsize='x-large')
ax = fig.add_subplot(1, 1, 1)
ax.tick_params(axis='both', labelsize='x-large')
#save figure
redfigname = directory  + '/leadingedge_neighbours_boxplot_red.png'
fig = boxplot.get_figure()
fig.savefig(redfigname) 

print 'Boxplots produced for 3 timepoints: ' 
print greenfigname 
print redfigname
print ' '
print 'Preparing 2d histograms...'
print ' '

print "number of NaNs:"

print 'ggdf', ggdf.isnull().sum().sum()
print 'grdf', ggdf.isnull().sum().sum()
print 'rrdf', ggdf.isnull().sum().sum()
print 'rgdf', ggdf.isnull().sum().sum()

ggdf['Neighbours'] =  ggdf['Neighbours'] * 1
grdf['Neighbours'] =  grdf['Neighbours'] * 1
rrdf['Neighbours'] =  rrdf['Neighbours'] * 1
rgdf['Neighbours'] =  rgdf['Neighbours'] * 1

print "datatypes"
print "ggdf", ggdf['Neighbours'].dtype
print "grdf", grdf['Neighbours'].dtype
print "rrdf", rrdf['Neighbours'].dtype
print "rgdf", rgdf['Neighbours'].dtype

fig = plt.figure()
gg = sns.jointplot('Frame', 'Neighbours', data=ggdf, kind="hex", space=0, color='#5EECFD')
gg.set_axis_labels('frame', 'neighbours', fontsize='x-large')
plt.tight_layout()
gg.savefig(directory  + '/gg_2dhistkde.png')
print directory  + '/gg_2dhistkde.png'

fig = plt.figure()
gr = sns.jointplot('Frame', 'Neighbours', data=grdf, kind="hex", space=0)
gr.set_axis_labels('frame', 'neighbours', fontsize='x-large')
plt.tight_layout()
gr.savefig(directory  + '/gr_2dhistkde.png')
print directory  + '/gr_2dhistkde.png'

fig = plt.figure()
rr = sns.jointplot('Frame', 'Neighbours', data=rrdf, kind="hex", space=0, color='#FF9090')
rr.set_axis_labels('frame', 'neighbours', fontsize='x-large')
plt.tight_layout()
rr.savefig(directory  + '/rr_2dhistkde.png')
print directory  + '/rr_2dhistkde.png'

fig = plt.figure()
rg = sns.jointplot('Frame', 'Neighbours', data=rgdf, kind="hex", space=0)
rg.set_axis_labels('frame', 'neighbours', fontsize='x-large')
plt.tight_layout()
rg.savefig(directory  + '/rg_2dhistkde.png')
print directory  + '/rg_2dhistkde.png'






