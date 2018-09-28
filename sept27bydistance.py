import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import pandas as pd
from optparse import OptionParser
from matplotlib.ticker import FormatStrFormatter

## OPTIONS PARSING
usage = '%prog options\n  > both files (-g, -r) are required\n'
parser = OptionParser(usage=usage)
parser.add_option('-g', type='string', dest='greenfile', help='path to green CSV file')
parser.add_option('-r', type='string', dest='redfile', help='path to red CSV file')
parser.add_option('-R', type='string', dest='percent_red', help='percentage of red cells in experiment')
parser.add_option('-G', type='string', dest='percent_green', help='percentage of green cells in experiment')
parser.add_option('-N', type='string', dest='n', help='number of bins')
parser.add_option('-t', type='string', dest='timeint', help='time interval between frames, in minutes')

(options, args) = parser.parse_args()
error = 0
if not(options.greenfile or options.redfile):
    error = 1
    print 'Error: no file input, use -g and/or -r'

def np_read_csv(file):
    da = np.genfromtxt(file, delimiter=',', names=True)  #make an array out of the .csv, called da. max_rows=5 is good for testing!
    return da.view((float, len(da.dtype.names))) #transform everything in the array to a float

try:
    greenfile, redfile = np_read_csv(options.greenfile), np_read_csv(options.redfile)
except KeyError:
    sys.exit('Error: could not parse CSV, use -C if file is for CellProfiler (MATLAB files are default)\n'
             'Use -h for options usage help')

n = int(options.n) #number of bins
percent_red = float(options.percent_red)
percent_green = float(options.percent_green)
timeint = float(options.timeint)

# Build directory to contain outputs
experiment_index = options.greenfile.find('EXP') 
experiment = options.greenfile[experiment_index:experiment_index+5] #experiment name to include in output figname

underscores = [pos for pos, char in enumerate(options.greenfile) if char == '_']
s = underscores[2]
e = underscores[7]
detail = options.greenfile[s+1:e] #max angle and radius to include in figname


dist = 1600/n  #number of horizontal pixels in each image segment
binedges = np.arange(0,1601,dist) #this will give the left edges of the bins to index

redtimeslist = redfile[:,0].tolist()  #list format of the red frame numbers (as many of each frame # as there are detected cells in it)
greentimeslist = greenfile[:,0].tolist()  #list format of the green frame numbers (as many of each frame # as there are detected cells in it)
max_time = int(greenfile[-1,0])


def timecollector(time):
	
	#sort red cells into bins based on horizontal position, being sure to count each cell only once
	redfirsttime = redtimeslist.index(time)  #index for first instance of frame number
	redlasttime = len(redtimeslist) - redtimeslist[::-1].index(time) - 1 #index for last instance of frame number
	uniqueindex_red = np.unique(redfile[redfirsttime:(redlasttime+1),1], return_index=True)[1] #this is an array of the unique red indices
	binlabels_red = np.digitize(redfile[redfirsttime:(redlasttime+1),4][uniqueindex_red],binedges) #put the cells into their bins

	#do the same thing for the green cells
	greenfirsttime = greentimeslist.index(time)  #index for first instance of frame number
	greenlasttime = len(greentimeslist) - greentimeslist[::-1].index(time) - 1 #index for last instance of frame number
	uniqueindex_green = np.unique(greenfile[greenfirsttime:(greenlasttime+1),1], return_index=True)[1] #this is an array of the unique green indices
	binlabels_green = np.digitize(greenfile[greenfirsttime:(greenlasttime+1),4][uniqueindex_green],binedges)  #outputs an array of bin indices in order of green_x; the first bin is index 1

	#initialize lists to keep the neighbours per cell within this time in
	rrneighbours = []
	rgneighbours = []
	ggneighbours = []
	grneighbours = []


	#now go through the red cells and add their neighbours to the initialized lists
	for m in (np.array(range(int(redfile[(redlasttime),1]))) + 1): #iterate through red cell ids in timeframe, beginning with 1
	        cell_indices = np.where(redfile[redfirsttime:(redlasttime+1),1] == m)[0] #find all the indices of that cell id
		rrneighbours.append(np.sum((redfile[redfirsttime:(redlasttime+1),2][cell_indices] > 0))) #sum up the red-red neighbours for cell m and append this to rrneighbours
		rgneighbours.append(np.sum((redfile[redfirsttime:(redlasttime+1),3][cell_indices] > 0))) #sum up the red-green neighbours for cell m and append this to rgneighbours

	#now do the same for the green cells
	for m in (np.array(range(int(greenfile[(greenlasttime),1]))) + 1): #iterate through green cell ids in timeframe, beginning with 1
	        cell_indices = np.where(greenfile[greenfirsttime:(greenlasttime+1),1] == m)[0] #find all the indices of that cell id
		ggneighbours.append(np.sum((greenfile[greenfirsttime:(greenlasttime+1),2][cell_indices] > 0))) #sum up the red-red neighbours for cell m and append this to rrneighbours
		grneighbours.append(np.sum((greenfile[greenfirsttime:(greenlasttime+1),3][cell_indices] > 0))) #sum up the red-green neighbours for cell m and append this to grneighbours

	rrneighbours = np.array(rrneighbours)
	ggneighbours = np.array(ggneighbours)
	rgneighbours = np.array(rgneighbours)
	grneighbours = np.array(grneighbours)

	#now sort the neighbour sums into the correct bins for each cell
	#in the end, for this timeblock, you will get four arrays, three columns each: from left to right: frame, bin, and number of neighbours per cell

	def bincompiler(neighboursarr, maincolor, b):

		if maincolor == 'green':
			ratio = percent_red / percent_green
			binlabels = binlabels_green
		else:
			ratio = percent_green / percent_red
			binlabels = binlabels_red

		if sum(neighboursarr[binlabels == b]) != 0:
			arrayb = np.zeros((neighboursarr[binlabels == int(b)].shape[0],3))
			arrayb[:,0] = time
			arrayb[:,1] = b	
			arrayb[:,2] = neighboursarr[binlabels == int(b)] * ratio
		else:
			arrayb = np.ones((1,3)) * np.NaN
		
		return arrayb
	
	#collect blocks for all bins at current time	
	gg = np.concatenate([bincompiler(ggneighbours, 'green', b) for b in (np.array(range(n)) + 1)], axis=0) 
	gr = np.concatenate([bincompiler(grneighbours, 'green', b) for b in (np.array(range(n)) + 1)], axis=0) 
	rr = np.concatenate([bincompiler(rrneighbours, 'red', b) for b in (np.array(range(n)) + 1)], axis=0) 
	rg = np.concatenate([bincompiler(rgneighbours, 'red', b) for b in (np.array(range(n)) + 1)], axis=0) 
	return gg, gr, rr, rg


frames = np.array(range(max_time + 1)) #array of frame numbers for x axis

alltimes_gg = np.concatenate([timecollector(x)[0] for x in range(max_time + 1)], axis=0) #get the array 'collector' for all times and put these together
alltimes_gr = np.concatenate([timecollector(x)[1] for x in range(max_time + 1)], axis=0) #get the array 'collector' for all times and put these together
alltimes_rr = np.concatenate([timecollector(x)[2] for x in range(max_time + 1)], axis=0) #get the array 'collector' for all times and put these together
alltimes_rg = np.concatenate([timecollector(x)[3] for x in range(max_time + 1)], axis=0) #get the array 'collector' for all times and put these together

cols = ('Frame','Bin','Neighbours')
color_idx = np.linspace(0, 1, n)

def make_figure(arr, kind, ax):

	df = pd.DataFrame(data=arr, index=arr[:,0], columns=cols)
	df = df.dropna(how='all') 
	meanarr = df.groupby(['Frame', 'Bin'],as_index=False)[['Neighbours']].mean().values
	stdarr = df.groupby(['Frame', 'Bin'],as_index=False)
	stds = stdarr['Neighbours'].apply(lambda x: x.sem()).values #get standard error (not actually stdev)
	df = pd.DataFrame(data=meanarr, index=meanarr[:,0], columns=('Frame','Bin','Mean'))
	df['StDev'] = stds

	for b, i in zip((np.array(range(n)) + 1), color_idx):
		means = df[df['Bin'] == float(b)]['Mean'].values
		stds = df[df['Bin'] == float(b)]['StDev'].values
		frames = df[df['Bin'] == float(b)]['Frame'].values * timeint / 60.0
		ax.errorbar(frames, means, yerr=stds, fmt='none', alpha=1.0, ecolor=plt.cm.gist_rainbow(i))

	for b, i in zip((np.array(range(n)) + 1), color_idx):
		means = df[df['Bin'] == float(b)]['Mean'].values
		frames = df[df['Bin'] == float(b)]['Frame'].values * timeint / 60.0
		ax.scatter(frames, means, s=10, label=('bin ' + str(int(b))), color=plt.cm.gist_rainbow(i))

	ax.set_xlabel('frame', fontsize='x-large')
	ax.set_ylabel('mean neighbours per cell', fontsize='x-large')
	ax.tick_params(axis='both', labelsize='x-large')
	ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
	ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
	ax.set_xlim(-1,(max_time * timeint / 60.0 + 1))
	if ax == ax2:
		ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize='x-large', ncol=1)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=False, sharey=True, figsize=(16,12))
make_figure(alltimes_gg, 'gg', ax1)
make_figure(alltimes_gr, 'gr', ax2)
make_figure(alltimes_rr, 'rr', ax3)
make_figure(alltimes_rg, 'rg', ax4)
	
figname = experiment + '_mean_neighbours_per_cell_by_distance_' + detail + '_' + str(n) + '_bins_1.png'
count = 1
while os.path.isfile(figname): #if the csv name already exists, make new files with _1, _2, _3 at the end to differentiate
   	count += 1
  	figname = experiment + '_mean_neighbours_per_cell_by_distance_' + detail + '_' + str(n) + '_bins_' + str(count) + '.png'
fig.savefig(figname, bbox_inches='tight') #bbox_inches prevents the legend from being cut off

print "Figure produced: " + figname



