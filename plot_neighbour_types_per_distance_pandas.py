import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import pandas as pd
from optparse import OptionParser

## OPTIONS PARSING
usage = '%prog options\n  > both files (-g, -r) are required\n'
parser = OptionParser(usage=usage)
parser.add_option('-g', type='string', dest='greenfile', help='path to green CSV file')
parser.add_option('-r', type='string', dest='redfile', help='path to red CSV file')
#parser.add_option('-R', type='string', dest='percent_red', help='percentage of red cells in experiment')
#parser.add_option('-G', type='string', dest='percent_green', help='percentage of green cells in experiment')
parser.add_option('-N', type='string', dest='n', help='number of bins')

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
#percent_red = float(options.percent_red)
#percent_green = float(options.percent_green)

# Build directory to contain outputs
experiment_index = options.greenfile.find('EXP') 
experiment = options.greenfile[experiment_index:experiment_index+5] #experiment name to include in output filenames

# Build directory to contain outputs
directory = experiment + '_neighbours_' + str(n) + '_bins'
directory_i = 1
while os.path.exists(directory):
    directory_i += 1
    directory = experiment + '_neighbours_' + str(n) + '_bins' + '_' + str(directory_i)
os.makedirs(directory)

dist = 1600/n  #number of horizontal pixels in each image segment
binedges = np.arange(0,1601,dist) #this will give the left edges of the bins to index

redtimeslist = redfile[:,0].tolist()  #list format of the red frame numbers (as many of each frame # as there are detected cells in it)
greentimeslist = greenfile[:,0].tolist()  #list format of the green frame numbers (as many of each frame # as there are detected cells in it)
max_time = int(greenfile[-1,0])


cols = ('Metadata_FrameNumber','ObjectNumber','Homotypic Neighbour IDs','Heterotypic Neighbour IDs','Main_Location_Center_X (px)','Main_Location_Center_Y (px)','Homo_Location_Center_X (px)','Homo_Location_Center_Y (px)','Hetero_Location_Center_X (px)','Hetero_Location_Center_Y (px)')

reddf = pd.DataFrame(data=redfile[:,:], index=redfile[:,0], columns=cols)
reddf[['Homotypic Neighbour IDs','Heterotypic Neighbour IDs']] = reddf[['Homotypic Neighbour IDs','Heterotypic Neighbour IDs']].astype(bool)
reddf[['Homotypic Neighbour IDs','Heterotypic Neighbour IDs']] = reddf[['Homotypic Neighbour IDs','Heterotypic Neighbour IDs']].astype(int)

rrdf = reddf.groupby(['Frame', 'Leading Cell ObjectNumber'],as_index=False)[['Homotypic Neighbour IDs']].sum()
rrdf = rrdf.drop(['Leading Cell ObjectNumber'], axis=1)
rrdf = rrdf.rename(index=str, columns={'Homotypic Neighbour IDs': "Neighbours"}) #rename 'Heterotypic Neighbour IDs' column as 'Neighbours'

print rrdf 

rgdf = reddf.groupby(['Frame', 'Leading Cell ObjectNumber'],as_index=False)[['Heterotypic Neighbour IDs']].sum()
rgdf = rgdf.drop(['Leading Cell ObjectNumber'], axis=1)
rgdf = rgdf.rename(index=str, columns={'Heterotypic Neighbour IDs': "Neighbours"}) #rename 'Heterotypic Neighbour IDs' column as 'Neighbours'

greendf = pd.DataFrame(data=greenfile[:,:], index=greenfile[:,0], columns=cols)
greendf[['Homotypic Neighbour IDs','Heterotypic Neighbour IDs']] = greendf[['Homotypic Neighbour IDs','Heterotypic Neighbour IDs']].astype(bool)
greendf[['Homotypic Neighbour IDs','Heterotypic Neighbour IDs']] = greendf[['Homotypic Neighbour IDs','Heterotypic Neighbour IDs']].astype(int)

ggdf = greendf.groupby(['Frame', 'Leading Cell ObjectNumber'],as_index=False)[['Homotypic Neighbour IDs']].sum()
ggdf = ggdf.drop(['Leading Cell ObjectNumber'], axis=1)
ggdf = ggdf.rename(index=str, columns={'Homotypic Neighbour IDs': "Neighbours"}) #rename 'Heterotypic Neighbour IDs' column as 'Neighbours'

grdf = greendf.groupby(['Frame', 'Leading Cell ObjectNumber'],as_index=False)[['Heterotypic Neighbour IDs']].sum()
grdf = grdf.drop(['Leading Cell ObjectNumber'], axis=1) #finally, remove the 'Leading Cell ObjectNumber' column from each dataframe
grdf = grdf.rename(index=str, columns={'Heterotypic Neighbour IDs': "Neighbours"}) #rename 'Heterotypic Neighbour IDs' column as 'Neighbours'

rrmeans = rrdf.groupby('Frame').sum().values #mean number of rrneighbours per cell per frame

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
		print m, "red-red"
		print (redfile[redfirsttime:(redlasttime+1),2][cell_indices] > 0)
		print np.sum((redfile[redfirsttime:(redlasttime+1),2][cell_indices] > 0))

	#now do the same for the green cells
	for m in (np.array(range(int(greenfile[(greenlasttime),1]))) + 1): #iterate through green cell ids in timeframe, beginning with 1
	        cell_indices = np.where(greenfile[greenfirsttime:(greenlasttime+1),1] == m)[0] #find all the indices of that cell id
		ggneighbours.append(np.sum((greenfile[greenfirsttime:(greenlasttime+1),2][cell_indices] > 0))) #sum up the red-red neighbours for cell m and append this to rrneighbours
		grneighbours.append(np.sum((greenfile[greenfirsttime:(greenlasttime+1),3][cell_indices] > 0))) #sum up the red-green neighbours for cell m and append this to grneighbours
		print m, "green-green"
		print (greenfile[greenfirsttime:(greenlasttime+1),2][cell_indices] > 0)
		print np.sum(greenfile[greenfirsttime:(greenlasttime+1),2][cell_indices] > 0)

	rrneighbours = np.array(rrneighbours)
	ggneighbours = np.array(ggneighbours)
	rgneighbours = np.array(rgneighbours)
	grneighbours = np.array(grneighbours)

	print "rrneighbours"
	print rrneighbours
	print "ggneighbours"
	print ggneighbours

	#now sort the neighbour sums into the correct bins for each cell
	#in the end, for this timeblock, you will get four arrays, three columns each: from left to right: frame, bin, and number of neighbours per cell

	def bincompiler(neighboursarr, maincolor, b):

		if maincolor == 'green':
			#ratio = percent_red / percent_green
			binlabels = binlabels_green
		else:
			#ratio = percent_green / percent_red
			binlabels = binlabels_red

		if sum(neighboursarr[binlabels == b]) != 0:
			arrayb = np.zeros((neighboursarr[binlabels == int(b)].shape[0],3))
			arrayb[:,0] = time
			arrayb[:,1] = b	
			arrayb[:,2] = neighboursarr[binlabels == int(b)] #* ratio
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

alltimes_gg = np.concatenate([timecollector(x)[0] for x in range(1)], axis=0) #get the array 'collector' for all times and put these together
alltimes_gr = np.concatenate([timecollector(x)[1] for x in range(1)], axis=0) #get the array 'collector' for all times and put these together
alltimes_rr = np.concatenate([timecollector(x)[2] for x in range(1)], axis=0) #get the array 'collector' for all times and put these together
alltimes_rg = np.concatenate([timecollector(x)[3] for x in range(1)], axis=0) #get the array 'collector' for all times and put these together

cols = ('Frame','Bin','Neighbours')
color_idx = np.linspace(0, 1, n)

def make_figure(arr, kind):

	df = pd.DataFrame(data=arr, index=arr[:,0], columns=cols)
	df = df.dropna(how='all') 
	meanarr = df.groupby(['Frame', 'Bin'],as_index=False)[['Neighbours']].mean().values
	stdarr = df.groupby(['Frame', 'Bin'],as_index=False)
	stds = stdarr['Neighbours'].apply(lambda x: x.sem()).values #get standard error (not actually stdev)
	df = pd.DataFrame(data=meanarr, index=meanarr[:,0], columns=('Frame','Bin','Mean'))
	df['StDev'] = stds

	fig = plt.figure()

	for b, i in zip((np.array(range(n)) + 1), color_idx):
		means = df[df['Bin'] == float(b)]['Mean'].values
		stds = df[df['Bin'] == float(b)]['StDev'].values
		frames = df[df['Bin'] == float(b)]['Frame'].values
		plt.errorbar(frames, means, yerr=stds, fmt='none', alpha=1.0, ecolor=plt.cm.gist_rainbow(i))

	for b, i in zip((np.array(range(n)) + 1), color_idx):
		means = df[df['Bin'] == float(b)]['Mean'].values
		frames = df[df['Bin'] == float(b)]['Frame'].values
		plt.scatter(frames, means, s=10, label=('bin ' + str(int(b))), color=plt.cm.gist_rainbow(i))

	plt.legend(loc=0, fontsize='x-large', ncol=2)
	plt.xlabel('frame', fontsize='x-large')
	plt.ylabel('mean neighbours', fontsize='x-large')
	ax = fig.add_subplot(1, 1, 1)
	ax.tick_params(axis='both', labelsize='x-large')
	plt.xlim(-20,max_time+20)
	
	figname = directory + '/mean_neighbours_by_distance_' + str(n) + '_bins_' + kind + '_1.png'
	count = 1
	while os.path.isfile(figname): #if the csv name already exists, make new files with _1, _2, _3 at the end to differentiate
    		count += 1
    		figname = directory + '/mean_neighbours_by_distance_' + str(n) + '_bins_' + kind + '_' + str(count) + '.png'
    	fig.savefig(figname, bbox_inches='tight') #bbox_inches prevents the legend from being cut off
	print figname

print "Directory produced:"
print directory
print ' '
print "Figures produced:"
make_figure(alltimes_gg, 'gg')
make_figure(alltimes_gr, 'gr')
make_figure(alltimes_rr, 'rr')
make_figure(alltimes_rg, 'rg')

# Create a figure with n subplots sharing the x axis
#plt.figure()
#f, axarr = plt.subplots(n, 1, sharex=True)

#for b in (np.array(range(n)) + 1):
 #   axarr[b - 1].scatter(frames, alltimes[(3*b - 3),:], s=5, label='green-green', color='c')
  #  axarr[b - 1].scatter(frames, alltimes[(3*b - 2),:], s=5, label='red-red', color='r')
   # axarr[b - 1].scatter(frames, alltimes[(3*b - 1),:], s=5, label='heterotypic', color='k')
    #axarr[b - 1].set_xticks(np.arange(0, (max_time + 60), 50))
    #axarr[b - 1].set_yticks(np.arange(0, 2.0, 0.5))
    #axarr[b - 1].set_xlim(-10, (max_time + 50))
    #axarr[b - 1].set_ylim(-0.5, 2.0)
    #plt.setp(axarr[b - 1].get_xticklabels(), rotation=30, horizontalalignment='right')
    #plt.setp(axarr[b - 1].get_yticklabels(), rotation=30)
    #axarr[b - 1].set_ylabel(("Bin " + str(b)),rotation=0, fontweight='bold')
    #axarr[b - 1].yaxis.set_label_position("right")

#axarr[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#f.subplots_adjust(hspace=0)
#f.text(0.5, 0.04, 'frame #', ha='center', va='center')
#f.text(0.0, 0.5, 'mean neighbours per cell', rotation=90, ha='center', va='center')
#plt.legend(loc='lower left')




