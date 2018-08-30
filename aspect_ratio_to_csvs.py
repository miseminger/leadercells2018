#!/usr/bin/env python

"""
Author: Madeline Iseminger
Last Modified Date: 21 April 2018

Read in (eg.) neighbours_2_green.csv, from my neighbour-finding code

Outputs:
  - max distance between two points (defines primary axis)
  - max perpendicular distance from primary axis in either direction
  - aspect ratio, always >= 1
  - orientation (in radians counterclockwise from positive x-axis)
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os
import sys
import pandas as pd
from optparse import OptionParser
from sklearn.cluster import DBSCAN

from datetime import datetime as dt
start_time = dt.now()
time_mark = start_time

## OPTIONS PARSING
usage = '%prog options\n  > both files (-g, -r) are required\n  > exactly one radius parameter (-m/-p) must be specified\n'
parser = OptionParser(usage=usage)
parser.add_option('-g', type='string', dest='greenfile', help='path to green CSV file')
parser.add_option('-r', type='string', dest='redfile', help='path to red CSV file')
parser.add_option('-m', type='string', dest='radius_m', help='radius in microns')
parser.add_option('-p', type='string', dest='radius_p', help='radius in pixels')

(options, args) = parser.parse_args()
error = 0
microns_per_pixel = 0.8
if options.radius_p:
    radius = float(options.radius_p) * microns_per_pixel
elif options.radius_m:
    radius = float(options.radius_m)
else:
    error = 1
    print 'Error: no radius specified, use -m or -p'
if not(options.greenfile or options.redfile):
    error = 1
    print 'Error: no file input, use -g and/or -r'
    
## SETUP
print 'Setting up...'
    
# Time functions
def timestring(time):
    time = str(time).split(':')
    seconds = '{:6.2f}'.format(round(float(time[2]), 2)) + 's'
    minutes = '{:>4}'.format(time[1].lstrip('0') + 'm')
    hours = '{:>3}'.format(time[0] + 'h')
    if float(time[0]) == 0:
        if float(time[1]) == 0:
            time = seconds
        else:
            time = minutes + seconds
    else:
        time = hours + minutes + seconds
    return '{:>14}'.format(time)

def mark():
    now = dt.now()
    total = timestring(now - start_time)
    last_block = timestring(now - time_mark)
    print '{:20}'.format(str(time)) + '{:20}'.format(str(total)) + '{:20}'.format(str(last_block))
    return now

def np_read_csv(file):
    da = np.genfromtxt(file, delimiter=',', names=True)  #make an array out of the .csv, called da. max_rows=5 is good for testing!
    return da.view((float, len(da.dtype.names))) #transform everything in the array to a float

try:
    greenfile, redfile = np_read_csv(options.greenfile), np_read_csv(options.redfile)
except KeyError:
    sys.exit('Error: could not parse CSV, use -C if file is for CellProfiler (MATLAB files are default)\n'
             'Use -h for options usage help')

data = np.concatenate((greenfile, redfile), axis=0)

def aspect_ratio(*points):
    
    def _distance(p, q):
        sum = 0
        i = 0
        while i < len(p):
            sum += (p[i] - q[i])**2
            i += 1
        return math.sqrt(sum)
    
    def _rotate(data, angle):
        data = np.transpose(data)
        tmx = np.matrix([[math.cos(angle),-math.sin(angle),0],
                     [math.sin(angle),math.cos(angle),0],
                     [0,0,1]])
        res = np.zeros([2, data.shape[1]])
        for i in range(0, data.shape[1]):
            omx = tmx * np.matrix([[data[0, i], ], [data[1, i], ], [1, ]])
            res[0, i] = omx[0, 0]
            res[1, i] = omx[1, 0]
        return np.transpose(res)
    
    def _line(p1, p2):
        ax = plt.gca()
        xmin, xmax = ax.get_xbound()

        if(p2[0] == p1[0]):
            xmin = xmax = p1[0]
            ymin, ymax = ax.get_ybound()
        else:
            ymax = p1[1]+(p2[1]-p1[1])/(p2[0]-p1[0])*(xmax-p1[0])
            ymin = p1[1]+(p2[1]-p1[1])/(p2[0]-p1[0])*(xmin-p1[0])

        l = mlines.Line2D([xmin,xmax], [ymin,ymax])
        ax.add_line(l)
        return l
    
    max_couple = [None, None]
    max_dist = 0
    max_indices = [None, None]
    
    pi = 0
    while pi < len(points):
        ni = pi + 1
        while ni < len(points):
            
            dist = _distance(points[pi], points[ni])
            if dist > max_dist:
                max_dist = dist
                max_couple = points[pi], points[ni]
                max_indices = [pi, ni]
            
            ni += 1
        pi += 1
        
    angle = math.atan((max_couple[0][1] - max_couple[1][1]) / (max_couple[0][0] - max_couple[1][0]))
    points = _rotate(points, -angle)
    
    max_couple = [points[i] for i in max_indices]
    
    maxh, minh = max_couple[0][1], max_couple[0][1]
    
    for p in points:
        if p[1] > maxh:
            maxh = p[1]
        if p[1] < minh:
            minh = p[1]
    
    min_dist = maxh - minh
    
    return np.array([max_dist, min_dist, (max_dist / min_dist), angle])

redtimeslist = redfile[:,0].tolist()  #list format of the red frame numbers (as many of each frame # as there are detected cells in it)
greentimeslist = greenfile[:,0].tolist()  #list format of the green frame numbers (as many of each frame # as there are detected cells in it)
max_time = greenfile[-1,0]
        
def timeblock(time):
    
    time_mark = dt.now()

    if time > 0 and time % 10 == 0:
        now = dt.now()
        total = timestring(now - start_time)
        last_block = timestring(now - time_mark)
        print '{:20}'.format(str(time)) + '{:20}'.format(str(total)) + '{:20}'.format(str(last_block))
        time_mark =  dt.now()

    redfirsttime = redtimeslist.index(time)  #index for first instance of frame number
    redlasttime = len(redtimeslist) - redtimeslist[::-1].index(time) - 1 #index for last instance of frame number
    
    greenfirsttime = greentimeslist.index(time)  #index for first instance of frame number
    greenlasttime = len(greentimeslist) - greentimeslist[::-1].index(time) - 1 #index for last instance of frame number

    greenids = greenfile[greenfirsttime:(greenlasttime + 1),1] #this is a list of the green cell ids from this time frame
    greenpositions = (greenfile[greenfirsttime:(greenlasttime + 1),4:6]) * microns_per_pixel #x and y positions of green cells in um
    ug, uniquegreenid_indices = np.unique(greenids, return_index=True) #take all the unique green cell ids
    cell_positions_green = greenpositions[uniquegreenid_indices,:] #and the positions associated with them
    
    redids = redfile[redfirsttime:(redlasttime + 1),1] #this is a list of the red cell ids from this time frame
    redpositions = (redfile[redfirsttime:(redlasttime + 1),4:6]) * microns_per_pixel #x and y positions of red cells in um
    ur, uniqueredid_indices = np.unique(redids, return_index=True) #take all the unique red cell ids
    cell_positions_red = redpositions[uniqueredid_indices,:] #and the positions associated with them
 
    cell_positions_total = np.concatenate((cell_positions_green, cell_positions_red), axis=0) #put the unique positions in one big array with green on top and red below
  
    #greentime and redtime are the chunks of the initial arrays for time t
    #heterotypic neighbour IDs are in column 4, homotypic neighbour IDs are in column 3
    greentime = greenfile[greenfirsttime:(greenlasttime + 1),:] #chunk of green array for time t
    redtime = redfile[redfirsttime:(redlasttime + 1),:] #chunk of red array for time t
    
    #arrays containing the id column of greentime and redtime
    greencellids = greentime[:,1].tolist()
    redcellids = redtime[:,1].tolist()
    
    #make empty lists for the neighbours
    homo_neighbours_per_green_cell = []
    hetero_neighbours_per_green_cell = []
    homo_neighbours_per_red_cell = []
    hetero_neighbours_per_red_cell = []

    #now within that time, go through and add the number of homo and hetero neighbours of each individual cell to 
    for i in np.arange(redcellids[0],redcellids[-1] + 1):
        redfirstid = redcellids.index(i)  #index for first instance of cell id i
        redlastid = len(redcellids) - redcellids[::-1].index(i) - 1 #index for last instance of cell id i
    #add up the number of neighbours (homotypic + heterotypic) for cell i and put them in separate lists
        homo_neighbours_per_red_cell.append(np.sum((redtime[redfirstid:(redlastid + 1),2] > 0)*1))
        hetero_neighbours_per_red_cell.append(np.sum((redtime[redfirstid:(redlastid + 1),3] > 0)*1))
   
    #repeat this whole thing for green cells
    for i in np.arange(greencellids[0],greencellids[-1] + 1):
        greenfirstid = greencellids.index(i)
        greenlastid = len(greencellids) - greencellids[::-1].index(i) - 1 
        homo_neighbours_per_green_cell.append(np.sum((greentime[greenfirstid:(greenlastid + 1),2] > 0)*1)) 
        hetero_neighbours_per_green_cell.append(np.sum((greentime[greenfirstid:(greenlastid + 1),3] > 0)*1))
    
    #now we want to be able to add up the neighbours for the cells in each cluster later on, so
    allhetero = np.array(hetero_neighbours_per_green_cell + hetero_neighbours_per_red_cell) #array of all hetero neighbours in time frame, by cell id: green first, red second
    allhomo = np.array(homo_neighbours_per_green_cell + homo_neighbours_per_red_cell) #array of all homo neighbours in time frame, by cell id: green first, red second
    
    #now combine the lists so you get one big list of all the numbers of neighbours had by every cell in the frame
    neighbours_per_red_cell = [x + y for x, y in zip(homo_neighbours_per_red_cell, hetero_neighbours_per_red_cell)] #add up homo and hetero
    neighbours_per_green_cell = [x + y for x, y in zip(homo_neighbours_per_green_cell, hetero_neighbours_per_green_cell)] #add up homo and hetero
    allneighbours = neighbours_per_red_cell + neighbours_per_green_cell #add up red and green lists
    #and find the mean number of neighbours possessed by any given cell in this frame
    mean_total_cell_neighbours_all = np.mean(allneighbours)
                         
    # Clustering
    min_pts = int(math.ceil(mean_total_cell_neighbours_all)) + 1
    db = DBSCAN(eps=radius, min_samples=min_pts).fit(cell_positions_total)
    labels = db.labels_     #Noisy samples are given the label -1

    framedata = np.zeros(((np.max(labels) + 1), 17))
    framedata[:,0] = time #frame number
    framedata[:,1] = np.sum(homo_neighbours_per_red_cell) #number of red-red neighbours in frame
    framedata[:,2] = np.sum(homo_neighbours_per_green_cell) #number of green-green neighbours in frame
    framedata[:,3] = np.sum(hetero_neighbours_per_red_cell + hetero_neighbours_per_green_cell) #number of heterotypic neighbours in frame
    framedata[:,4] = np.max(labels) #number of clusters in frame

    #get the points of the cells in each cluster
    for m in range(np.max(labels)): #start at m=0
        framedata[m,5] = m + 1 # cluster ID--but make them start at 1 instead of 0 for aesthetic reasons :)
        cluster_indices = np.where(labels == m)[0] #locations in array of all cells belonging to cluster m

        #find the number of cells in each cluster that are red or green
        framedata[m,12] = sum(cluster_indices <= (cell_positions_green.shape[0] - 1)) #number of green cells in cluster m
        framedata[m,13] = sum(cluster_indices > (cell_positions_green.shape[0] - 1)) #number of red cells in cluster m

        #find the number of neighbours in each cluster that are red-red, red-green and green-green
	framedata[m,14] = np.sum(allhomo[cluster_indices[cluster_indices > (cell_positions_green.shape[0] - 1)]]) #red-red neighbours in cluster m
	framedata[m,15] = np.sum(allhomo[cluster_indices[cluster_indices <= (cell_positions_green.shape[0] - 1)]]) #green-green neighbours in cluster m
	framedata[m,16] = np.sum(allhetero[cluster_indices]) #red-green neighbours in cluster m; note that I didn't divide by 2 here

        cluster0 = []
        
        transposedpoints = np.transpose(cell_positions_total) #cell_positions_total is already set up for just the time
        
        for ci in cluster_indices:
            cluster0.append([dim[ci] for dim in transposedpoints])
        
        if len(cluster0) == 1:  #aspect ratios can't be calculated for a cluster of 1, so skip it and leave zeros in the row
            continue
        
        else:
            framedata[m,6] = np.mean([item[0] for item in cluster0]) / microns_per_pixel # mean x value of cluster m (px)
            framedata[m,7] = np.mean([item[1] for item in cluster0]) / microns_per_pixel # mean y value of cluster m (px)
            framedata[m,8:12] = aspect_ratio(*cluster0)
        
    #in case any rows were missed (and contain zeros where the aspect ratio should be), here's how to clean up the array:
    rows_to_delete = np.where(framedata[:,6] == 0)[0]    #identify rows missed
    rows_to_delete = np.flip(rows_to_delete, 0) 

    for row in rows_to_delete:      #delete the missed rows from framedata
        framedata = np.delete(framedata, row, 0)

    framedata[:,5] = np.arange(framedata.shape[0]) + 1  #reorder the cluster numbers to go continuously from 1 to n

    return framedata


print 'Processing frames...'
print '{:20}'.format('Frames processed') + '{:20}'.format('Total runtime') + '{:20}'.format('Block runtime')    
alltimeblocks = np.concatenate([timeblock(x) for x in np.arange(0,(max_time + 1))])

columnlabels = ['Metadata_FrameNumber', 'Total Red-Red Neighbours', 'Total Green-Green Neighbours', 'Total Red-Green Neighbours', 'Total Clusters', 'Cluster ID','Mean_X (px)','Mean_Y (px)','Max_Dist (um)', 'Min_Dist (um)', 'Max_Dist / Min_Dist', 'Rotation Angle (rad)', 'Green Cells in Cluster', 'Red Cells in Cluster', 'Total Red-Red Neighbours in Cluster', 'Total Green-Green Neighbours in Cluster', 'Total Red-Green Neighbours in Cluster']
savearray = pd.DataFrame(alltimeblocks, columns=columnlabels)

experiment_index = options.greenfile.find('EXP') 
experiment = options.greenfile[experiment_index:experiment_index+5] #experiment name to include in output figname

underscores = [pos for pos, char in enumerate(options.greenfile) if char == '_']
s = underscores[2]
e = underscores[7]
detail = options.greenfile[s+1:e] #max angle and radius to include in figname

csv_name = experiment + '_clusteringdata_' + detail + '_1.csv'
count = 1
while os.path.isfile(csv_name): #if the csv name already exists, make new files with _1, _2, _3 at the end to differentiate
    count += 1
    csv_name = experiment + '_clusteringdata_' + detail + '_' + str(count) + '.csv'

savearray.to_csv(csv_name, index=False, header=True, sep=',')

# Script completion text
print '\n' + str(int(max_time + 1)) + ' frames processed'
print 'CSV produced: ' + csv_name
print 'Total runtime: ' + timestring(dt.now() - start_time) + '\n'
#
#
