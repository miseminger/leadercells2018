#!/usr/bin/python2

"""
Author: Madeline Iseminger
Last Modified Date: 31 Mar 2018

Processes CellProfiler or MATLAB .csv files and determines homotypic and heterotypic cell neighbour IDs and positions for each cell at all frames
use .csv files of the type: AllMyExpt_MyCells_Red.csv, AllMyExpt_MyCells_green.csv
Arguments:
    -g, -r: path to green and red CSV files (at least one file must be provided)
    -m or -p: radius in microns or pixels
Outputs:
    Two CSV files (one for red, one for green) with one row for each cell with frame number, object number, number of green-green neighbours, red-red neighbours, green-red neighbours, and the x and y positions of each of these
    columns of final CSV: 'Metadata_FrameNumber', 'ObjectNumber', 'Homotypic Neighbour IDs', 'Heterotypic Neighbour IDs','Main_Location_Center_X (px)', 'Main_Location_Center_Y (px)','Homo_Location_Center_X (px)', 'Homo_Location_Center_Y (px)','Hetero_Location_Center_X (px)', 'Hetero_Location_Center_Y (px)'
    """

from datetime import datetime as dt
start_time = dt.now()
time_mark = start_time

import os
import sys
import math
import numpy as np
import numpy.linalg as la
import pandas as pd
from optparse import OptionParser

## OPTIONS PARSING
usage = '%prog options\n  > both files (-g, -r) and number of timeframes (-t) are required\n  > exactly one radius parameter (-m/-p) must be specified\n  > -C is an optional flag'
parser = OptionParser(usage=usage)
parser.add_option('-g', type='string', dest='greenfile', help='path to green CSV file')
parser.add_option('-r', type='string', dest='redfile', help='path to red CSV file')
parser.add_option('-m', type='string', dest='radius_m', help='radius in microns')
parser.add_option('-p', type='string', dest='radius_p', help='radius in pixels')
parser.add_option('-C', action='store_true', dest='cp',
                  help='use CellProfiler file instead of MATLAB')
parser.add_option('-a', type='string', dest='angle', help='maximum angle for true neighbours in degrees')

(options, args) = parser.parse_args()
error = 0

if options.cp:
    microns_per_pixel = 0.8
    start_count = 0
    data_cols = ['Metadata_FrameNumber', 'ObjectNumber', 'Location_Center_X', 'Location_Center_Y']
else:
    microns_per_pixel = 0.8
    start_count = 1
    data_cols = ['Frame', 'CentroidX', 'CentroidY']

if not(options.greenfile and options.redfile):
    error = 1
    print 'Error: insufficient file input, both files must be provided using -g and -r'

if options.angle:
    max_angle = float(options.angle)
else:
    max_angle = float(0.0)


if options.radius_p:
    radius = float(options.radius_p) * microns_per_pixel
elif options.radius_m:
    radius = float(options.radius_m)
else:
    error = 1
    print 'Error: no radius specified, use -m or -p'
if error:
    sys.exit('Use -h for options usage help')


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


## EXTRACT DATA
def np_read_csv(file):
    da = np.genfromtxt(file, delimiter=',', names=True, usecols=data_cols)  #make an array out of the .csv, called da. max_rows=5 is good for testing!
    return da.view((float, len(da.dtype.names))) #transform everything in the array to a float

try:
    green, red = np_read_csv(options.greenfile), np_read_csv(options.redfile)
except KeyError:
    sys.exit('Error: could not parse CSV, use -C if file is for CellProfiler (MATLAB files are default)\n'
             'Use -h for options usage help')
max_time = int(green[-1, 0])
    
# Process frames
time_mark = dt.now()
print 'Setup complete: ' + timestring(time_mark - start_time) + '\n'

total_frames = max_time - start_count + 1 #total_frames is the total number of frames :)

def py_ang(v1, v2):
	""" Returns the angle in radians between vectors 'v1' and 'v2'    """
	cosang = np.dot(v1, v2)
	sinang = la.norm(np.cross(v1, v2))
	return math.degrees(np.arctan2(sinang, cosang))


#FIND NEIGHBOURS
def find_neighbours(primary, secondary):

    print 'Processing frames...'
    print '{:20}'.format('Frames processed') + '{:20}'.format('Total runtime') + '{:20}'.format('Block runtime')
        
#    time_mark = dt.now()
    
    def timecollector(time):
        time_mark = dt.now()

        if time > 0 and time % 10 == 0:
            now = dt.now()
            total = timestring(now - start_time)
            last_block = timestring(now - time_mark)
            print '{:20}'.format(str(time)) + '{:20}'.format(str(total)) + '{:20}'.format(str(last_block))
            time_mark =  dt.now()

        timeslist = primary[:,0].tolist()  #list format of the frame numbers (as many of each frame # as there are cells in it)
        firsttime = timeslist.index(time)  #index for first instance of frame number
        lasttime = len(timeslist) - timeslist[::-1].index(time) - 1 #index for last instance of frame number
      
        def makecellblock(i):
    
            homo_neighbour_ids = []
            homo_neighbour_x = []
            homo_neighbour_y = []
    
            hetero_neighbour_ids = []
            hetero_neighbour_x = []
            hetero_neighbour_y = []
    
            x, y = primary[i,2] * microns_per_pixel, primary[i,3] * microns_per_pixel

        #now go through and find all the green neighbours of cell i in that same timeframe (these are called ni)
        #find homotypic neighbours
            homo_ni_array = np.arange(firsttime,(lasttime + 1))

            for ni in homo_ni_array: 
                nx, ny = primary[ni,2] * microns_per_pixel, primary[ni,3] * microns_per_pixel
                distance = math.sqrt((x - nx)**2 + (y - ny)**2)
                
		#first collect all ids and positions of cells within the radius into lists   
                if distance < float(radius) and ni != i:
                        ni_id = primary[ni,1]
                        homo_neighbour_ids.append(ni_id)
                        homo_neighbour_x.append(primary[ni,2])
                        homo_neighbour_y.append(primary[ni,3])

	    false_homo_neighbour_ids = []
	    false_homo_neighbour_x = []
	    false_homo_neighbour_y = []

	    for m in range(len(homo_neighbour_ids)): #for each potential neighbour
		m_id = homo_neighbour_ids[m]
		m_x = homo_neighbour_x[m] * microns_per_pixel
		m_y = homo_neighbour_y[m] * microns_per_pixel
		m_distance = math.sqrt((x - m_x)**2 + (y - m_y)**2) #distance neighbouring cell m is from cell i
		vm = (m_x, m_y)
		for j in range(len(homo_neighbour_ids)): #go through all the other neighbouring cells; if the angle is less than max_angle, count it!
		    j_id = homo_neighbour_ids[j]
		    j_x = homo_neighbour_x[j] * microns_per_pixel
		    j_y = homo_neighbour_y[j] * microns_per_pixel
		    j_distance = math.sqrt((x - j_x)**2 + (y - j_y)**2) #distance neighbouring cell m is from cell i
		    vj = (j_x, j_y)

		    angle = py_ang(vj, vm) #returns the angle between the vectors in degrees
		    if angle <= max_angle and j_id != m_id:
			if m_distance > j_distance and (m_id not in false_homo_neighbour_ids): #if m is closer, it is the true neighbour
			    false_homo_neighbour_ids.append(m_id)
			    false_homo_neighbour_x.append(m_x)
			    false_homo_neighbour_y.append(m_y)
			elif m_distance <= j_distance and (j_id not in false_homo_neighbour_ids): #if j is closer, it is the true neighbour
			    false_homo_neighbour_ids.append(j_id)
			    false_homo_neighbour_x.append(j_x)
			    false_homo_neighbour_y.append(j_y)

	    find = lambda searchList, elem: [[i for i, x in enumerate(searchList) if x == e] for e in elem]
	    indices_to_delete = find(homo_neighbour_ids, false_homo_neighbour_ids) #returns a nested list of indices of the false neighbour ids within homo_neighbour_ids
	    indices_to_delete = np.array([item for sublist in indices_to_delete for item in sublist]) #flattens the list

	    homo_neighbour_ids = list(np.delete(np.array(homo_neighbour_ids), indices_to_delete))
	    homo_neighbour_x = list(np.delete(np.array(homo_neighbour_x), indices_to_delete))
	    homo_neighbour_y = list(np.delete(np.array(homo_neighbour_y), indices_to_delete))

            if len(homo_neighbour_ids) == 0:	#if that cell has no homotypic neighbours
                homo_neighbour_ids.append(0)
                homo_neighbour_x.append(0)
                homo_neighbour_y.append(0)

	#now do the exact same thing for the heterotypic neighbours

        #find heterotypic neighbour ids
            secondary_timeslist = secondary[:,0].tolist()
            ni_firsttime = secondary_timeslist.index(time)  #index for first instance of frame number
            ni_lasttime = len(secondary_timeslist) - secondary_timeslist[::-1].index(time) - 1 #index for last instance of frame number
            hetero_ni_array = np.arange(ni_firsttime,(ni_lasttime + 1))

            for ni in hetero_ni_array: 
                nx, ny = secondary[ni,2] * microns_per_pixel, secondary[ni,3] * microns_per_pixel
                distance = math.sqrt((x - nx)**2 + (y - ny)**2)
                    
                if distance < float(radius):
                        ni_id = secondary[ni,1]
                        hetero_neighbour_ids.append(ni_id)
                        hetero_neighbour_x.append(secondary[ni,2])
                        hetero_neighbour_y.append(secondary[ni,3])

#now go through all those and find the false neighbours
	    false_hetero_neighbour_ids = []
	    false_hetero_neighbour_x = []
	    false_hetero_neighbour_y = []

	    for m in range(len(hetero_neighbour_ids)): #for each potential neighbour
		m_id = hetero_neighbour_ids[m]
		m_x = hetero_neighbour_x[m] * microns_per_pixel
		m_y = hetero_neighbour_y[m] * microns_per_pixel
		m_distance = math.sqrt((x - m_x)**2 + (y - m_y)**2) #distance neighbouring cell m is from cell i
		vm = (m_x, m_y)
		for j in range(len(hetero_neighbour_ids)): #go through all the other neighbouring cells; if the angle is less than max_angle, count it!
		    j_id = hetero_neighbour_ids[j]
		    j_x = hetero_neighbour_x[j] * microns_per_pixel
		    j_y = hetero_neighbour_y[j] * microns_per_pixel
		    j_distance = math.sqrt((x - j_x)**2 + (y - j_y)**2) #distance neighbouring cell m is from cell i
		    vj = (j_x, j_y)

		    angle = py_ang(vj, vm) #returns the angle between the vectors in degrees
		    if angle <= max_angle and j_id != m_id:
			if m_distance > j_distance and (m_id not in false_hetero_neighbour_ids): #if j is closer, it is the true neighbour
			    false_hetero_neighbour_ids.append(m_id)
			    false_hetero_neighbour_x.append(m_x)
			    false_hetero_neighbour_y.append(m_y)
			elif m_distance <= j_distance and (j_id not in false_hetero_neighbour_ids): #if m is closer, it is the true neighbour
			    false_hetero_neighbour_ids.append(j_id)
			    false_hetero_neighbour_x.append(j_x)
			    false_hetero_neighbour_y.append(j_y)

	    find = lambda searchList, elem: [[i for i, x in enumerate(searchList) if x == e] for e in elem] #returns list of indices at which to delete
	    indices_to_delete = find(hetero_neighbour_ids, false_hetero_neighbour_ids) #returns a nested list of indices of the false neighbour ids within homo_neighbour_ids
	    indices_to_delete = [item for sublist in indices_to_delete for item in sublist] #flattens the list

	    hetero_neighbour_ids = list(np.delete(np.array(hetero_neighbour_ids), indices_to_delete))
	    hetero_neighbour_x = list(np.delete(np.array(hetero_neighbour_x), indices_to_delete))
	    hetero_neighbour_y = list(np.delete(np.array(hetero_neighbour_y), indices_to_delete))

            if len(hetero_neighbour_ids) == 0:	#if that cell has no homotypic neighbours
                hetero_neighbour_ids.append(0)
                hetero_neighbour_x.append(0)
                hetero_neighbour_y.append(0)

            homo_neighbour_ids = np.transpose(np.asarray(homo_neighbour_ids))
            homo_neighbour_ids.shape = (homo_neighbour_ids.shape[0],1)
            homo_neighbour_x = np.transpose(np.asarray(homo_neighbour_x))
            homo_neighbour_x.shape = (homo_neighbour_x.shape[0],1)
            homo_neighbour_y = np.transpose(np.asarray(homo_neighbour_y))
            homo_neighbour_y.shape = (homo_neighbour_y.shape[0],1)
    
            hetero_neighbour_ids = np.transpose(np.asarray(hetero_neighbour_ids))
            hetero_neighbour_ids.shape = (hetero_neighbour_ids.shape[0],1)
            hetero_neighbour_x = np.transpose(np.asarray(hetero_neighbour_x))
            hetero_neighbour_x.shape = (hetero_neighbour_x.shape[0],1)
            hetero_neighbour_y = np.transpose(np.asarray(hetero_neighbour_y))
            hetero_neighbour_y.shape = (hetero_neighbour_y.shape[0],1)
    
            #now make a block for the first cell ID, primary[i,1]    
            cellblocklength = max(homo_neighbour_ids.shape[0],hetero_neighbour_ids.shape[0])
    
            #make the arrays the same length by adding zeros on the end of the shorter one
            difference = np.abs(homo_neighbour_ids.shape[0] - hetero_neighbour_ids.shape[0])
            addzeros = np.zeros((difference,1))
            if homo_neighbour_ids.shape[0] < hetero_neighbour_ids.shape[0]:
                homo_neighbour_ids = np.concatenate((homo_neighbour_ids, addzeros))
                homo_neighbour_x = np.concatenate((homo_neighbour_x, addzeros))
                homo_neighbour_y = np.concatenate((homo_neighbour_y, addzeros))
            else:
                hetero_neighbour_ids = np.concatenate((hetero_neighbour_ids, addzeros))
                hetero_neighbour_x = np.concatenate((hetero_neighbour_x, addzeros))
                hetero_neighbour_y = np.concatenate((hetero_neighbour_y, addzeros))
            
            cellblock = np.zeros((cellblocklength,2))
            cellblock[:,0] = time  #first column is time
            cellblock[:,1] = primary[i,1]  #second column is cell ID
    
            #after you make all the timeblocks the right length, append the neighbour ids!
            cellblock = np.concatenate((cellblock, homo_neighbour_ids), axis=1)  #third column is homotypic neighbour IDs
            cellblock = np.concatenate((cellblock, hetero_neighbour_ids), axis=1)  #third column is homotypic neighbour IDs
            
            mainpos = np.zeros((cellblocklength,2))
            mainpos[:,0] = primary[i,2]  #first column is x position of main cell in px
            mainpos[:,1] = primary[i,3]  #second column is y position of main cell in px
    
            cellblock = np.concatenate((cellblock, mainpos), axis=1)  #fifth column is x position of main cell, sixth is y position
            cellblock = np.concatenate((cellblock, homo_neighbour_x), axis=1)  #seventh column is x position of homo neighbour
            cellblock = np.concatenate((cellblock, homo_neighbour_y), axis=1)  #eighth column is y position of homo neighbour
            cellblock = np.concatenate((cellblock, hetero_neighbour_x), axis=1)  #ninth column is x position of hetero neighbour
            cellblock = np.concatenate((cellblock, hetero_neighbour_y), axis=1)  #tenth column is y position of hetero neighbour
    
            return cellblock

        allcellblocks = np.concatenate([makecellblock(x) for x in np.arange(firsttime,(lasttime + 1))])
        return allcellblocks

    alltimes = np.concatenate([timecollector(x) for x in np.arange(start_count,(max_time + 1))])
    return alltimes


#now loop through and find the neighbours for everything

#define column labels
columnlabels = ['Metadata_FrameNumber', 'ObjectNumber', 'Homotypic Neighbour IDs', 'Heterotypic Neighbour IDs','Main_Location_Center_X (px)', 'Main_Location_Center_Y (px)','Homo_Location_Center_X (px)', 'Homo_Location_Center_Y (px)','Hetero_Location_Center_X (px)', 'Hetero_Location_Center_Y (px)']

print "Finding neighbours of red cells..."
np_neighbours_red = find_neighbours(red, green)
np_neighbours_red = pd.DataFrame(np_neighbours_red, columns=columnlabels)

print "Finding neighbours of green cells..."
np_neighbours_green = find_neighbours(green, red)    
np_neighbours_green = pd.DataFrame(np_neighbours_green, columns=columnlabels)

#save csvs
experiment_index = options.greenfile.find('EXP') #get experiment name to include in output filenames
experiment = options.greenfile[experiment_index:experiment_index+5]

green_csv_name = experiment + '_allneighbours_green_'  + str(int(radius)) + 'um_radius_' + str(max_angle) + 'deg_max_angle_1.csv'
red_csv_name = experiment + '_allneighbours_red_'  + str(int(radius)) + 'um_radius_' + str(max_angle) + 'deg_max_angle_1.csv'

count = 1
while os.path.isfile(red_csv_name): #if the csv name already exists, make new files with _1, _2, _3 at the end to differentiate
    count += 1
    red_csv_name = experiment + '_allneighbours_red_'  + str(int(radius)) + 'um_radius_' + str(max_angle) + 'deg_max_angle_' + str(count) + '.csv'

count = 1
while os.path.isfile(green_csv_name): #if the csv name already exists, make new files with _1, _2, _3 at the end to differentiate
    count += 1
    green_csv_name = experiment + '_allneighbours_green_'  + str(int(radius)) + 'um_radius_' + str(max_angle) + 'deg_max_angle_' + str(count) + '.csv'

np_neighbours_red.to_csv(red_csv_name, index=False, header=True, sep=',')
np_neighbours_green.to_csv(green_csv_name, index=False, header=True, sep=',')

# Script completion text
print '\n' + str(int(total_frames)) + ' frames processed'
print ' '
print 'CSVs produced:'
print green_csv_name 
print red_csv_name
print ' '
print 'Total runtime: ' + timestring(dt.now() - start_time) + '\n'

