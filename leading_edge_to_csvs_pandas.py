#!/usr/bin/python2

"""
Author: Madeline Iseminger
Last Modified Date: 31 Mar 2018

Implements a binary mask technique to find the leading edge of the cell mass in each frame, then increments by 35um along the leading edge 
and matches the position of each increment to the nearest cell found in the CellProfiler output .csv.  The IDs and positions of 
these leading edge cells are recorded in arrays (one for green cells, one for red) and the IDs and positions of the 
homotypic and heterotypic neighbours of the leading cells are then found and recorded in output .csv files.

Use .csv files of the type: AllMyExpt_MyCells_Red.csv, AllMyExpt_MyCells_green.csv, and phase images.

Arguments:
    -g, -r: path to green and red CSV files (at least one file must be provided)
    -m or -p: radius in microns or pixels
Outputs:
    Two CSV files (one for red, one for green) with one row for each cell with frame number, object number, number of green-green neighbours, red-red neighbours, green-red neighbours, and the x and y positions of each of these
    columns of final CSV: 'Metadata_FrameNumber', 'ObjectNumber', 'Homotypic Neighbour IDs', 'Heterotypic Neighbour IDs','Main_Location_Center_X (px)', 'Main_Location_Center_Y (px)','Homo_Location_Center_X (px)', 'Homo_Location_Center_Y (px)','Hetero_Location_Center_X (px)', 'Hetero_Location_Center_Y (px)'
    """
    
#note:  for the CellProfiler data, (0,0) is in the lower left-hand corner of the image (http://forum.cellprofiler.org/t/xy-center-coordinates-for-objects-and-all-their-neighbors/627/3); Python indexes it as the upper lefthand corner.

from datetime import datetime as dt
start_time = dt.now()
time_mark = start_time

import os
import sys
import math
import numpy as np
import pandas as pd
from optparse import OptionParser

from scipy import ndimage
from scipy import misc
#from imageio import imread
import imageio
from scipy.spatial import cKDTree
from skimage.filters.rank import entropy
from skimage.morphology import disk

from graphviz import Graph
import subprocess
import shlex
import matplotlib.pyplot as plt
from PIL import Image, ImageMath

## OPTIONS PARSING
usage = '%prog options\n  > both files (-g, -r) and number of timeframes (-t) are required\n  > exactly one radius parameter (-m/-p) must be specified\n  > -C is an optional flag'
parser = OptionParser(usage=usage)
#parser.add_option('-i', type='string', dest='imgfile', help='path to any of the ch00 phase images')
parser.add_option('-g', type='string', dest='greenfile', help='path to green CSV file')
parser.add_option('-r', type='string', dest='redfile', help='path to red CSV file')
parser.add_option('-m', type='string', dest='radius_m', help='radius in microns')
parser.add_option('-p', type='string', dest='radius_p', help='radius in pixels')
parser.add_option('-C', action='store_true', dest='cp',
                  help='use CellProfiler file instead of MATLAB')
parser.add_option('-t', type='string', dest='timearray', help='frame numbers to make check images for, if desired.  Input frame numbers separated by commas, or the start, end (inclusive) and step, separated by colons')
parser.add_option('-o', type='string', dest='overlaydir', help='overlay node-edge graphs on microscopy images, argument is path to background image file directory ')

(options, args) = parser.parse_args()
error = 0

#if options.cp:
microns_per_pixel = 0.8
start_count = 0

if not(options.greenfile and options.redfile):
    error = 1
    print 'Error: insufficient file input, both files must be provided using -g and -r'

if options.radius_p:
    radius = float(options.radius_p) * microns_per_pixel
elif options.radius_m:
    radius = float(options.radius_m) 
else:
    error = 1
    print 'Error: no radius specified, use -m or -p'
if error:
    sys.exit('Use -h for options usage help')

if options.overlaydir:  #put background image file names from directory into a pandas series to search later
	overlaydir = options.overlaydir
	onlyfiles = [f for f in os.listdir(overlaydir) if os.path.isfile(os.path.join(overlaydir, f))]
	Myseries = pd.Series(onlyfiles)

experiment_index = options.greenfile.find('EXP') 
experiment = options.greenfile[experiment_index:experiment_index+5] #experiment name to include in output figname

underscores = [pos for pos, char in enumerate(options.greenfile) if char == '_']
s = underscores[2]
e = underscores[7]
detail = options.greenfile[s+1:e] #max angle and radius to include in figname

directory = experiment + 'leadingedge_csvs_and_checks_' + detail
directory_i = 1
while os.path.exists(directory):
    directory_i += 1
    directory = experiment + 'leadingedge_csvs_and_checks_' + detail + '_' + str(directory_i)
os.makedirs(directory)

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
    da = np.genfromtxt(file, delimiter=',', names=True) #, usecols=data_cols)  #make an array out of the .csv, called da. max_rows=5 is good for testing!
    return da.view((float, len(da.dtype.names))) #transform everything in the array to a float

try:
    greenfile, redfile = np_read_csv(options.greenfile), np_read_csv(options.redfile)
except KeyError:
    sys.exit('Error: could not parse CSV, use -C if file is for CellProfiler (MATLAB files are default)\n'
             'Use -h for options usage help')

columnlabels = ['Metadata_FrameNumber', 'ObjectNumber', 'Homotypic Neighbour IDs', 'Heterotypic Neighbour IDs','Main_Location_Center_X (px)', 'Main_Location_Center_Y (px)','Homo_Location_Center_X (px)', 'Homo_Location_Center_Y (px)','Hetero_Location_Center_X (px)', 'Hetero_Location_Center_Y (px)']

greendf = pd.DataFrame(data=greenfile, index=greenfile[:,0], columns=columnlabels)
reddf = pd.DataFrame(data=redfile, index=redfile[:,0], columns=columnlabels)

max_time = int(greenfile[-1, 0])    #max_time is the final frame number

if options.timearray:
	if ':' in options.timearray:
		timearray = np.asarray(options.timearray.split(':')).astype(int) #make user-inputted times into an array
		timearray = np.arange(timearray[0], timearray[1] + 1, timearray[2])
	else:
		timearray = np.asarray(options.timearray.split(',')).astype(int) #make user-inputted times into an array
else: 
	timearray = np.asarray(range(max_time+1))
    
# Process frames
time_mark = dt.now()
print 'Setup complete: ' + timestring(time_mark - start_time) + '\n'

total_frames = max_time - start_count + 1 #total_frames is the total number of frames :)



def leadingarray_time():
    
    print 'Finding leading edge cells...'
    print '{:20}'.format('Frames processed') + '{:20}'.format('Total runtime') + '{:20}'.format('Block runtime')
    
    def timeblock(time):

    	#timename = '%03d' % time  #add leading zeros to time for filename search
        timename = '%04d' % time #experiment 9

        time_mark = dt.now()
    
        if time > 0 and time % 10 == 0:
                now = dt.now()
                total = timestring(now - start_time)
                last_block = timestring(now - time_mark)
                print '{:20}'.format(str(time)) + '{:20}'.format(str(total)) + '{:20}'.format(str(last_block))
                time_mark =  dt.now()
                
        #filename = imgfile[:(imgfile[0:-3].rindex('t') + 1)] + str("%03d" % (time,)) + imgfile[imgfile.rindex('_'):]
	finder = Myseries.str.contains(('t' + timename), regex=False)
	filename = overlaydir + '/' + Myseries[finder.index[finder][0]] #filepath for background image for the current frame
        
        greentimeslist = greenfile[:,0].tolist()  #list format of the frame numbers (as many of each frame # as there are cells in it)
        greenfirsttime = greentimeslist.index(time)  #index for first instance of frame number
        greenlasttime = len(greentimeslist) - greentimeslist[::-1].index(time) - 1 #index for last instance of frame number
        
        redtimeslist = redfile[:,0].tolist()  #list format of the frame numbers (as many of each frame # as there are cells in it)
        redfirsttime = redtimeslist.index(time)  #index for first instance of frame number
        redlasttime = len(redtimeslist) - redtimeslist[::-1].index(time) - 1 #index for last instance of frame number
        
        leadingedge_green = []        
        leadingedge_red = []
               
        im = Image.open(filename)
        image_array = np.array(im)
    
        max_y = image_array.shape[0]  #y is in the vertical direction
        
        #apply entropy mask    
        ent = entropy(image_array, disk(5))
        threshold = ent.mean()
        entmask = ent > threshold  #gives boolean array
        entmask = entmask * np.uint8(255) #cell area is mostly 255 (white), background is mostly 0 (black)
        
        #despeckle
        despeckled = ndimage.median_filter(entmask, size=50)  
        
        ##save image if needed for check images
	if time in timearray:
		picname = directory + '/binarymask_t' + ('%03d' % time) + '.png'
		imageio.imwrite(picname, despeckled)

	#get arrays of the unique positions of all the red and green cells in the frame
	allredxy = reddf[reddf['Metadata_FrameNumber']==float(time)][['ObjectNumber','Main_Location_Center_X (px)', 'Main_Location_Center_Y (px)']]
	uniqueredxy = allredxy.drop_duplicates('ObjectNumber')	
	combined_red_x_y_arrays = uniqueredxy.values[:,1:]
	uniqueredids = uniqueredxy.values[:,0]

	allgreenxy = greendf[greendf['Metadata_FrameNumber']==float(time)][['ObjectNumber','Main_Location_Center_X (px)', 'Main_Location_Center_Y (px)']]
	uniquegreenxy = allgreenxy.drop_duplicates('ObjectNumber')
	combined_green_x_y_arrays = uniquegreenxy.values[:,1:]
	uniquegreenids = uniquegreenxy.values[:,0]

        #increment every 35-40um along the y axis
        stepsize = 10 #step size in um
	numpoints = int(max_y/(float(stepsize) / microns_per_pixel))
            
	for y in np.linspace(0,max_y,numpoints,endpoint=False, dtype=int):
            yslice = despeckled[y,:].tolist()  #list all the pixels for a given y, starting with y=0
            try:
                x = yslice[::1].index(0) - 1  #x is the white pixel right before the first black pixel when coming from the left
            except ValueError:
                continue   #if the wound has healed in that row, skip it and find the leading cell in the next row
                
            #may want a bit of code that counts up the rows in which the wound has healed for each frame
            
            # Shoe-horn existing data for entry into KDTree routines
            points_list = [x,y]
        
            def do_kdtree(combined_x_y_arrays,points):
                mytree = cKDTree(combined_x_y_arrays)
                dist, indexes = mytree.query(points, k=1)  #find the closest point
                return dist, indexes
            
            greenresults = do_kdtree(combined_green_x_y_arrays,points_list)
            green_dist = greenresults[0]
            green_indexes = greenresults[1]
            
            redresults = do_kdtree(combined_red_x_y_arrays,points_list)
            red_dist = redresults[0]
            red_indexes = redresults[1]
	                
            #if that position is within one radius of the white pixel, put that cell id in a list of leading edge cells (if it isn't already there)
		#changed to: if red position is closer than green, or vice versa (ignore "one radius" rule)--August 10, 2018
            if red_dist < green_dist and uniqueredids[red_indexes] not in leadingedge_red:
                leadingedge_red.append(uniqueredids[red_indexes])
            elif green_dist < red_dist and uniquegreenids[green_indexes] not in leadingedge_green:
                leadingedge_green.append(uniquegreenids[green_indexes])

	greentimechunk = greendf[greendf['Metadata_FrameNumber']==float(time)]
	redtimechunk = reddf[reddf['Metadata_FrameNumber']==float(time)]

	def catcher(x,timechunk): #returns the chunk of data corresponding to leading cell id x
		leadcell = timechunk[timechunk['ObjectNumber']==float(x)]
		return leadcell

	allleadspertimegreen = pd.concat([catcher(x, greentimechunk) for x in leadingedge_green])
	allleadspertimered = pd.concat([catcher(x, redtimechunk) for x in leadingedge_red])

	return allleadspertimegreen, allleadspertimered

    #stack data for all frames
    allgreens = pd.concat([timeblock(x)[0] for x in np.arange(0,(max_time+1))])
    allreds = pd.concat([timeblock(x)[1] for x in np.arange(0,(max_time+1))])
    return allgreens, allreds

colors = leadingarray_time()
lgdf = colors[0]
lrdf = colors[1]

#rename 'ObjectNumber' column to 'Leading Cell ObjectNumber'
#lgdf = lgdf.rename(index=str, columns={"ObjectNumber": "Leading Cell ObjectNumber"})
#lrdf = lrdf.rename(index=str, columns={"ObjectNumber": "Leading Cell ObjectNumber"})

green_csv_name = directory + '/' + experiment + '_leadingedge_green_' + detail + '_1.csv'
red_csv_name = directory + '/' + experiment + '_leadingedge_red_' + detail + '_1.csv'

count = 1
while os.path.isfile(red_csv_name): #if the csv name already exists, make new files with _1, _2, _3 at the end to differentiate
    count += 1
    red_csv_name = directory + '/' + experiment + '_leadingedge_red_' + detail + '_' + str(count) + '.csv'

count = 1
while os.path.isfile(green_csv_name): #if the csv name already exists, make new files with _1, _2, _3 at the end to differentiate
    count += 1
    green_csv_name = directory + '/' + experiment + '_leadingedge_green_' + detail + '_' + str(count) + '.csv'

lrdf.to_csv(red_csv_name, index=False, header=True, sep=',')
lgdf.to_csv(green_csv_name, index=False, header=True, sep=',')

#now make the check images from the csvs that just got made
print 'CSV files complete!  Making check images...'

def makecheckimages(time):

    	timename = '%04d' % time  #add leading zeros to time for filename search
	
	#get names and positions for all leading and neighbour cells
	redleadids = lrdf[lrdf['Metadata_FrameNumber']==float(time)]['ObjectNumber'].values #red leading cell ids
	redhomoids = lrdf[lrdf['Metadata_FrameNumber']==float(time)]['Homotypic Neighbour IDs'].values #red leading cell ids
	nonzerohomo = np.nonzero(redhomoids)
	redhomoids = redhomoids[nonzerohomo] #remove zeros
	redheteroids = lrdf[lrdf['Metadata_FrameNumber']==float(time)]['Heterotypic Neighbour IDs'].values #red leading cell ids
	nonzerohetero = np.nonzero(redheteroids)
	redheteroids = redheteroids[nonzerohetero] #remove zeros

	leadrednames = np.array(['R' + str(int(x)) for x in redleadids])	#red cell names (R1, R2, ...): there will be doubles in here but it doesn't matter
	redhomonames = np.array(['R' + str(int(x)) for x in redhomoids])
	redheteronames = np.array(['G' + str(int(x)) for x in redheteroids])

	leadredpos = lrdf[lrdf['Metadata_FrameNumber']==float(time)][['Main_Location_Center_X (px)', 'Main_Location_Center_Y (px)']]
	leadredpos = leadredpos.drop_duplicates().values
	#array of red cell positions (x and y, in px)
	redhomopos = lrdf[lrdf['Metadata_FrameNumber']==float(time)][['Homo_Location_Center_X (px)', 'Homo_Location_Center_Y (px)']].values
	redhomopos = redhomopos[nonzerohomo]
	redheteropos = lrdf[lrdf['Metadata_FrameNumber']==float(time)][['Hetero_Location_Center_X (px)', 'Hetero_Location_Center_Y (px)']].values
	redheteropos = redheteropos[nonzerohetero]

	greenleadids = lgdf[lgdf['Metadata_FrameNumber']==float(time)]['ObjectNumber'].values #red leading cell ids
	greenhomoids = lgdf[lgdf['Metadata_FrameNumber']==float(time)]['Homotypic Neighbour IDs'].values #red leading cell ids
	nonzerohomog = np.nonzero(greenhomoids)
	greenhomoids = greenhomoids[nonzerohomog] #remove zeros
	greenheteroids = lgdf[lgdf['Metadata_FrameNumber']==float(time)]['Heterotypic Neighbour IDs'].values #red leading cell ids
	nonzeroheterog = np.nonzero(greenheteroids)
	greenheteroids = greenheteroids[nonzeroheterog] #remove zeros

	leadgreennames = np.array(['G' + str(int(x)) for x in greenleadids])	#red cell names (R1, R2, ...): there will be doubles in here but it doesn't matter
	greenhomonames = np.array(['G' + str(int(x)) for x in greenhomoids])
	greenheteronames = np.array(['R' + str(int(x)) for x in greenheteroids])

	leadgreenpos = lgdf[lgdf['Metadata_FrameNumber']==float(time)][['Main_Location_Center_X (px)', 'Main_Location_Center_Y (px)']]#array of red cell positions (x and y, in px)
	leadgreenpos = leadgreenpos.drop_duplicates().values
	greenhomopos = lgdf[lgdf['Metadata_FrameNumber']==float(time)][['Homo_Location_Center_X (px)', 'Homo_Location_Center_Y (px)']].values
	greenheteropos = lgdf[lgdf['Metadata_FrameNumber']==float(time)][['Hetero_Location_Center_X (px)', 'Hetero_Location_Center_Y (px)']].values
	greenhomopos = greenhomopos[nonzerohomog]
	greenheteropos = greenheteropos[nonzeroheterog]


	greenneighbourpos = np.concatenate((greenhomopos,redheteropos), axis=0)
	greenneighbournames = np.concatenate((greenhomonames,redheteronames), axis=0)

	redneighbourpos = np.concatenate((redhomopos,greenheteropos), axis=0)
	redneighbournames = np.concatenate((redhomonames,greenheteronames), axis=0)
	#also need nested lists of neighbour ids so the program knows where to make edges
	rrneighbours = []
	rgneighbours = []
	ggneighbours = []
	grneighbours = []

	greentimechunk = lgdf[lgdf['Metadata_FrameNumber']==float(time)]
	greencellids = np.unique(greenleadids)

	redtimechunk = lrdf[lrdf['Metadata_FrameNumber']==float(time)]
	redcellids = np.unique(redleadids)


	for m in redcellids:
		rrn = redtimechunk[redtimechunk['ObjectNumber']==float(m)]['Homotypic Neighbour IDs'].values
		rrn = rrn[np.nonzero(rrn)] #remove zero elements
		rrneighbours.append(np.array(['R' + str(int(x)) for x in rrn]))

		rgn = redtimechunk[redtimechunk['ObjectNumber']==float(m)]['Heterotypic Neighbour IDs'].values
		rgn = rgn[np.nonzero(rgn)] #remove zero elements
		rgneighbours.append(np.array(['G' + str(int(x)) for x in rgn]))

	for m in greencellids:
		ggn = greentimechunk[greentimechunk['ObjectNumber']==float(m)]['Homotypic Neighbour IDs'].values
		ggn = ggn[np.nonzero(ggn)] #remove zero elements
		ggneighbours.append(np.array(['G' + str(int(x)) for x in ggn]))

		grn = greentimechunk[greentimechunk['ObjectNumber']==float(m)]['Heterotypic Neighbour IDs'].values
		grn = grn[np.nonzero(grn)] #remove zero elements
		grneighbours.append(np.array(['R' + str(int(x)) for x in grn]))

	# Generate graph image
	scale = 0.05 

	redf = '#FF9090' #red
	redl = '#660000' #dark red
	greenf = '#5EECFD'  #cyan
	greenl = '#008080' #teal
	yellow = '#FFFF00' #yellow to outline leading edge cells

	g = Graph(name='cellgraph', format='png', engine='neato', strict=True)
	g.attr(size='8,6')

	g.attr('node', pin='true', shape='circle', width='0.5', fixedsize='true', style='filled', penwidth='3',
       colorscheme='set312', fontname='courier', fontsize='10.0')
	g.attr('edge', fontsize='10.0', penwidth='2', fontname='courier')

#first put in all the nodes
    	for i in range(int(redneighbourpos.shape[0])):
    	#go through each red cell
    		x, y = redneighbourpos[i, 0], redneighbourpos[i, 1] #get x and y positions of cells
    		position = str(x * scale) + ',' + str(y * scale * -1)

    		note = ''

    		g.node(redneighbournames[i], label=redneighbournames[i] + note, pos=position, color=redf, fillcolor=redf)


    	for i in range(int(greenneighbourpos.shape[0])):
    	#go through each green cell
    		x, y = greenneighbourpos[i, 0], greenneighbourpos[i, 1] #get x and y positions of cells
    		position = str(x * scale) + ',' + str(y * scale * -1)

    		note = ''

    		g.node(greenneighbournames[i], label=greenneighbournames[i] + note, pos=position, color=greenf, fillcolor=greenf)


    	for i in range(int(leadredpos.shape[0])):
    	#go through each red cell
    		x, y = leadredpos[i, 0], leadredpos[i, 1] #get x and y positions of cells
    		position = str(x * scale) + ',' + str(y * scale * -1)

    		note = ''

    		g.node(leadrednames[i], label=leadrednames[i] + note, pos=position, color=yellow, fillcolor=redf)


    	for i in range(int(leadgreenpos.shape[0])):
    	#go through each green cell
    		x, y = leadgreenpos[i, 0], leadgreenpos[i, 1] #get x and y positions of cells
    		position = str(x * scale) + ',' + str(y * scale * -1)

    		note = ''

    		g.node(leadgreennames[i], label=leadgreennames[i] + note, pos=position, color=yellow, fillcolor=greenf)


	#then add the edges
    	for i in range(int(leadredpos.shape[0])):
    	#go through each red cell

    		#then go through the red neighbours of that cell and make red edges
    		for ni in rrneighbours[i]:  
        		lcolor = redl
        		label = ''
        		g.edge(leadrednames[i], ni, color=redf, label=label)
    
        
    		#then go through the green neighbours of that cell and make black edges
    		for ni in rgneighbours[i]: 
        		lcolor = 'black'
        		label = ''
        		g.edge(leadrednames[i], ni, color=lcolor, label=label)


    	for i in range(int(leadgreenpos.shape[0])):
   	 #go through each green cell
	
   	 	#then go through the green neighbours of that cell and make green edges
   	 	for ni in ggneighbours[i]:  
     		   	lcolor = greenl
     		   	label = ''
        		g.edge(leadgreennames[i], ni, color=lcolor, label=label)
    	
        
    		#then go through the red neighbours of that cell and make black edges
    		for ni in grneighbours[i]: 
        		lcolor = 'black'
        		label = ''
        		g.edge(leadgreennames[i], ni, color=lcolor, label=label)


	g.render(directory + '/cellgraph_t' + timename,  cleanup=True)

	#stick the graph on a white rectangle to make all the graphs the same dimensions for movie-making
	graph = Image.open(directory + '/cellgraph_t' + timename + '.png')
	width, height = graph.size
	dim = graph.mode

	new_img = Image.new(dim, (int(height * 4. / 3), height), "white")
	new_img.paste(graph,(0,0))

	new_img.save((directory + '/cellgraph_t' + timename + '.png'),"PNG") #replace original .png with the fixed width one


	#now overlay the graph onto its background image

	finder = Myseries.str.contains(('t' + timename), regex=False)
	background = overlaydir + '/' + Myseries[finder.index[finder][0]] #filepath for background image for the current frame
	       	
	background = Image.open(background, 'r') #load background image
	background = background.convert("RGBA")	#convert to RGBA

	binarymask = Image.open((directory + '/binarymask_t' + ('%03d' % time) + '.png'),'r')
	binarymask = binarymask.convert("RGBA")	#convert to RGBA
	background = Image.blend(background, binarymask, 0.5)
	
	graph = Image.open((directory + '/cellgraph_t' + timename + '.png'), 'r') #open resized .png, should be in RGBA format
	graph_rgb = np.array(graph)[:,:,0:3]

	sat_mask = (graph_rgb==255) #saturated r,g, or b values
	white_mask = np.invert(np.sum(sat_mask,axis=-1)==3) #which pixels are not white
	rgba_mask = np.dstack([white_mask]*4) #we want to work in RGBA - broadcasting would be more efficient...
	rgba_mask = Image.fromarray(np.uint8(255*rgba_mask)) #PIL masking requires images, not boolean arrays

	ratio = float(background.size[0]) / background.size[1]
	background = background.resize((int(round(graph.size[1] * ratio)), graph.size[1])) #resize background to fit graph image

	new_img = Image.new('RGBA', np.shape(background)[0:2][::-1], (0, 0, 0, 0)) 
	new_img.paste(background,(0,0)) #put in the background image
	new_img.paste(graph,(0,0),mask=rgba_mask) #put in the overlaying image 

    	new_img.save(directory + '/leadingedge_checkimg_t' + timename + '.png')	#save overlay image



#make check images: node-edge graphs overlaid on background images
if options.timearray:
	[makecheckimages(x) for x in timearray]

# Script completion text
print '\n' + str(int(total_frames)) + ' frames processed'
print 'Directory produced: ' + directory
print 'CSVs produced:'
print green_csv_name 
print red_csv_name
print 'Check images produced for frames ' + str(timearray)
print 'Total runtime: ' + timestring(dt.now() - start_time) + '\n'

