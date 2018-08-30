#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 14:31:08 2018

@author: madeline
"""

from datetime import datetime
start_time = datetime.now()

import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image, ImageMath
from optparse import OptionParser
from graphviz import Graph
import subprocess
import shlex

## OPTIONS PARSING
usage = '%prog options\n  > both files (-g, -r) are required\n'
parser = OptionParser(usage=usage)
parser.add_option('-g', type='string', dest='greenfile', help='path to green CSV file')
parser.add_option('-r', type='string', dest='redfile', help='path to red CSV file')
parser.add_option('-o', type='string', dest='overlaydir', help='overlay node-edge graphs on microscopy images, argument is path to background image file directory ')
parser.add_option('-t', type='string', dest='timearray', help='frame numbers to make graphs for; if not set, default is all frame numbers.  Input frame numbers separated by commas, or the start, end (inclusive) and step, separated by colons')
parser.add_option('-m', action='store_true', dest='movie', help='make a movie out of the overlay images (or graphs if the overlay option is not set)')

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

max_time = int(greenfile[-1,0])

if options.overlaydir:  #put background image file names from directory into a pandas series to search later
	overlaydir = options.overlaydir
	onlyfiles = [f for f in os.listdir(overlaydir) if os.path.isfile(os.path.join(overlaydir, f))]
	Myseries = pd.Series(onlyfiles)

if options.timearray:
	if ':' in options.timearray:
		timearray = np.asarray(options.timearray.split(':')).astype(int) #make user-inputted times into an array
		timearray = np.arange(timearray[0], timearray[1] + 1, timearray[2])
	else:
		timearray = np.asarray(options.timearray.split(',')).astype(int) #make user-inputted times into an array
else: 
	timearray = np.asarray(range(max_time+1))

# Build directory to contain outputs
experiment_index = options.greenfile.find('EXP') 
experiment = options.greenfile[experiment_index:experiment_index+5] #experiment name to include in output figname

underscores = [pos for pos, char in enumerate(options.greenfile) if char == '_']
s = underscores[2]
e = underscores[7]
detail = options.greenfile[s+1:e] #max angle and radius to include in figname

directory = experiment + '_nodeedgegraphs_' + detail
directory_i = 0
while os.path.exists(directory):
    directory_i += 1
    directory = experiment + '_nodeedgegraphs_' + detail + '_' + str(directory_i)
os.makedirs(directory)

if options.overlaydir:
	directory2 = experiment + '_nodeedgegraph_overlays_' + detail
	directory_i = 0
	while os.path.exists(directory2):
	    directory_i += 1
	    directory2 = experiment + '_nodeedgegraph_overlays_' + detail + '_' + str(directory_i)
	os.makedirs(directory2)


redtimeslist = redfile[:,0].tolist()  #list format of the red frame numbers (as many of each frame # as there are detected cells in it)
greentimeslist = greenfile[:,0].tolist()  #list format of the green frame numbers (as many of each frame # as there are detected cells in it)



def timecollector(time):

    	timename = '%03d' % time  #add leading zeros to time for filename search
	
	redfirsttime = redtimeslist.index(time)  #index for first instance of frame number
	redlasttime = len(redtimeslist) - redtimeslist[::-1].index(time) - 1 #index for last instance of frame number
	uniqueindex_red = np.unique(redfile[redfirsttime:(redlasttime+1),1], return_index=True)[1] #this is an array of the unique red indices

	greenfirsttime = greentimeslist.index(time)  #index for first instance of frame number
	greenlasttime = len(greentimeslist) - greentimeslist[::-1].index(time) - 1 #index for last instance of frame number
	uniqueindex_green = np.unique(greenfile[greenfirsttime:(greenlasttime+1),1], return_index=True)[1] #this is an array of the unique green indices
	
	rednames = np.array(['R' + str(int(x)) for x in redfile[redfirsttime:(redlasttime+1),1][uniqueindex_red]])	#red cell names (R1, R2, ...)
	greennames = np.array(['G' + str(int(x)) for x in greenfile[greenfirsttime:(greenlasttime+1),1][uniqueindex_green]])  #green cell names

	redpos = redfile[redfirsttime:(redlasttime+1),4:6][uniqueindex_red]	#array of red cell positions (x and y, in px)
	greenpos = greenfile[greenfirsttime:(greenlasttime+1),4:6][uniqueindex_green]	#array of green cell positions (x and y, in px)

	rrneighbours = []
	rgneighbours = []
	ggneighbours = []
	grneighbours = []


	#now go through the red cells and add their neighbours to the initialized lists
	for m in (np.array(range(int(redfile[(redlasttime),1]))) + 1): #iterate through red cell ids in timeframe, beginning with 1
	        cell_indices = np.where(redfile[redfirsttime:(redlasttime+1),1] == m)[0] #find all the indices of that cell id
		rrn = (redfile[redfirsttime:(redlasttime+1),2][cell_indices]) #this is the homotypic neighbours column for red cell m
		rrn = rrn[np.nonzero(rrn)] #remove zero elements
		rrneighbours.append(np.array(['R' + str(int(x)) for x in rrn]))

		rgn = (redfile[redfirsttime:(redlasttime+1),3][cell_indices]) #this is the heterotypic neighbours column for red cell m
		rgn = rgn[np.nonzero(rgn)] #remove zero elements
		rgneighbours.append(np.array(['G' + str(int(x)) for x in rgn]))



	#now do the same for the green cells
	for m in (np.array(range(int(greenfile[(greenlasttime),1]))) + 1): #iterate through green cell ids in timeframe, beginning with 1
	        cell_indices = np.where(greenfile[greenfirsttime:(greenlasttime+1),1] == m)[0] #find all the indices of that cell id

		ggn = (greenfile[greenfirsttime:(greenlasttime+1),2][cell_indices]) #this is the homotypic neighbours column for green cell m
		ggn = ggn[np.nonzero(ggn)] #remove zero elements
		ggneighbours.append(np.array(['G' + str(int(x)) for x in ggn]))

		grn = (greenfile[greenfirsttime:(greenlasttime+1),3][cell_indices]) #this is the heterotypic neighbours column for green cell m
		grn = grn[np.nonzero(grn)] #remove zero elements
		grneighbours.append(np.array(['R' + str(int(x)) for x in grn]))

	# Generate graph image
	scale = 0.05 

	redf = '#FF9090' #red
	redl = '#660000' #dark red
	greenf = '#5EECFD'  #cyan
	greenl = '#008080' #teal

	g = Graph(name='cellgraph', format='png', engine='neato', strict=True)
	g.attr(size='8,6')

	g.attr('node', pin='true', shape='circle', width='0.5', fixedsize='true', style='filled', penwidth='3',
       colorscheme='set312', fontname='courier', fontsize='10.0')
	g.attr('edge', fontsize='10.0', penwidth='2', fontname='courier')


#first put in all the nodes
    	for i in range(int(redpos.shape[0])):
    	#go through each red cell
    		x, y = redpos[i, 0], redpos[i, 1] #get x and y positions of cells
    		position = str(x * scale) + ',' + str(y * scale * -1)

    		note = ''

    		g.node(rednames[i], label=rednames[i] + note, pos=position, color=redf, fillcolor=redf)


    	for i in range(int(greenpos.shape[0])):
    	#go through each green cell
    		x, y = greenpos[i, 0], greenpos[i, 1] #get x and y positions of cells
    		position = str(x * scale) + ',' + str(y * scale * -1)

    		note = ''

    		g.node(greennames[i], label=greennames[i] + note, pos=position, color=greenf, fillcolor=greenf)


	#then add the edges
    	for i in range(int(redpos.shape[0])):
    	#go through each red cell

    		#then go through the red neighbours of that cell and make red edges
    		for ni in rrneighbours[i]:  
        		lcolor = redl
        		label = ''
        		g.edge(rednames[i], ni, color=redf, label=label)
    
        
    		#then go through the green neighbours of that cell and make black edges
    		for ni in rgneighbours[i]: 
        		lcolor = 'black'
        		label = ''
        		g.edge(rednames[i], ni, color=lcolor, label=label)


    	for i in range(int(greenpos.shape[0])):
   	 #go through each green cell
	
   	 	#then go through the green neighbours of that cell and make green edges
   	 	for ni in ggneighbours[i]:  
     		   	lcolor = greenl
     		   	label = ''
        		g.edge(greennames[i], ni, color=lcolor, label=label)
    	
        
    		#then go through the red neighbours of that cell and make black edges
    		for ni in grneighbours[i]: 
        		lcolor = 'black'
        		label = ''
        		g.edge(greennames[i], ni, color=lcolor, label=label)


	g.render(directory + '/cellgraph_t' + timename,  cleanup=True)

	#stick the graph on a white rectangle to make all the graphs the same dimensions for movie-making
	graph = Image.open(directory + '/cellgraph_t' + timename + '.png')
	width, height = graph.size
	dim = graph.mode

	new_img = Image.new(dim, (int(height * 4. / 3), height), "white")
	new_img.paste(graph,(0,0))

	new_img.save((directory + '/cellgraph_t' + timename + '.png'),"PNG") #replace original .png with the fixed width one

	#now overlay the graph onto its background image if called for
	if options.overlaydir:

		#pad time integer with zeros
		finder = Myseries.str.contains(('t' + timename), regex=False)
		background = overlaydir + '/' + Myseries[finder.index[finder][0]] #filepath for background image for the current frame
	       	
		background = Image.open(background, 'r') #load background image
		background = background.convert("RGBA")	#convert to RGBA

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

    		new_img.save(directory2 + '/cellgraph_overlay_t' + timename + '.png')	#save overlay image


for t in timearray:
	timecollector(t) 

#make a movie out of the saved images
#movie will save in the directory of overlay or graph images
if options.movie:

	if options.overlaydir:

		os.chdir(directory2) #switch to directory2
		firstname = 'cellgraph_overlay_t'
		moviename = 'overlaymovie'
		#create movie using Ceylin's script
		string  = '/zfs/scratch/keshetlab/RoskelleyLab/Scripts/make_movie.sh ' + str(timearray[0]) + ' ' + firstname + ' .png ' + moviename
		subprocess.call(shlex.split(string))
		#resize movie to 360p
		resizestring = 'ffmpeg -i ' + moviename + '-' + str(timearray[0]) + '.mp4 -vf scale=720:-2 -crf 18 ' + moviename + '-' + str(timearray[0]) + '_720p.mp4'
		#delete large copy of movie
		delstring = 'rm ' + moviename + '-' + str(timearray[0]) + '.mp4'
		subprocess.call(shlex.split(delstring))

	else:
		os.chdir(directory) #switch to directory
		firstname = 'cellgraph_t'
		moviename = 'graphmovie'
		#create movie using Ceylin's script
		string  = '/zfs/scratch/keshetlab/RoskelleyLab/Scripts/make_movie.sh ' + str(timearray[0]) + ' ' + firstname + ' .png ' + moviename
		subprocess.call(shlex.split(string))
		#resize movie to 360p
		resizestring = 'ffmpeg -i ' + moviename + '-' + str(timearray[0]) + '.mp4 -vf scale=720:-2 -crf 18 ' + moviename + '-' + str(timearray[0]) + '_720p.mp4'
		subprocess.call(shlex.split(resizestring))
		#delete large copy of movie
		delstring = 'rm ' + moviename + '-' + str(timearray[0]) + '.mp4'
		subprocess.call(shlex.split(delstring))


# Script completion text
print str(timearray.shape[0]) + ' frames processed'
if options.overlaydir:
	print 'Directories produced: ' + str(directory) + ' and ' + str(directory2)
	if options.movie:
		print 'Movie produced: ' + str(directory2) + '/' + str(moviename) + '-' + str(timearray[0]) + '_720p.mp4'
else:
	print 'Directory produced: ' + str(directory)
	if options.movie:
		print 'Movie produced: ' + str(directory) + '/' + str(moviename) + '-' + str(timearray[0]) + '_720p.mp4'

