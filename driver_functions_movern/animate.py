#!/usr/bin/env python

import matplotlib.pylab as plt
import matplotlib.animation as anim
import numpy as np
import sys,os

#global frame variables
fig = 0
ax1 = 0
im = 0

#global datasets
time = 0
raster = 0
title = 0

def load_memory(root):
	raster_frames = []
	for i in xrange(get_num_frames(root)):
		raster_frames.append(np.loadtxt("%s%d.asc" % (root, i+1), skiprows=6))
	title, h, data = load_data("."+root+"_frame_metadata")
	time = data[:, h.index("Time")]
	return np.array(raster_frames), time

def get_num_frames(root):
	n = 0
	while 1:
		filename = "%s%d.asc" % (root, n+1)
		if os.path.isfile(filename):
			n += 1
		else:
			return n

def load_data(filename):
	f = open(filename)
	title = f.next()
	h = f.next().strip('\n').split('\t')
	f.close()
	data = np.loadtxt(filename, skiprows=2)

	return title, h, data

def run(root, ui=False):
	global fig
	global ax1
	global im
	global raster,time
	global title

	#setup figure
	fig = plt.figure()
	title = fig.suptitle("Time 0 kyr")
	ax1 = plt.axes()
	ax1.set_xticks([])
	ax1.set_yticks([])

	raster,time = load_memory(root)

	im = ax1.imshow(raster[0,:,:], cmap='gray')
	im.set_clim(0, np.max(raster))

	ani = anim.FuncAnimation(fig, func, frames=get_num_frames(root), interval=30, blit=False)

	if ui:
		while 1:
			save = raw_input("Would you like to save the figure? (Y/N) ")
			if save.lower() in ["y", "yes"]:
				print "Saving as mp4"
				ani.save(root + ".mp4", fps = 10)

				print "Converting to gif"
				print "Isolating frames"
				os.popen("ffmpeg -i %s.mp4 -r 10 .output%%05d.png" % root)
				print "Compiling gif"
				os.popen("convert .output*.png %s.gif" % root)
				print "Cleaning up temporary png files"
				os.popen("rm .output*.png")
				break
			elif save.lower() in ["n", "no"]:
				break
			else:
				print "Invalid response"
	plt.show()
	if ui:
		while 1:
			delete = raw_input("Would you like to clear frames from memory? (Y/N) ")
			if delete.lower() in ["y", "yes", "aye"]:
				os.remove(".%s_frame_metadata" % root)
				i = 1
				while 1:
					if os.path.isfile("%s%d.asc" % (root, i)):
						os.remove("%s%d.asc" % (root, i))
						try:
							os.remove("%s%d_erosion.asc" % (root, i))
						except:
							pass
						try:
							os.remove("%s%d_sa" %(root, i))
						except:
							pass
						i+=1
					else:
						break
				break
			elif delete.lower() in ["n", "no"]:
				break
			else:
				print "Invalid response"

def func(n):
    #fig.suptitle("Time %d kyr" % (time[n]/1000.))
    title.set_text("Time %d kyr" % (time[n]/1000.))
    im.set_data(raster[n,:,:])

    return im


if __name__ == "__main__":
	if len(sys.argv) > 1:
		root = sys.argv[1]
	else:
		root = raw_input("Enter root name: ")
	run(root, ui=True)
