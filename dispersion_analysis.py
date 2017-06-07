"""
Written by Daniel Green
Created on Tue Jun 06 10:59:26 2017

Purpose:  Process and analyze image data from CSDX experiment.
"""
###############################################################################
import dispersion_plots as disp
import sys
"""
import numpy
import mytools
import prime
import pycine
import scipy.ndimage
import time
"""
###############################################################################
def main():
    # determine filename
    """
    vidfile = raw_input('Enter name of video file:  ')
    """
    if len(sys.argv) == 2:      # check if filename was included in cmdline
        vidfile = sys.argv[1]
    elif len(sys.argv) == 1:    # otherwise request filename
        vidfile = input('Enter name of video file:  ')
    else:
        print('ERROR:  Too many command line arguments.')
        sys.exit()
        
    # read in movie file
    t, images = disp.read_movie(vidfile, framelimits=(0,10000))
    
    # average blocks to determine 2D FFT spectral estimate
    fHz, kpix, power = disp.FFT_map_2D(t, images, center=(63.5,63.5), df=500)
    
    # cut off front and back of vidfile name
    front = vidfile.find('/')
    end = vidfile.find('f0t')
    pref = "CSDX_Plots_frame0-10k_df500/" + vidfile[front+1:end]
    
    # plot the data from 2D FFT dispersion estimate at r = .5,1,1.5,2,... cm
    # save dispersion plot as jpg
    for i in range(1,13):
        print i*.5,'cm'
        ax, cb, im = disp.plot_FFT_2D_dispersion(fHz, kpix, power,
                                                 radius=i*.5, angular=True,
                                                 fileprefix=pref)
        
main()