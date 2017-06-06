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
    t, images = disp.read_movie(vidfile)
    
    # average blocks to determine 2D FFT spectral estimate
    fHz, kpix, power = disp.FFT_map_2D(t, images, center=(63.5,63.5))
    
    # plot the data from 2D FFT dispersion estimate
    ax, cb, im = disp.plot_FFT_2D_dispersion(fHz, kpix, power,
                                             radius=disp.image_rcal()*10)
       
    
    
main()