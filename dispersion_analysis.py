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
    if len(sys.argv) == 2:      # check if filename was included in cmdline
        vidfile = sys.argv[1]
    elif len(sys.argv) == 1:    # otherwise request filename
        vidfile = input('Enter name of video file:  ')
    else:
        print('ERROR:  Too many command line arguments.')
        sys.exit()
    """
    files = ["B390_1200G_f0t10000.h5", "B450_1400G_f0t10000.h5",
             "B520_1600G_f0t10000.h5", "B580_1800G_f0t10000.h5",
             "B610_1900G_f0t10000.h5", "B650_2000G_f0t10000.h5"]
    
    for k in range(len(files)):
        vidfile = "CSDX_vids/20151104_" + files[k]
        
        # read in movie file
        t, images = disp.read_movie(vidfile, framelimits=(0,10000))
        
        # convert frames to polar
        p_images,nr,ntheta = disp.FFT_polar_conv(images,center=(63.5,63.5))
        
        delf = [100,200,500,1000]
        for j in delf:
            # average blocks to determine 2D FFT spectral estimate
            fHz, kpix, power = disp.FFT_map_2D(t,p_images,nr,ntheta,df=j)
        
            # cut off front and back of vidfile name
            front = vidfile.find('/')
            end = vidfile.find('f0t')
            pref = "CSDXplots_smallR_negKwide/df" + str(j) + vidfile[front:end]
            
            # plot the data from 2D FFT dispersion estimate at r = .5,1,1.5,2,... cm
            # save dispersion plot as jpg
            for i in range(1,6):
            #for i in range(1,13):
                print i*.5,'cm'
                ax, cb, im = disp.plot_FFT_2D_dispersion(fHz, kpix, power,kmax=1000,
                                                         radius=i*.5, angular=True,
                                                         fileprefix=pref)
            
main()