"""
Written by Daniel Green
Created on Thu Jun 29 16:24:02 2017

Purpose:  Identify and plot max frequency values at each wave number.
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
    # list of all video files
    files = ["B390_1200G_f0t10000.h5", "B450_1400G_f0t10000.h5",
             "B520_1600G_f0t10000.h5", "B580_1800G_f0t10000.h5",
             "B610_1900G_f0t10000.h5", "B650_2000G_f0t10000.h5"]
    
    # iterate through video files (only first 3; fluctation most pronounced)
    for k in range(4):
    #for k in range(3):
        vidfile = "CSDX_vids/20151104_" + files[k]

        # read in movie file
        t, images = disp.read_movie(vidfile, framelimits=(0,10000))
        #t, images = disp.read_movie(vidfile, framelimits=(0,10000))
        
        # determine center of image with center of mass
        xcom, ycom = disp.get_imaging_center(images)
        print "CoM: ", xcom, ycom
        
        # convert frames to polar
        p_images,nr,ntheta = disp.FFT_polar_conv(images,center=(xcom,ycom))
        
        dfreq = [200,500,1000]
        for j in dfreq:
            # average blocks to determine 2D FFT spectral estimate
            fHz, kpix, power = disp.FFT_map_2D(t, p_images, nr, ntheta, df=j)
            
            # cut off front and back of vidfile name
            front = vidfile.find('/')
            end = vidfile.find('f0t')
            pref = "CSDXplots_maxf_4pts_peaks/df" + str(j) + vidfile[front:end]
    
            # plot the data from 2D FFT dispersion estimate at r = .6,.7,...,2.0 cm
            # save dispersion plot as jpg
            for i in range(6,21):
            #for i in range(8,9):
                print i*.1,'cm'
                ax,cb,im = disp.plot_FFT_2D_dispersion(fHz, kpix, power, kmax=1000,
                                                       fmax=75e3, radius=i*.1,
                                                       angular=True, numpts=4,
                                                       filepref=pref)
        
    
main()