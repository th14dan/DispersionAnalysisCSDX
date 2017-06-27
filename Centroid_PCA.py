"""
Written by Daniel Green
Created on Mon Jun 26 15:02:59 2017

Purpose:  Identify centroid and principle component of disperstion plots.
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
    #for k in range(1):
    for k in range(4):
        vidfile = "CSDX_vids/20151104_" + files[k]

        # read in movie file
        #t, images = disp.read_movie(vidfile, framelimits=(0,1000))
        t, images = disp.read_movie(vidfile, framelimits=(0,10000))
        
        # determine center of image with center of mass
        xcom, ycom = disp.get_imaging_center(images)
        print "CoM: ", xcom, ycom
        
        # convert frames to polar
        p_images,nr,ntheta = disp.FFT_polar_conv(images,center=(xcom,ycom))
        
        # average blocks to determine 2D FFT spectral estimate
        fHz, kpix, power = disp.FFT_map_2D(t,p_images,nr,ntheta,df=200)
        
        # cut off front and back of vidfile name
        front = vidfile.find('/')
        end = vidfile.find('f0t')
        pref = "CSDXplots_PCA/df" + str(200) + vidfile[front:end]

        # plot the data from 2D FFT dispersion estimate at r = .5,.6,...,1.5 cm
        # save dispersion plot as jpg
        for i in range(6,21):
            print i*.1,'cm'
            ax,cb,im = disp.plot_FFT_2D_dispersion(fHz, kpix, power, kmax=1000,
                                                   fmax=75e3, radius=i*.1,
                                                   angular=True, pca=True,
                                                   fileprefix=pref)
    


    
main()