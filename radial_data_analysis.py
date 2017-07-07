"""
Written by Daniel Green
Created on Thu Jul 06 13:16:53 2017

Purpose:  Read in CSDX data from excel file and plot it.
"""
###############################################################################

import pandas
import numpy as np
import matplotlib.pyplot as plt

###############################################################################

def main():
    # identify excel files that will be converted into numpy arrays
    file1 = "CSDX_radial_data/Te_CSDX_2015.xlsx"
    file2 = "CSDX_radial_data/v_ExB_vs_B_from_LIF_OCT_2015.xlsx"
    
    # file1 lists electron temperatures at various radii and B-field strengths
    eTempData = pandas.read_excel(file1)
    eTempArray = eTempData.as_matrix()
    
    # 2x18 arrays of electron temp at various radii
    # each array is for different B-field strengths
    """
    eTemp800 = np.array([eTempArray[:,0], eTempArray[:,1]])
    eTemp1000 = np.array([eTempArray[:,0], eTempArray[:,2]])
    """
    eTemp1200 = np.array([eTempArray[:,0], eTempArray[:,3]])
    eTemp1400 = np.array([eTempArray[:,0], eTempArray[:,4]])
    
    i1200Data = pandas.read_excel(file2, sheetname=3)
    i1200Array = i1200Data.as_matrix()
    i1200Array[-1,1:] = np.array([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
    iDens1200 = np.array([i1200Array[:,0],i1200Array[:,1]])
    iTemp1200 = np.array([i1200Array[:,0],i1200Array[:,2]])

    i1400Data = pandas.read_excel(file2, sheetname=5)
    i1400Array = i1400Data.as_matrix()
    iDens1400 = np.array([i1400Array[:,0],i1400Array[:,1]])
    iTemp1400 = np.array([i1400Array[:,0],i1400Array[:,2]])
    
    plt.plot(iDens1200[0,:-1], iDens1200[1,:-1], 'ro')
    plt.title("Ion Density vs Radius (1200 G)")
    plt.xlabel("Radius (cm)")
    plt.ylabel("Ion Density")
    rad_plot = "CSDX_radial_plots/ion_dens_v_rad_1200G.jpg"
    plt.savefig(rad_plot)
    plt.show()
    
    plt.plot(iTemp1200[0,:-1], iTemp1200[1,:-1],'ro')
    plt.title("Ion Temperature vs Radius (1200 G)")
    plt.xlabel("Radius (cm)")
    plt.ylabel("Ion Temperature")
    rad_plot = "CSDX_radial_plots/ion_temp_v_rad_1200G.jpg"
    plt.savefig(rad_plot)
    plt.show()
    
    plt.plot(iDens1400[0], iDens1400[1], 'bo')
    plt.title("Ion Density vs Radius (1400 G)")
    plt.xlabel("Radius (cm)")
    plt.ylabel("Ion Density")
    rad_plot = "CSDX_radial_plots/ion_dens_v_rad_1400G.jpg"
    plt.savefig(rad_plot)
    plt.show()
    
    plt.plot(iTemp1400[0], iTemp1400[1],'bo')
    plt.title("Ion Temperature vs Radius (1400 G)")
    plt.xlabel("Radius (cm)")
    plt.ylabel("Ion Temperature")
    rad_plot = "CSDX_radial_plots/ion_temp_v_rad_1400G.jpg"
    plt.savefig(rad_plot)
    plt.show()
    
###############################################################################
   
main()