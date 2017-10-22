#! /usr/bin/env python
#@author: Emre Havazli

import os
import sys
from numpy import *
import matplotlib
import glob
import h5py
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from osgeo import gdal


def main(argv):
    try:
        directory = argv[1]
    except:
        print '''
    *******************************************

       Usage: hist.py [directory]

              directory : location of the TS folders (where syn folders are)

    *******************************************
    '''
        sys.exit(1)


    ts_list = glob.glob(directory+'/syn*')
    velocities = {}
    for n in ts_list:
        spl = n.split('syn')
        pick_yr = spl[-1]
        year,test = pick_yr.split('-')
        year = int(year)
        if year in velocities: pass
        else:
            print 'Working on '+str(year)+' year long time series'
            spl_lst = glob.glob(directory+'/syn'+str(year)+'*')
            data_to_plot = []
            for i in spl_lst:
                velocity_file = i +'/velocity_sim.h5'
                f = h5py.File(velocity_file,'r')
                dset = f['velocity'].get('velocity')
                data_to_plot.extend(dset)
                dset_hist = (asarray(dset)*1000.0)
                bin_values = arange(start=-2, stop=2, step=0.001)
                plt.hist(dset_hist.flatten(),bins=bin_values,normed=1,color='blue',histtype='stepfilled')
                plt.ylabel('PDF',fontsize=14)
                plt.xlabel('Velocity (mm/yr)',fontsize=14)
                plt.savefig(i+'/hist_'+str(i[-4:])+'.tiff', bbox_inches='tight', dpi = 300)
                plt.close()
            velocities[year] = asarray(data_to_plot)
            plt.hist(data_to_plot.flatten,bins=bin_values,normed=1,color='blue',histtype='stepfilled')
            plt.ylabel('PDF',fontsize=14)
            plt.xlabel('Velocity (mm/yr)',fontsize=14)
            plt.savefig('./hist_'+str(year)+'_combined.tiff', bbox_inches='tight', dpi = 300)
            plt.close()
if __name__ == '__main__':
    main(sys.argv[:])
