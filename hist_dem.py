#! /usr/bin/env python
#@author: Emre Havazli

import os
import sys
from numpy import *
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from osgeo import gdal
import h5py

ds = gdal.Open("socorro_test.tif")
dem_tif = array(ds.GetRasterBand(1).ReadAsArray())
dem_inv = (dem_tif)
norm_max = float(amax(dem_inv))
norm_min = float(amin(dem_inv))
dem_inv_norm = dem_inv/(norm_max-norm_min)
normalizer = norm_min/(norm_max-norm_min)
strato_norm = dem_inv_norm-normalizer
strato = 1.0-strato_norm

bin_values = arange(start=1000, stop=3500, step=1)
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.hist(dem_inv.flatten(),bins=bin_values,normed=1,color='grey',histtype='stepfilled')
plt.ylabel('PDF',fontsize=14)
plt.xlabel('Height',fontsize=14)

velocity_file ='./velocity_sim.h5'
f = h5py.File(velocity_file,'r')
dset = f['velocity'].get('velocity')
dset_hist = (asarray(dset)*1000.0)
bin_values = arange(start=-3, stop=3, step=0.001)
ax2 = fig.add_subplot(212)
ax2.hist(dset_hist.flatten(),bins=bin_values,normed=1,color='blue',histtype='step')
plt.ylabel('PDF',fontsize=14)
plt.xlabel('Velocity',fontsize=14)
fig.tight_layout()
plt.savefig('./hist_dem.tiff', bbox_inches='tight', dpi = 300)
plt.close()
