#! /usr/bin/env python
#@author: Emre Havazli

import os
import sys
from numpy import *
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from osgeo import gdal

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
#                plt.hist(dset_hist.flatten(),normed=1,color='blue',histtype='stepfilled')
#plt.hist(dem_inv.flatten(),bins=bin_values,normed=1,color='blue',histtype='stepfilled')
plt.hist(dem_inv.flatten(),normed=1,color='blue',histtype='stepfilled')
plt.ylabel('PDF',fontsize=14)
plt.xlabel('Height',fontsize=14)
plt.savefig('./hist_dem.tiff', bbox_inches='tight', dpi = 300)
plt.close()

bin_values = arange(start=0, stop=1, step=0.001)
plt.hist(strato.flatten(),normed=1,bins=bin_values,color='blue',histtype='stepfilled')
plt.ylabel('PDF',fontsize=14)
plt.xlabel('Height',fontsize=14)
plt.savefig('./hist_strato.tiff', bbox_inches='tight', dpi = 300)
plt.close()
