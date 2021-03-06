#! /usr/bin/env python2
#@author: Emre Havazli

import os
import sys
import h5py
import glob
from numpy import *
from pysar._writefile import write

directory='./'
ts_list = glob.glob(directory+'/syn*')
velocities = {}
atr={}
atr['FILE_LENGTH'] = 1000
atr['WIDTH'] = 1000
#atr['XMIN'] = 0
#atr['YMIN'] = 0
#atr['XMAX'] = 999
#atr['YMAX'] = 999
#atr['WAVELENGTH'] = 0.0562356467937372
atr['FILE_TYPE'] = 'velocity'

for n in ts_list:
    spl = n.split('syn')
    pick_yr = spl[-1]
    year,test = pick_yr.split('-')
    year = int(year)
    summ = zeros((1000,1000))
    if year in velocities: pass
    else:
        print 'Working on '+str(year)+' year long time series'
        velocities[year] = zeros((1000,1000))
        spl_lst = glob.glob(directory+'/syn'+str(year)+'*')
        for i in spl_lst:
            velocity_file = i +'/velocity_sim.h5'
            f = h5py.File(velocity_file,'r')
            dset = f['velocity'].get('velocity')
            vel = asarray(dset)
            velocities[year]=velocities[year]+vel
        average = velocities[year]/float(len(spl_lst))
#        average = velocities[year]
        filename = 'average_vel_'+str(year)+'_years.h5'
        try:
            os.remove(filename)
        except OSError:
            pass
            write(average,atr,filename)
