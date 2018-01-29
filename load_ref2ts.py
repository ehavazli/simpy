#! /usr/bin/env python
#@author: Emre Havazli

import sys
import os
from numpy import *
import h5py
import glob
import datetime

def main(argv):
    try:
        directory = argv[1]+'/'
    except:
        print '''
    *******************************************

       Usage: load_ref2ts.py [directory]
      
              directory: directory to read .ref files

    *******************************************
    '''
        sys.exit(1)

#directory = './syn30/' 
    dateList = glob.glob(directory+'*.ref')
    h5file = (directory+'timeseries_sim.h5')
    f = h5py.File(h5file,'w')
    d_2 = []
    gg = f.create_group('timeseries')
    print dateList
    for date in dateList:
        d = date[-10:-4]
        d_2 = str(datetime.datetime.strptime(d,'%y%m%d').strftime('%Y%m%d'))
        if not os.path.basename(date) in f:
            print 'Adding ' + d_2
#        group = gg.create_group(os.path.basename(d_2))
            file_open = open(date,'rb')
            file = fromfile(file_open,dtype='float32')
            shape = file.shape[0]
            ax_size = int(sqrt(shape))
            file = reshape(file,(ax_size,ax_size))
            dset = gg.create_dataset(d_2, data=file, compression='gzip')
            gg.attrs['FILE_LENGTH'] = ax_size
            gg.attrs['WIDTH'] = ax_size
            gg.attrs['XMIN'] = 0
            gg.attrs['YMIN'] = 0
            gg.attrs['XMAX'] = 999
            gg.attrs['YMAX'] = 999
            gg.attrs['WAVELENGTH'] = 0.0562356467937372
            gg.attrs['DATE12'] = '20'+str((dateList[0][-10:-4])+'-'+d_2)
            gg.attrs['DATE'] = '20'+str(dateList[0][-10:-4])
            gg.attrs['ref_date'] = '20'+str(dateList[0][-10:-4])
            gg.attrs['ref_x'] = '100'
            gg.attrs['ref_y'] = '100'
#######################################
if __name__ == '__main__':
    main(sys.argv[:])
