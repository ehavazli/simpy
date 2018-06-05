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
        print ('''
    *******************************************

       Usage: load_ref2ts.py [directory]

              directory: directory to read .ref files

    *******************************************
    ''')
        sys.exit(1)

#directory = './syn30/'
    dateList = glob.glob(directory+'*.ref')
    h5file = (directory+'timeseries_sim.h5')
    f = h5py.File(h5file,'w')
    d_2 = []
#    gg = f.create_group('timeseries')
    print (dateList)
    dates = []
    data = []
    for date in dateList:
        d = date[-10:-4]
        d_2 = str(datetime.datetime.strptime(d,'%y%m%d').strftime('%Y%m%d'))
        dates.append(d_2)
        if not os.path.basename(date) in f:
            print ('Adding ' + d_2)
            file_open = open(date,'rb')
            file = fromfile(file_open,dtype='float32')
            shape = file.shape[0]
            ax_size = int(sqrt(shape))
            file = asarray(reshape(file,(ax_size,ax_size)))
            data.append(file)
    data = asarray(data)
    print(data.shape)
    dset = f.create_dataset('timeseries', data=data, chunks=True)
    f.attrs['FILE_LENGTH'] = ax_size
    f.attrs['WIDTH'] = ax_size
    f.attrs['XMIN'] = 0
    f.attrs['YMIN'] = 0
    f.attrs['XMAX'] = 999
    f.attrs['YMAX'] = 999
    f.attrs['WAVELENGTH'] = 0.0562356467937372
    f.attrs['DATE12'] = '20'+str((dateList[0][-10:-4])+'-'+d_2)
    f.attrs['date'] = '20'+str(dateList[0][-10:-4])
    f.attrs['REF_DATE'] = '20'+str(dateList[0][-10:-4])
    f.attrs['ref_x'] = '100'
    f.attrs['ref_y'] = '100'
    f.attrs['PROCESSOR'] = 'isce'
    dates = array(dates, dtype=string_)
    dset = f.create_dataset('date', data=dates)


#######################################
if __name__ == '__main__':
    main(sys.argv[:])
