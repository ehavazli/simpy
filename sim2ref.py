#! /usr/bin/env python
#@author: Emre Havazli

import os
import sys
from numpy import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import glob
from skimage.restoration import unwrap_phase

def wrap(unw):
    wrapped = unw - around(unw/(2*pi)) * 2*pi
    return wrapped

def main(argv):
    try:
        directory = argv[1]+'/'
        y = int(argv[2])
        x = int(argv[3])
    except:
        print ('''
    *******************************************

       Usage: sim2ref.py [directory] [y] [x]

              directory: directory to read .syn files
              y: seeding pixel y coordinate
              x: seeding pixel x coordinate


    *******************************************
    ''')
        sys.exit(1)
    factor = -1*float(0.0562356467937372)/(4.*pi)
    date=[]
    #directory = './syn30/'
    files_list = glob.glob(directory+'*.syn')
    #y = int(100)
    #x = int(100)

    for day in sorted(files_list):
        date.append(day[-10:-4])
#     date.append(day[-14:-8])
    for i in range(len(date)):
        file1 = open(directory+str(date[0])+'.syn','r')
        file2 = open(directory+str(date[i])+'.syn','r')
        sar1 = fromfile(file1, dtype='float32')
        sar2 = fromfile(file2, dtype='float32')
        sar1 = reshape(sar1,(1000,1000))
        sar2 = reshape(sar2,(1000,1000))
    ##SEEDING##
        sar1_seed = sar1 - sar1[y][x]
        sar2_seed = sar2 - sar2[y][x]
    ##NO SEEDING##
#    sar1_seed = sar1
#    sar2_seed = sar2
    ####
#    print sar1[100][100]
#    print sar2[100][100]
        ts_test = sar1_seed-sar2_seed
#    print ts_test[100][100]
#    ts_test = wrap(ts_test)
#    ts_test_shape = reshape(ts_test,(1000,1000))
#    ts_test_unw = unwrap_phase(ts_test_shape)
#    ts_test_unw = unwrap_phase(ts_test)
#    ts_test_range = ts_test_unw * factor
#    ts_test_range = ts_test * factor

        filename = str(date[i]+'.ref')
        f = open(os.path.join(directory,filename), 'wb')
        ts_test.astype('float32').tofile(f)
#    ts_test_range.astype('float32').tofile(f)
        f.close()

        fig = plt.figure()
        ax = fig.add_axes([0.1,0.1,0.8,0.8])
        title=fig.suptitle(str(date[i]))
#    im = ax.imshow(ts_test_unw,cmap='jet',vmin=-pi,vmax=pi)
        im = ax.imshow(ts_test,cmap='jet',origin='lower')

        cb = fig.colorbar(im)
        fig.savefig(directory+filename+'.png')
        plt.close()
        print ('Referenced syn saved: '+ filename)

#######################################
if __name__ == '__main__':
    main(sys.argv[:])
