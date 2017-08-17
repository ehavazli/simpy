#! /usr/bin/env python2
#@author: Emre Havazli

import sys
from numpy import *
import h5py
import pysar
import pysar._readfile as readfile
import pysar._pysar_utilities as ut
import csv

def main(argv):
    try:
        velocityfile = argv[0]
        y = int(argv[1])
        x = int(argv[2])
        snum_pix = int(argv[3])
        enum_pix = int(argv[4])
    except:
        print '''
    *******************************************

       Usage: average_circle.py h5file y x start_num_pix end_num_pix

    *******************************************
    '''
        sys.exit(1)

    f = h5py.File(velocityfile,'r')
    dset = f['velocity'].get('velocity')
    velocity = asarray(dset)
    ref_vel = (velocity[y][x])*1000.0 #convert velocity values from meters to mm
#    print 'REF VEL: ' + str(ref_vel)
    atr = readfile.read_attribute(velocityfile)
    ts_year = velocityfile.split('_')[2]
    csv_file = open('summary_'+str(ts_year)+'years_mean_std.csv',"wb")
    writer = csv.writer(csv_file)
    writer.writerow(('Distance (km)', 'Time Series Length (years)','Average (mm/yr)', 'Standard Deviation (mm/yr)'))

    for i in range(snum_pix,(enum_pix+1)):
        dist= (i/10.0) #convert distance to km
        crc_mask = ut.circle_index(atr,(y,x,i))
        masked_vel = (velocity[crc_mask])*1000.0 #convert velocity values from meters to mm
    #    print 'MASKED VEL: ' +str(masked_vel)
        diff_masked_vel = (ref_vel)-(masked_vel)
    #    print 'DIFF MASKED VEL: '+str(diff_masked_vel)
        average_crc = mean(diff_masked_vel,dtype=float64)
        std_masked = std(diff_masked_vel,dtype=float64)
        writer.writerow((dist,ts_year,average_crc,std_masked))
    #    print 'Masked Vel: ' + str(masked_vel)
    #    print 'Referenced Vel: ' +str(diff_masked_vel)
    #    print 'Average of referenced vel in the circle: ' +str(average_crc)+' mm/yr'
    #    print 'Standard Dev of referenced vel in the Circle: ' +str(std_masked)
    csv_file.close()
###########################
if __name__ == '__main__':
    main(sys.argv[1:])
