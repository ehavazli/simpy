#! /usr/bin/env python
#@author: Emre Havazli

import os
import sys
import h5py
import glob
import math
#import matplotlib.pyplot as plt

##########################
def main(argv):
    try:
        start_yr = int(argv[1])
        end_yr = int(argv[2])
        n_days=int(argv[3])
        noft = int(argv[4])
        signal = argv[5]
        multip = float(argv[6])
    except:
        print ('''
    *******************************************

       Usage: auto_pilot.py [number of years to start] [number of years to end] [acquisition interval(days)] [number of tests] [signals to include] [multiplier]

              number of years to start: minimum length of time series
              number of years to end: maximum length of time series
              acquisition interval: number of days between acquisitions
              number of tests: number of time series to be simulated and run per year
              signals to include: type of the signal to include in the simulated acquisitions (strato, turbulent, deformation, atmosphere_all, all)
              multiplier : 1 for standard turbulent 0.03 for standard stratified. Use stratified multiplier for atmosphere_all option
   *******************************************
    ''')
        sys.exit(1)

    for i in xrange(start_yr,end_yr+1):
        n_im = int(math.ceil((365.0/n_days)*i))
        for n in xrange(0,noft):
            x = os.path.basename('TS'+str(i)+'-'+str(n))
            print ('Number of images: '+str(n_im)+' scenes')
            print ('Acquisition interval: '+str(n_days)+' days')

            fname=x+'.process'
            f = open(fname,'w')
            f.write('#! /bin/csh')
            f.write('\n#BSUB -J '+x)
            f.write('\n#BSUB -o '+x+'.o%J')
            f.write('\n#BSUB -e '+x+'.e%J')
            f.write('\n#BSUB -W 1:00')
            f.write('\n#BSUB -q general')
            f.write('\n#BSUB -n 1')
            f.write('\n#BSUB -B')
            f.write('\n#BSUB -N')
#            f.write('\n#BSUB -u e.havazli@rsmas.miami.edu')
            f.write('\ns.c')
            f.write('\ncd /projects/scratch/sinkhole/ehavazli/TSSAR/')
            f.write('\nsimulate_test.py '+str(n_im)+' '+str(n_days)+' ./syn'+str(i)+'-'+str(n)+'/ '+ signal+' '+str(multip))
            f.write('\nsim2ref.py ./syn'+str(i)+'-'+str(n)+' 100 100')
            f.write('\ncd ./syn'+str(i)+'-'+str(n))
            f.write('\nload_ref2ts.py .')
            f.write('\ntimeseries2velocity.py timeseries_sim.h5 -o velocity_sim.h5')

            f.close()
            jobCmd='bsub -P sinkhole < '+fname
            os.system(jobCmd)
#######################################
if __name__ == '__main__':
    main(sys.argv[:])
