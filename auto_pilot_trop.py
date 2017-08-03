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
    except:
        print '''
    *******************************************

       Usage: auto_pilot.py [number of years to start] [number of years to end] [acquisition interval(days)] [number of tests] [signals to include]

              number of years to start: minimum length of time series
              number of years to end: maximum length of time series
              acquisition interval: number of days between acquisitions
              number of tests: number of time series to be simulated and run per year
              signals to include: type of the signal to include in the simulated acquisitions (strato, turbulent, deformation, atmosphere_all, all)

    *******************************************
    '''
        sys.exit(1)
    
    for i in xrange(start_yr,end_yr+1):
        n_im = int(math.ceil((365.0/35.0)*i))
        for n in xrange(0,noft):
            x = os.path.basename('TS'+str(i)+'-'+str(n))
#            gname = x+'.template'
#            g = open(gname,'w')
#            g.write('pysar.inputFiles     = /nethome/ehavazli/SIM/synth_data_'+str(i)+str(n)+'/*.unw')
#            g.write('\npysar.corFiles       = /nethome/ehavazli/SIM/synth_data_'+str(i)+str(n)+'/*.cor')
#            g.write('\npysar.seed.yx     = 100 , 100')
#            g.close()

            fname=x+'.process'
            f = open(fname,'w')

            f.write('#! /bin/csh')
            f.write('\n#BSUB -J '+x)
            f.write('\n#BSUB -o '+x+'.o%J')
            f.write('\n#BSUB -e '+x+'.e%J')
            f.write('\n#BSUB -W 24:00')
#            f.write('\n#BSUB -W '+Walltime)
            f.write('\n#BSUB -q general')
            f.write('\n#BSUB -n 1')
  #          f.write('\n#BSUB -R "span[hosts=1]"')
            f.write('\n#BSUB -B')
            f.write('\n#BSUB -N')
#            f.write('\n#BSUB -u e.havazli@rsmas.miami.edu')
            f.write('\ncd /projects/scratch/sinkhole/ehavazli/TSSAR/')
            f.write('\nsimulate_test.py '+str(n_im)+' '+str(n_days)+' ./syn'+str(i)+'-'+str(n)+'/ '+ signal+' 1')
            f.write('\nsim2ref.py ./syn'+str(i)+'-'+str(n)+' 100 100')
            f.write('\ncd ./syn'+str(i)+'-'+str(n))
            f.write('\nload_ref2ts.py .')
            f.write('\ntimeseries2velocity.py timeseries_sim.h5 -o velocity_sim.h5')
#            f.write('\nsimulate.py '+str(n_im)+' '+str(n_days)+' ./syn'+str(i)+str(n)+'/ '+ signal)
#            f.write('\ninterfero.py '+str(i)+str(n))
#            f.write('\nload_data_sim.py ~/SIM/TS'+str(i)+'-'+str(n)+'.template')
#            f.write('\nmv Pairs_'+str(i)+str(n)+'.png $TSSARDIR/TS'+str(i)+'-'+str(n))
#            f.write('\nmv syn'+str(i)+str(n)+' $TSSARDIR/TS'+str(i)+'-'+str(n))
#            f.write('\nmv synth_data'+str(i)+str(n)+' $TSSARDIR/TS'+str(i)+'-'+str(n))
#            f.write('\npysarApp.py ~/SIM/TS'+str(i)+'-'+str(n)+'.template '+ '--dir $TSSARDIR/TS'+str(i)+'-'+str(n))
#            f.write('\nmv Pairs_'+str(i)+str(n)+'.png $TSSARDIR/TS'+str(i)+'-'+str(n))
#            f.write('\nmv syn'+str(i)+str(n)+' $TSSARDIR/TS'+str(i)+'-'+str(n))
#            f.write('\ntimeseries2velocity.py -f $TSSARDIR/TS'+str(i)+'-'+str(n)+'/timeseries.h5 -o $TSSARDIR/TS'+str(i)+'-'+str(n)+'/velocity.h5')

            f.close()
            jobCmd='bsub -P sinkhole < '+fname
            os.system(jobCmd)
#######################################
if __name__ == '__main__':
    main(sys.argv[:])
