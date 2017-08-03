#! /usr/bin/env python
#@author: Emre Havazli

import os
import sys
import h5py
import glob
import math
import matplotlib.pyplot as plt

##########################
def main(argv):
    try:
        start_yr = int(argv[1])
        end_yr = int(argv[2])
        n_days=int(argv[3])
        noft = int(argv[4])
        std_orb = int(argv[5])
    except:
        print '''
    *******************************************

       Usage: auto_pilot.py [number of years to start] [number of years to end] [acquisition interval(days)] [number of tests] [standard dev for orb]

              number of years to start: minimum length of time series
              number of years to end: maximum length of time series
              acquisition interval: number of days between acquisitions
              number of tests: number of time series to be simulated and run per year
              standard dev for orb: standard deviation for random orbital location generation

    *******************************************
    '''
        sys.exit(1)
    
    for i in xrange(start_yr,end_yr+1):
        n_im = int(math.ceil((365.0/35.0)*i))
        for n in xrange(0,noft):
            x = os.path.basename('TS'+str(i)+'-'+str(n))
            gname = x+'.template'
            g = open(gname,'w')
            g.write('pysar.unwrapFiles     = $TSSARDIR/Test_'+str(i)+'_'+str(n)+'/synth_'+str(i)+'_'+str(n_days)+'/*.unw')
            g.write('\npysar.corFiles       = $TSSARDIR/Test_'+str(i)+'_'+str(n)+'/synth_'+str(i)+'_'+str(n_days)+'/*.cor')
            g.write('\npysar.reference.yx     = 100 , 100')
            g.close()
            mk_test_dir = 'mkdir $TSSARDIR/Test_'+str(i)+'_'+str(n)
            mv_temp ='mv '+str(gname)+' '+'$TSSARDIR/Test_'+str(i)+'_'+str(n)
            os.system(mk_test_dir)
            os.system(mv_temp)

            fname=x+'.process'
            f = open(fname,'w')

            f.write('#! /bin/csh')
            f.write('\n#BSUB -J '+x)
            f.write('\n#BSUB -o '+x+'.o%J')
            f.write('\n#BSUB -e '+x+'.e%J')
            f.write('\n#BSUB -W 24:00')
#            f.write('\n#BSUB -W '+Walltime)
            if i < 7:
                f.write('\n#BSUB -q general')
            else:
                f.write('\n#BSUB -q bigmem')
            f.write('\n#BSUB -n 1')
  #          f.write('\n#BSUB -R "span[hosts=1]"')
            f.write('\n#BSUB -B')
            f.write('\n#BSUB -N')
#            f.write('\ncd ~/SIM')
            f.write('\ncd $TSSARDIR')
#            f.write('\nmkdir Test_'+str(i)+'_'+str(n))
            f.write('\ncd Test_'+str(i)+'_'+str(n))
#            f.write('\nsimulate.py '+str(n_im)+' '+str(n_days)+' ./syn'+str(i)+str(n)+'/ '+ signal+' 0.03')
#            f.write('\nsimulate.py '+str(n_im)+' '+str(n_days)+' ./syn'+str(i)+str(n)+'/ '+ signal)
            f.write('\ninterfero_decor.py '+str(i)+' '+str(n_days)+' '+str(std_orb))
#            f.write('\nload_data_sim.py ./TS'+str(i)+'-'+str(n)+'.template')
#            f.write('\nmv Pairs_'+str(i)+str(n)+'.png $TSSARDIR/TS'+str(i)+'-'+str(n))
#            f.write('\nmv syn'+str(i)+str(n)+' $TSSARDIR/TS'+str(i)+'-'+str(n))
#            f.write('\nmv synth_data'+str(i)+str(n)+' $TSSARDIR/TS'+str(i)+'-'+str(n))
            f.write('\npysarApp_sim.py ./TS'+str(i)+'-'+str(n)+'.template '+ '--dir ./TS'+str(i)+'-'+str(n))
#            f.write('\nmv Pairs_'+str(i)+str(n)+'.png $TSSARDIR/TS'+str(i)+'-'+str(n))
#            f.write('\nmv syn'+str(i)+str(n)+' $TSSARDIR/TS'+str(i)+'-'+str(n))
#            f.write('\ntimeseries2velocity.py -f $TSSARDIR/TS'+str(i)+'-'+str(n)+'/timeseries.h5 -o $TSSARDIR/TS'+str(i)+'-'+str(n)+'/velocity.h5')

            f.close()
            jobCmd='bsub -P sinkhole < '+fname
            os.system(jobCmd)
#######################################
if __name__ == '__main__':
    main(sys.argv[:])
