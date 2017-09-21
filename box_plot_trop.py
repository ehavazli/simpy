#! /usr/bin/env python
#@author: Emre Havazli

import os
import sys
import glob
import h5py
from numpy import *
from operator import itemgetter
import collections
import matplotlib.pyplot as plt

def main(argv):
    try:
        directory = argv[1]
        signal_type = argv[2]
    except:
        print '''
    *******************************************

       Usage: box_plot.py [directory]

              directory : location of the TS folders

    *******************************************
    '''
        sys.exit(1)
        

    ts_list = glob.glob(directory+'/syn*')
    labels = []
    velocities = {}
    for n in ts_list:
        spl = n.split('syn')
        pick_yr = spl[-1]
        year,test = pick_yr.split('-')
        year = int(year)
        if year in velocities: pass
        else:
            print 'Working on '+str(year)+' year long time series'
            spl_lst = glob.glob(directory+'/syn'+str(year)+'*')
            data_to_plot = []
            average = zeros(1000,1000) 
            for i in spl_lst:
                velocity_file = i +'/velocity_sim.h5'  
                f = h5py.File(velocity_file,'r')
                dset = f['velocity'].get('velocity')
                data_to_plot.extend(dset)
                average += dset
            velocities[year] = asarray(data_to_plot) 

# Create a figure instance
#    fig, ax = plt.subplots(1)
#    plt.ylabel('mm/yr')
#    plt.xlabel('Time Series Length')
    sorted_velocities = collections.OrderedDict(sorted(velocities.items()))
    vel_values = []

    for key, values in sorted_velocities.iteritems():
        vel_values.append(values*1000.0) #Convert from meters to milimeters
        labels.append(key)

# Create a figure instance
    fig, ax = plt.subplots(1)
    plt.ylabel('mm/yr')
    plt.xlabel('Time Series Length')
    plt.ylim(-10,10)
    ax.tick_params(labelsize=20)
#    plt.xlim([(amin(vel_values)-2),(amax(vel_values)+2)])
    bp = ax.boxplot(vel_values,labels=labels, showfliers=False,whis=[2.5, 97.5])
    plt.setp(bp['boxes'], color='blue')
    plt.setp(bp['whiskers'], color='blue', linestyle='--')
    plt.setp(bp['medians'], color = 'red')
    whis = [item.get_ydata()[1] for item in bp['whiskers']]
    n = -2 
    unc = []
    for i in range(0,(len(whis)/2)):
        n=n+2
        unc.append((abs(whis[n])) + (abs(whis[n+1])))

    ax2 = ax.twiny()
    ax1Xs = ax.get_xticks()
    ax2.set_xticks(ax1Xs)
    ax2.set_xbound(ax.get_xbound())
    unc_label=[]
    for i in range(0,(len(unc))):
        unc_label.append('('+r'$\pm$'+str(round((unc[i]/2),2))+')')
    ax2.set_xticklabels(unc_label,position=(0.945,0.945))
    ax2.tick_params(direction='in',length=0)
    if signal_type == 'strato':
        ax2.set_xlabel("Vertical Stratification Uncertainties")
    elif signal_type = 'turbulent':
        ax2.set_xlabel("Turbulence Mixing Uncertainties")
    elif signal_type = 'combined':
        ax2.set_xlabel("Combined Tropospheric Uncertainties")
    else:
                print '''
    *******************************************

       Usage: box_plot_trop.py [directory] [signal type]

              directory : location of the TS folders
              signal type: strato, turbulent, combined
    *******************************************
    '''
        sys.exit(1)


# Save the figure
    fig.savefig('box_plot.tiff', bbox_inches='tight', dpi = 300)
    plt.close()
#######################################
if __name__ == '__main__':
    main(sys.argv[:])  
