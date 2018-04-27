#! /usr/bin/env python
#@author: Emre Havazli

import os
import sys
import glob
import h5py
from numpy import *
from scipy import stats
from operator import itemgetter
import collections
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

def main(argv):
    try:
        directory = argv[1]
        signal_type = argv[2]
    except:
        print '''
    *******************************************

       Usage: box_plot.py [directory] [signal type]

              directory : location of the TS folders
              signal type: strato, turbulent, combined

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
            for i in spl_lst:
                velocity_file = i +'/velocity_sim.h5'
                f = h5py.File(velocity_file,'r')
                dset = f['velocity'].get('velocity')
                data_to_plot.extend(dset)
            velocities[year] = asarray(data_to_plot)

    sorted_velocities = collections.OrderedDict(sorted(velocities.items()))
    vel_values = []

    for key, values in sorted_velocities.iteritems():
        vel_values.append(values*1000.0) #Convert from meters to milimeters
        labels.append(key)
# Create a figure instance
    fig, ax = plt.subplots(1)
    plt.ylabel('Propagated Error (mm/yr)',fontsize=14)
    plt.xlabel('Time Series Length (years)',fontsize=14)
    plt.ylim(-6,6)
    ax.tick_params(labelsize=12)

    medianprops = dict(linestyle=None, linewidth=0,color = 'red')

    bp = ax.boxplot(vel_values,labels=labels,meanline=True,showmeans=True,showfliers=False,medianprops=medianprops,whis=[2.5, 97.5])
    plt.setp(bp['boxes'], color='blue')
    plt.setp(bp['whiskers'], color='blue', linestyle='--')
    plt.setp(bp['means'], linestyle='-',color = 'red')

    whis = [item.get_ydata()[1] for item in bp['whiskers']]
    n = -2
    unc = []
    for i in range(0,(len(whis)/2)):
        n=n+2
        unc.append(abs((whis[n]) - (whis[n+1])))

    ax2 = ax.twiny()
    ax1Xs = ax.get_xticks()
    ax2.set_xticks(ax1Xs)
    ax2.set_xbound(ax.get_xbound())
    unc_label=[]
    for i in range(0,(len(unc))):
        unc_label.append('('+r'$\pm$'+str(round((unc[i]/2),2))+')')
    ax2.set_xticklabels(unc_label,position=(0.945,0.945))
    ax2.tick_params(direction='in',length=0,labelsize=12)
    if signal_type == 'strato':
        ax2.set_xlabel("Vertical Stratification Only",fontsize=16)
    elif signal_type == 'turbulent':
        ax2.set_xlabel("Turbulence Mixing Only",fontsize=16)
    elif signal_type == 'combined':
        ax2.set_xlabel("Combined Tropospheric Delay",fontsize=16)
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
    fig.savefig('box_plot_VelNorm.tiff', bbox_inches='tight', dpi = 300)
    plt.close()

##Make histograms
#    for i in range(0,len(vel_values)):
#        fig,ax =  plt.subplots(1)
#        bin_values = arange(start=-4, stop=4, step=0.01)
##        vel_values=asarray(vel_values[i])
#        plt.hist(vel_values[i].flatten(),bins=bin_values,normed=1,color='blue',histtype='stepfilled')
#        plt.ylabel('PDF',fontsize=14)
#        plt.xlabel('Velocity (mm/yr)',fontsize=14)
#        fig.savefig('hist_all_'+str(labels[i])+'_years.tiff', bbox_inches='tight', dpi = 300)
#        plt.close()



#######################################
if __name__ == '__main__':
    main(sys.argv[:])
