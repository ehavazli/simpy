#! /usr/bin/env python
#@author: Emre Havazli

import os
import sys
import glob
import h5py
from numpy import *
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

       Prepare box plots using the velocity_simStd.h5 files
       Calculation is SD = sqrt(sd^2+...+sd^2)/100
       Usage: box_plot_trop_std.py [directory] [signal type]

              directory : location of the TS folders
              signal type: strato, turbulent, combined

    *******************************************
    '''
        sys.exit(1)

    ts_list = glob.glob(directory+'/syn*')
    labels = []
    standard_dev = {}
    for n in ts_list:
        spl = n.split('syn')
        pick_yr = spl[-1]
        year,test = pick_yr.split('-')
        year = int(year)
        if year in standard_dev: pass
        else:
            print 'Working on '+str(year)+' year long time series'
            spl_lst = glob.glob(directory+'/syn'+str(year)+'*')
            data_to_plot = []
            SD_1 = zeros((1000,1000))
            for i in spl_lst:
                velocity_file = i +'/velocity_simStd.h5'
                f = h5py.File(velocity_file,'r')
                dset = f['velocity'].get('velocity')
                sd_sq = asarray(dset)**2
                data_to_plot.extend(dset)
                SD_1 += sd_sq
            SD = sqrt(SD_1)/len(spl_lst)
            standard_dev[year] = SD

            h5file = 'average_std_'+str(year)+'_years.h5'
            g = h5py.File(h5file,'w')
            hh = g.create_group('velocity')
            dset = hh.create_dataset('velocity', data=SD, compression='gzip')

            for key, value in f['velocity'].attrs.iteritems():
                hh.attrs[key] = value
            g.close()

    sorted_standard_dev = collections.OrderedDict(sorted(standard_dev.items()))
    standard_dev_values = []

    for key, values in sorted_standard_dev.iteritems():
        standard_dev_values.append(values*1000.0) #Convert from meters to milimeters
        labels.append(key)
# Create a figure instance
    fig, ax = plt.subplots(1)
    plt.ylabel('Uncertainty (mm/yr)',fontsize=14)
    plt.xlabel('Time Series Length (years)',fontsize=14)
    plt.ylim(-10,10)
    ax.tick_params(labelsize=12)

    medianprops = dict(linestyle=None,linewidth=0,color = 'red')
    bp = ax.boxplot(standard_dev_values,labels=labels,meanline = True,showmeans=True,showfliers=False,medianprops=medianprops,whis=[2.5, 97.5])
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
        ax2.set_xlabel("Vertical Stratification Uncertainties",fontsize=16)
    elif signal_type == 'turbulent':
        ax2.set_xlabel("Turbulence Mixing Uncertainties",fontsize=16)
    elif signal_type == 'combined':
        ax2.set_xlabel("Combined Tropospheric Uncertainties",fontsize=16)
    else:
        print '''
    *******************************************

       Prepare box plots using the velocity_simStd.h5 files
       Calculation is SD = sqrt(sd^2+...+sd^2)/100
       Usage: box_plot_trop_std.py [directory] [signal type]

              directory : location of the TS folders
              signal type: strato, turbulent, combined

    *******************************************
    '''
        sys.exit(1)


# Save the figure
    fig.savefig('box_plot_AvgStd.tiff', bbox_inches='tight', dpi = 300)
    plt.close()

#######################################
if __name__ == '__main__':
    main(sys.argv[:])
