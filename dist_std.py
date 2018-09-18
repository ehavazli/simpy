#! /usr/bin/env python2
#@author: Emre Havazli

import os
import sys
from numpy import *
from scipy.optimize import leastsq
import glob
import h5py
import random
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

def main(argv):
    try:
        directory = argv[0]
        e_y = int(argv[1])
        e_x = int(argv[2])
        modelFit = argv[3]
    except:
        modelFit = ''
        print '''
    *******************************************

       Usage: dist_std.py [directory] [end_y] [end_x] [modelFit]
       e_y = end iteration at this y
       e_x = end iteration at this x
       modelFit = log (logarothmic), pow (power-law)
    *******************************************
    '''
        sys.exit(1)

    files = sorted(glob.glob(directory+'/average_std_*_years.h5'))
    x_res = float(183.25)
    y_res = float(333)
    std_pix = {}
    for k in files:
        file_name = k.split('_')[-2:]
        years = file_name[0]
        print 'Working on '+str(years)+' years long time series'
        f = h5py.File(k,'r')
        dset = f['velocity'].get('velocity')
        stda = asarray(dset)
        std_lst = []
        dist = []
        for i in xrange(0,10):
            s_y = random.randint(0,999)
            s_x = random.randint(0,999)
            ref_std = stda[s_y][s_x]
            # print str(s_y)+' '+str(s_x)
            for y in xrange(0,e_y):
                for x in xrange(0,e_x):
                    std_val = abs(stda[y][x]-ref_std)
                    std_lst.append(std_val)
                    dist.append(sqrt(((s_y-y)*y_res)**2+((s_x-x)*x_res)**2))
        # print years+' '+str(len(std_lst))+' '+str(len(dist))
        std_pix[int(years)] = [std_lst,dist]
    x_first = arange(0,50000,5000)
    x_sec = arange(40000,100000,10000)
    x_thrd = arange(100000,360000,20000)
    x = concatenate((x_first,x_sec[1:],x_thrd),axis=0)
    n_range = len(x) - 1
    color=iter(cm.jet_r(linspace(0,1,8)))

    for key, value in sorted(std_pix.iteritems()):
        print 'Working on plots of '+str(key)+' years data'
        c=next(color)
        avg = []
        avg_std=[]
        avg_er = []
        for n in xrange(0,n_range):
            print 'Distances between '+str(x[n])+' and '+str(x[n+1])
            plt.ylabel('Uncertainty (mm/yr)')
            plt.xlabel('Distance (km)')
            for i, q in enumerate(value[1]):
                if x[n] <= value[1][i] <= x[n+1]:
                    avg_std.append(value[0][i])
            avg.append(nanmean(avg_std))
            avg_er.append(std(avg_std))
#        avg.insert(0,0)
#        avg_er.insert(0,0)
        x_plt = x[1:]
        # print 'X: '+str(x_plt)
        # print 'AVG: '+str(len(avg))
        # print 'PLOTTED: '+str(array(avg)*100000.0)
        # print 'AVG_ER: '+str(avg_er)
        # plt.scatter(x/1000.0,array(avg)*10000.0,c=c,label=(str(key)+' years'))
        if modelFit == 'log':
            testY = (array(avg)*100000.0)
            testX = x_plt/1000.0
            testYerr = array(avg_er)*100000.0
            logFit = polyfit(log(testX),testY,1)
            testFit = logFit[0]*log(testX)+logFit[1]
            plt.errorbar(testX,testY,testYerr,c=c,marker='o',xerr = None,ls='none',label=(str(key)+' years'))
            plt.plot(testX,testFit,c=c,label=(str(round(logFit[1],2))+'+'+str(round(logFit[0],2))+'log(D)'))
        elif modelFit == 'pow':
            # Define function for calculating a power law
            powerlaw = lambda x, amp, index: amp * (x**index)
            logx = log10(x_plt/1000.0)
#            logx = logx[isnan(logx)]
            logy = log10(array(avg)*100000.0)
#            logy = logy[isnan(logy)]
            yerr = array(avg_er)*100000.0
            xdata = x_plt/1000.0
            ydata = array(avg)*100000.0
            logyerr = yerr / ydata

            # define our (line) fitting function
            fitfunc = lambda p, x: p[0] + p[1] * x
            errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

            pinit = [1.0, -1.0]
            out = leastsq(errfunc, pinit,
                       args=(logx, logy, logyerr), full_output=1)

            pfinal = out[0]
            covar = out[1]
            print pfinal
            print covar

            index = pfinal[1]
            amp = 10.0**pfinal[0]

            indexErr = sqrt( covar[1][1] )
            ampErr = sqrt( covar[0][0] ) * amp
            ##########
            # Plotting data
            ##########
            plt.plot(xdata, powerlaw(xdata, amp, index),c=c)     # Fit
            plt.errorbar(xdata, ydata, yerr=yerr, c=c,marker='o',xerr = None,ls='none',label=(str(key)+' years'))  # Data
            # plt.text(5, 6.5, 'Ampli = %5.2f +/- %5.2f' % (amp, ampErr))
            # plt.text(5, 5.5, 'Index = %5.2f +/- %5.2f' % (index, indexErr))
        else:
            plt.errorbar(x_plt/1000.0,array(avg)*100000.0,array(avg_er)*100000.0,c=c,marker='o',xerr = None,ls='none',label=(str(key)+' years'))
    plt.legend(ncol=2,loc=2)
    plt.axis([0, 350,0,10])
    plt.savefig(directory+'/STDvsDIST.png',bbox_inches="tight",dpi=600)





#     color=iter(cm.jet(linspace(0,1,8)))
#     for key, value in sorted(std_pix.iteritems()):
#         c=next(color)
#         plt.ylabel('Uncertainty (mm/yr)')
#         plt.xlabel('Distance (km)')
# #        plt.axis([0, 350,0,2])
#         plots = plt.plot(array(value[1])/1000.0,array(value[0])*10000.0, '*',c=c,label=(str(key)+' years'))
# #        plt.plot(sol,c=c,label=(str(key)+' years'))
#         plt.legend()
#     plt.savefig(directory+'/STDvsDIST.png',bbox_inches="tight",dpi=600)
# #    plt.close()

            # out_file = open(directory+'std_dst_'+years+'_'+str(y)+'_'+str(s_y)+'-''+str(s_x)+'.txt', 'w')
            # std_pix_ar = array(std_pix)*float(10000.0)
            # for item in std_pix_ar:
            #     out_file.write("%s\n" % item)
            # out_file.close




###########################
if __name__ == '__main__':
    main(sys.argv[1:])
