#! /usr/bin/env python

#@author: Emre Havazli

import os
import sys
from numpy import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import datetime
import glob
import scipy.io
import shutil
from osgeo import gdal

############################################################################################################
def mogi(data, coeffs):
    """
    Computes surface displacements Ux, Uy, Uz in meters from a point spherical
    pressure source in an elastic half space [3].
    evaluate a single Mogi peak over a 2D (2 by N) numpy array of evalpts,
    where coeffs = (x0,y0,z0,dV)
    coeffs = [x0,y0,z0,rs,dp,G]
    """
    dx = data[0,:] - coeffs[0]
    dy = data[1,:] - coeffs[1]
    dz = coeffs[2]
#    c = dV * 3. / 4. * pi
    c = (0.75)*((coeffs[3]**3*coeffs[4])/coeffs[5])
    # or equivalently c= (3/4) a^3 dP / rigidity
    # where a = sphere radius, dP = delta Pressure
    R = sqrt(dx*dx + dy*dy + dz*dz)
    Ux = c*(dx/R**3)
    Uy = c*(dy/R**3)
    Uz = c*(dz/R**3)
    return array((Ux,Uy,Uz))

def wrap(unw):
    wrapped =unw - around(unw/(2*pi)) * 2*pi
    return wrapped

############################################################################################################
def main(argv):
    try:
        n_im = int(argv[1])
        interval = int(argv[2])
        directory = argv[3]+'/'
        inc_signal = argv[4]
    except:
        print ('''
    *******************************************

       Usage: simulate.py [number of images] [image acq interval] [directory] [signals to include] [strato noise level (in meters)/ turbulent multiplier]

       signals to include: strato, turbulent, deformation, atmosphere_all, all

    *******************************************
    ''')
        sys.exit(1)

    syn_file_list = glob.glob(directory+'*.syn')
    orb_file_list = glob.glob(directory+'*.orb')

    if os.path.exists(directory):
        shutil.rmtree(directory, ignore_errors=True)
        os.makedirs(directory)
    else:
        os.makedirs(directory)
############################################################################################################

#####SAR IMAGE#################################################################################################
    x1=linspace(0.,999.,num=1000, endpoint=True) # 1-D space
    y1=linspace(0.,999.,num=1000, endpoint=True)
    [X,Y]=meshgrid(x1,y1) # 2-D space

    range2phase=4*pi/float(0.0562356467937372)         #double-way, 2*2*pi/lamda
    rand_phase = random.random((len(x1),len(x1)))
    rand_phase = rand_phase*range2phase
############################################################################################################

############################################################################################################
# deformation model
#Create grid size
#xx= [linspace(0.,999.,num=1000, endpoint=True)] # 1-D space
#yy = [linspace(0.,999.,num=1000, endpoint=True)] # 1-D space
#[x,y] = meshgrid(xx,yy)

    [x,y] = meshgrid(x1,y1)
    a = len(x)**2
    b = len(y)**2
    coord_x = x.reshape((1,a))
    coord_y = y.reshape((1,b))
    coord = vstack((coord_x,coord_y))
############################################################################################################

############################################################################################################
#create and populate btemp

#    n_im=int(35) ## number of images to be generated
#    interval = int(105) #time interval between images

    btemp= []
    n=0
    A=0
    while n<n_im:
        n+=1
        A+=interval
        btemp.append(A)
############################################################################################################

####IMPORT DEM to GENERATE STRATO###########################################################################

    ds = gdal.Open("socorro_test.tif")
    dem_tif = array(ds.GetRasterBand(1).ReadAsArray())
    dem_inv = (dem_tif)
    norm_max = float(amax(dem_inv))
    norm_min = float(amin(dem_inv))
    dem_inv_norm = dem_inv/(norm_max-norm_min)
    normalizer = norm_min/(norm_max-norm_min)
    strato_norm = dem_inv_norm-normalizer
    strato = 1.0-strato_norm
    #strato = (dem_tif/1000.0)*(-1)

    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    title=fig.suptitle('Inverted and Normalized DEM')
    im = ax.imshow(strato,cmap='terrain',origin='lower')
    cb = fig.colorbar(im)
    #cb.set_label('meters')
    fig.savefig('strato.png')
    plt.close()

    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    title=fig.suptitle('Resampled DEM')
    im = ax.imshow(dem_tif,cmap='gray',origin='lower')
    cb = fig.colorbar(im)
    cb.set_label('meters')
    fig.savefig('dem.png')
    plt.close()

############################################################################################################

    n=0
    for b in btemp:
    #   par = [rx,ry,z,radius,pressure change,modulus(shear or youngs)]
        par = [500,500,100,50,8,32000] ## 100m depth, 50 meters radius,
                                            ## 8(for 2mm/yr, 1200 for 35cm/yr) MPa pressure change, 32GPa shear modulus
                                            ## 190 depth, 500 radius (Socorro values for 100 meter pixel resolution)
                                            ## 0.025 pressure change (1.9 mm/yr pysar estimation, 2mm/yr matlab calculation)
        unit_T = array([-0.38,0.07,-0.92])
        unit = unit_T.reshape((1,3))
        U = mogi(coord,par)
        rangechange = dot(unit,U)
        pred = rangechange.reshape((len(x1),len(x1)))
        deformation = pred.copy()
        deformation = (deformation/float(365))*float(b)
        deformation = range2phase*deformation
        if n <= n_im:
            n += 1

            if inc_signal == 'strato':
    #            random_strat = random.uniform((0),(4*pi))
                strat_lvl = float(argv[5])
                random_strat = random.normal(0,strat_lvl)
#                atmo_noise_strat = random_strat * strato * range2phase
                atmo_noise_strat = random_strat * strato
#                image = wrap((atmo_noise_strat))
                image = atmo_noise_strat
            elif inc_signal == 'turbulent':
                random_turb = int(random.uniform(1,700))
                atmo_turb = scipy.io.loadmat('./atm_all/atmo_surf_turb_'+str(random_turb)+'.mat')
#                atmo_noise_turb = atmo_turb['fsurf'] * 0.1 * range2phase
#                image = wrap(rand_phase+atmo_noise_turb)
                random_multp =random.uniform(0,1)
                turb_lvl = float(argv[5])
                print (random_multp)
                image = ((atmo_turb['fsurf'] * 0.1*random_multp*turb_lvl))
            elif inc_signal == 'deformation':
#                image = wrap(rand_phase+deformation)
                image = (deformation)
            elif inc_signal == 'atmosphere_all':
                random_turb = int(random.uniform(1,700))
                strat_lvl = float(argv[5])
                random_strat = random.normal(0,strat_lvl)
#                atmo_noise_strat = random_strat * strato * range2phase
                atmo_noise_strat = (random_strat * strato)
                atmo_turb = scipy.io.loadmat('./atm_all/atmo_surf_turb_'+str(random_turb)+'.mat')
#                atmo_noise_turb = atmo_turb['fsurf'] * 0.1 * range2phase
                random_multp_turb =random.uniform(0,1)
                atmo_noise_turb = atmo_turb['fsurf'] * 0.1 * random_multp_turb
#                image = wrap(rand_phase+atmo_noise_strat+atmo_noise_turb)
                image = atmo_noise_turb+atmo_noise_strat
            elif inc_signal == 'all':
                strat_lvl = argv[5]
                random_turb = int(random.uniform(1,700))
                random_strat = random.normal(0,strat_lvl)
                atmo_noise_strat = random_strat * strato * range2phase
                atmo_turb = scipy.io.loadmat('./atm_all/atmo_surf_turb_'+str(random_turb)+'.mat')
                atmo_noise_turb = atmo_turb['fsurf'] * 0.1 * range2phase
#                image = wrap(rand_phase+deformation+(atmo_noise_strat)+(atmo_noise_turb))
                image = (deformation+(random_strat * strato)+(atmo_turb['fsurf'] * 0.1))
            else:
                print ('''
    *******************************************

       Usage: simulate.py [number of images] [image acq interval] [directory] [signals to include] [strato noise level (in meters)]

       signals to include: strato, turbulent, deformation, atmosphere_all, all

    *******************************************
    ''')
                sys.exit(1)

############################################################################################################

############################################################################################################
####SAVE DATA####
            start_date = '000101'
            date_1 = datetime.datetime.strptime(start_date, "%y%m%d")
            date_2 = str(date_1 + datetime.timedelta(days=(b)))
            yy = str(date_2[2:4])
            mm = str(date_2[5:7])
            dd = str(date_2[8:10])
            date_2 = yy+mm+dd

            f = open(directory+date_2+'.syn', 'wb')
            image.astype('float32').tofile(f)
            f.close()
#            g.close()

#################
            print (date_2+'.syn saved')
####PLOT####
            if inc_signal == 'strato':
                fig = plt.figure()
                ax = fig.add_axes([0.1,0.1,0.8,0.8])
#                title=fig.suptitle(str(date_2)+" Vertical Stratification Delay Multiplier: "+ str(random_strat))
                # title=fig.suptitle("Vertical Stratification Delay Multiplier: "+ str(round(random_strat*1000.))
                # im = ax.imshow(atmo_noise_strat*1000.,cmap='rainbow',origin='lower')
                title=fig.suptitle("Vertical Stratification Delay Multiplier: "+ str(round(random_strat*1000.,1)) +' mm')
                im = ax.imshow(atmo_noise_strat*1000.,cmap='rainbow',origin='lower')
    #        im = ax.imshow(atmo_noise_strat,cmap='jet',origin='lower')
                cb = fig.colorbar(im)
                cb.set_label('mm')
                fig.savefig(directory+date_2+'.atmo_noise_stra.png')
                plt.close()
#####################
                fig = plt.figure()
                ax = fig.add_axes([0.1,0.1,0.8,0.8])
                title=fig.suptitle(date_2)
#                im = ax.imshow(image,cmap='jet',vmin=-pi,vmax=pi,origin='lower')
                im = ax.imshow(image,cmap='rainbow',origin='lower')
                cb = fig.colorbar(im)
                cb.set_label('Meters')
                fig.savefig(directory+date_2+'.syn.png')
                plt.close()
#####################
#                fig = plt.figure()
#                ax = fig.add_axes([0.1,0.1,0.8,0.8])
#                title=fig.suptitle(date_2)
#                im = ax.imshow(image_unw,cmap='jet',origin='lower')
#                cb = fig.colorbar(im)
#                cb.set_label('[meters]')
#                fig.savefig(directory+date_2+'_unw.syn.png')
#                plt.close()
######################
            elif inc_signal == 'turbulent':
                fig = plt.figure()
                ax = fig.add_axes([0.1,0.1,0.8,0.8])
                # ax.tick_params(labelsize=20)
#                title=fig.suptitle(str(date_2)+" Turbulence Mixing Delay Multiplier: "+str(random_multp))
                title=fig.suptitle("Turbulence Mixing Delay Multiplier: "+str(round(0.1*random_multp,2)))
                im = ax.imshow(image*1000.,cmap='rainbow',origin='lower')
                cb = fig.colorbar(im)
                cb.set_label('mm')
                fig.savefig(directory+date_2+'.atmo_noise_turb.png')
                plt.close()

                fig = plt.figure()
                ax = fig.add_axes([0.1,0.1,0.8,0.8])
                ax.tick_params(labelsize=20)
                title=fig.suptitle(date_2)
#                im = ax.imshow(image,cmap='jet',vmin=-pi,vmax=pi,origin='lower')
                im = ax.imshow(image,cmap='rainbow',origin='lower')
                cb = fig.colorbar(im)
                cb.set_label('mm')
                fig.savefig(directory+date_2+'.syn.png')
                plt.close()
######################

#####################
            elif inc_signal == 'deformation':
                fig = plt.figure()
                ax = fig.add_axes([0.1,0.1,0.8,0.8])
                title=fig.suptitle(str(date_2)+" Deformation")
                im = ax.imshow(deformation,cmap='rainbow',origin='lower')
                cb = fig.colorbar(im)
                cb.set_label('[meters]')
                fig.savefig(directory+date_2+'.deformation.png')
                plt.close()

                fig = plt.figure()
                ax = fig.add_axes([0.1,0.1,0.8,0.8])
                title=fig.suptitle(date_2)
#                im = ax.imshow(image,cmap='jet',vmin=-pi,vmax=pi,origin='lower')
                im = ax.imshow(image,cmap='rainbow',origin='lower')
                cb = fig.colorbar(im)
                cb.set_label('[meters]')
                fig.savefig(directory+date_2+'.syn.png')
                plt.close()
######################
            elif inc_signal == 'atmosphere_all':
                fig = plt.figure()
                ax = fig.add_axes([0.1,0.1,0.8,0.8])
#                title=fig.suptitle(str(date_2)+" Stratified Atmospheric Noise: "+ str(random_strat))
                title=fig.suptitle("Vertical Stratification Delay Multiplier: "+ str(("%.3f" % random_strat)))
                # im = ax.imshow(atmo_noise_turb,cmap='rainbow',origin='lower')
                im = ax.imshow(atmo_noise_strat,cmap='rainbow',origin='lower',vmin = -0.05, vmax = 0.05)
    #        im = ax.imshow(atmo_noise_strat,cmap='jet',origin='lower')
                cb = fig.colorbar(im)
                cb.set_label('Meters')
                fig.savefig(directory+date_2+'.atmo_noise_stra.png')
                plt.close()

                fig = plt.figure()
                ax = fig.add_axes([0.1,0.1,0.8,0.8])
                title=fig.suptitle("Turbulence Mixing Delay Multiplier: "+str(("%.3f" % random_multp_turb)))
                # im = ax.imshow(atmo_noise_turb,cmap='rainbow',origin='lower')
                im = ax.imshow(atmo_noise_turb,cmap='rainbow',origin='lower',vmin = -0.05, vmax = 0.05)
                cb = fig.colorbar(im)
                cb.set_label('Meters')
                fig.savefig(directory+date_2+'.atmo_noise_turb.png')
                plt.close()

                fig = plt.figure()
                ax = fig.add_axes([0.1,0.1,0.8,0.8])
                title=fig.suptitle(date_2)
#                im = ax.imshow(image,cmap='jet',vmin=-pi,vmax=pi,origin='lower')
                # im = ax.imshow(image,cmap='rainbow',origin='lower')
                im = ax.imshow(image,cmap='rainbow',origin='lower',vmin = -0.05, vmax = 0.05)

                cb = fig.colorbar(im)
                cb.set_label('Meters')
                fig.savefig(directory+date_2+'.syn.png')
                plt.close()

            elif inc_signal == 'all':
                fig = plt.figure()
                ax = fig.add_axes([0.1,0.1,0.8,0.8])
                title=fig.suptitle(str(date_2)+" Stratified Atmospheric Noise: "+ str(random_strat))
                im = ax.imshow(atmo_noise_strat,cmap='rainbow',origin='lower')
    #        im = ax.imshow(atmo_noise_strat,cmap='jet',origin='lower')
                cb = fig.colorbar(im)
                cb.set_label('Phase[radians]')
                fig.savefig(directory+date_2+'.atmo_noise_stra.png')
                plt.close()

                fig = plt.figure()
                ax = fig.add_axes([0.1,0.1,0.8,0.8])
                title=fig.suptitle(str(date_2)+" Turbulent Atmospheric Noise")
                im = ax.imshow(atmo_noise_turb,cmap='rainbow',origin='lower')
                cb = fig.colorbar(im)
                cb.set_label('Phase[radians]')
                fig.savefig(directory+date_2+'.atmo_noise_turb.png')
                plt.close()

                fig = plt.figure()
                ax = fig.add_axes([0.1,0.1,0.8,0.8])
                title=fig.suptitle(date_2)
#                im = ax.imshow(image,cmap='jet',vmin=-pi,vmax=pi,origin='lower')
                im = ax.imshow(image,cmap='rainbow',origin='lower')
                cb = fig.colorbar(im)
                cb.set_label('Phase[radians]')
                fig.savefig(directory+date_2+'.syn.png')
                plt.close()

                fig = plt.figure()
                ax = fig.add_axes([0.1,0.1,0.8,0.8])
                title=fig.suptitle(str(date_2)+" Deformation")
                im = ax.imshow(deformation,cmap='rainbow',origin='lower')
                cb = fig.colorbar(im)
                cb.set_label('Phase[radians]')
                fig.savefig(directory+date_2+'.deformation.png')
                plt.close()
######################



        else:
            break

############################################################################################################

############################################################################################################
####Generate Orbit files####
#        files_list = glob.glob(directory+'*.syn')
#        for day in sorted(files_list):
#            orb = int(random.normal(0,300))
#            f = open((directory+day[-10:-4]+'.orb'), 'w')
#        #    f.write('0')
#            f.write(str(orb))
#            f.close()
#            print day[-10:-4]+'.orb saved'
#######################################
if __name__ == '__main__':
    main(sys.argv[:])
