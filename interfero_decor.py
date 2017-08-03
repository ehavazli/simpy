#! /usr/bin/env python
#@author: Emre Havazli

import os
import sys
from numpy import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
#import matplotlib.delaunay as md
import matplotlib.tri as trian
import glob
import datetime
from random import randint
from skimage.restoration import unwrap_phase
import shutil
from operator import itemgetter
#################################
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

def all_pairs(date,location):
        result = []
        bperp = []
        source = zip(date,location)
        for p1 in range(len(source)):
                for p2 in range(p1+1,len(source)):
                    result.append([source[p1][0],source[p2][0]])
                    bperp.append(source[p1][1]-source[p2][1])
        return result, bperp


def delaunay_pairs (date,location,perp_max,temp_max,n_im,interval):       
    triang = trian.Triangulation(date,location)
    edges = triang.edges

#    centers,edges,tri,neighbors = md.delaunay(date,location)
    pairs = edges.tolist()
#    print edges

    bperp = []
    for idx in xrange(len(pairs)):
        if pairs[idx][0] > pairs[idx][1]:
            index1=pairs[idx][1]
            index2=pairs[idx][0]
            pairs[idx][0]=index1
            pairs[idx][1]=index2
    pairs_in=sorted(pairs, key=itemgetter(0))
    pairs_in2 = []    
    [pairs_in2.append(i) for i in pairs_in if not i in pairs_in2]    
    pairs = []
    for idx in xrange(len(pairs_in2)):
        perp = location[pairs_in2[idx][1]]-location[pairs_in2[idx][0]]
        abs_perp = abs(perp)
        date_1 = datetime.datetime.strptime(date[pairs_in2[idx][1]], "%y%m%d")
        date_2 = datetime.datetime.strptime(date[pairs_in2[idx][0]], "%y%m%d")
        temp_base = date_2 - date_1
        btemp = abs(temp_base.days)
        if abs_perp > perp_max:
            pass
        elif btemp > temp_max:
            pass
        else:
            bperp.append(perp)
            pairs.append(pairs_in2[idx])
            
    fig = plt.figure()
    ax = fig.add_subplot(111)
    datelist = []
    dates_plot= []
    for i in xrange(len(date)):
        datelist.append(datetime.datetime.strptime(date[i],'%y%m%d'))

    location_base = array(location) - location[0]
#    location_base = array(location)

    ax.plot([datelist],[location_base],'ro')
#    ax.plot([dates_plot],[location_base],'ro')
    for k,j in pairs:
        date_1 = datelist[k]
        date_2 = datelist[j]
#        date_1 = dates_plot[k]
#        date_2 = dates_plot[j]
#        print date_1
        ax.plot([date_1,date_2],[location_base[k],location_base[j]],'b')
#    ax.set_xticklabels(dates_plot)
    fig.suptitle('Baselines')
    plt.ylabel('Bperp')
    plt.xlabel('Time (years)')
#    plt.xticks(rotation = 45)

    plt.xlim([(datelist[0]-(datetime.timedelta(90))),(datelist[-1]+(datetime.timedelta(90)))])
    fig.canvas.draw()
    labels=ax.get_xticklabels()
    new_labels=[]
    for label in labels:
        d_l = str(label.get_text())
        print d_l
        yy = str(d_l[2:4])
        mm = str(d_l[5:7])
#        dd = str(d_l[8:10])
        new_labels.append(yy+mm)
    ax.set_xticklabels(new_labels)
    fig.savefig('Pairs_'+str(n_im)+'_'+str(interval)+'.png')
    
    print 'PAIRS: ' + str(len(pairs))
    return pairs,bperp
    
def wrap(unw):
    wrapped = unw - around(unw/(2*pi)) * 2*pi
    return wrapped

def bperp_decor(baseline,critical_baseline):
    bperp_decor = 1-abs(baseline/critical_baseline)
    if bperp_decor < 0:
        bperp_decor = 0
    else:pass

    return bperp_decor

def btemp_decor(t, tau):
    x = (t/float(tau))*(-1)
    a = (float(0.8)*(exp(x)))-((99999999999999999999**-1)*(t))+0.2
    return a

#def acq_orb(size,orb_std):
#    A = zeros((size,size))
#    O = int(random.normal(0,orb_std))
#    return A,O
#################################
def main(argv):
    try:
        st_yr = int(argv[1])
        interval = int(argv[2])
#        width_length = int(argv[3])
        width_length = 1000
        orb_std = int(argv[3])
        n_im = int(ceil((st_yr * 365.0)/interval)) 
        if os.path.exists('./synth_'+str(st_yr)+'_'+str(interval)+'/'):
                directory = './synth_'+str(st_yr)+'_'+str(interval)+'/'
                print 'Overwriting to directory: '+ str(directory)
        else: 
            directory = './synth_'+str(st_yr)+'_'+str(interval)+'/'
            print 'Creating directory: ' + str(directory)
    except:
        print '''
    *******************************************

       Usage: interfero_decor.py [starting year] [acquisition interval] [standard dev for random orbital location]

              interfero.py 3 35 300 (3 years length, 35 days acq interval, 1000x1000(pixels) image, 300 meters standard dev for random orbital location from a normally distributed set) 

    *******************************************
    '''
        sys.exit(1)
        
    range2phase=4*pi/float(0.0562356467937372)         #double-way, 2*2*pi/lamda
    critical_baseline = float(1100)
#    directory = './syn/'
    date = []
    start_date = datetime.datetime.strptime('000101', "%y%m%d")
    for i in xrange(0,n_im):
        d = str(start_date + datetime.timedelta(days=(i*interval)))
        yy = str(d[2:4])
        mm = str(d[5:7])
        dd = str(d[8:10])
        day = yy+mm+dd
        date.append(day)
              
    location = []
    for i in xrange(0,n_im):
        location.append(int(random.normal(0,orb_std)))
        
#    files_list = glob.glob(directory+'*.syn')
#    for day in sorted(files_list):
#        date.append(day[-10:-4])
#        f=open((directory+day[-10:-4]+'.orb'),'r')
#        orb = int(f.readline())
#        location.append(orb)
        
    #pairs_all,bperp = all_pairs(date,location)
    pairs_del,bperp = delaunay_pairs(date,location,50000,50000,n_im,interval)
    
#    path = 'synth_data_'+str(n_yr)+'/'
    if os.path.exists(directory):
        shutil.rmtree(directory, ignore_errors=True)
        os.makedirs(directory)
    else: 
        os.makedirs(directory)
        
    try:
        pairs_all
    except: 
        NameError
        print 'Delaunay method was selected'
        pairs = pairs_del
    else:
        pairs = []
        for k,j in pairs_all:
            a = date.index(k)
            b = date.index(j)
            pairs.append([a,b])        
    
    print str(len(pairs))+' interferograms are going to be generated'
            
    n=-1
    for k,j in pairs:
#        file1 = open(directory+str(date[k])+'.syn','r')
#        file2 = open(directory+str(date[j])+'.syn','r')    
#        sar1 = fromfile(file1, dtype='float32')
#        width_length = int(sqrt(sar1.shape))
#        sar2 = fromfile(file2, dtype='float32')    
#        interfero = sar2-sar1
        date_1 = datetime.datetime.strptime(date[k], "%y%m%d")
        date_2 = datetime.datetime.strptime(date[j], "%y%m%d")
        delta_t = (date_2-date_1).days
        delta_t = abs(int(delta_t))

        x1=linspace(0.,999.,num=1000, endpoint=True) # 1-D space
        y1=linspace(0.,999.,num=1000, endpoint=True)
        [x,y] = meshgrid(x1,y1)
        a = len(x)**2
        b = len(y)**2
        coord_x = x.reshape((1,a))
        coord_y = y.reshape((1,b))
        coord = vstack((coord_x,coord_y))

        par = [500,500,500,100,215,32000] ## 100m depth, 50 meters radius, 
                                          ## 8(for 2mm/yr, 1200 for 35cm/yr) MPa pressure change, 32GPa shear modulus
                                          ## 190 depth, 500 radius (Socorro values for 100 meter pixel resolution)
                                          ## 0.025 pressure change (1.9 mm/yr pysar estimation, 2mm/yr matlab calculation)
                                          ## 500m depth, 100m radius, 215 MPa delP
        unit_T = array([-0.38,0.07,-0.92])
        unit = unit_T.reshape((1,3))
        U = mogi(coord,par)
        rangechange = dot(unit,U)
        pred = rangechange.reshape((len(x1),len(x1)))
        deformation = pred.copy()
        deformation = (deformation/float(365))*float(delta_t)
        deformation = range2phase*deformation



        interfero = zeros((width_length,width_length))
        n=n+1
#        gamma_bperp = bperp_decor(abs(bperp[n]),critical_baseline)
        gamma_bperp = 0.1
    
#        print 'bperp= '+str(bperp[n])
#        print 'gamma_bperp= '+str(gamma_bperp)
#        date_1 = datetime.datetime.strptime(date[k], "%y%m%d")
#        date_2 = datetime.datetime.strptime(date[j], "%y%m%d")
#        delta_t = (date_2-date_1).days
#        delta_t = abs(int(delta_t))
#        print 'btemp= '+str(delta_t)
#        gamma_btemp = btemp_decor(delta_t,1000)
        gamma_btemp = 0.1
#        print 'gamma_btemp= '+str(gamma_btemp)
        gamma = gamma_bperp*gamma_btemp
#        print 'gamma = '+str(gamma)
        
        wn = random.normal(0,(1.00000000001-gamma),(width_length,width_length))# white noise; 0 mean, ((1-gamma)) std
#        wn = wn * range2phase*0.01 #last multiplier is actually to convert it to from meters to cm
        wn = wn*range2phase*0.01
#        interfero = reshape(interfero,(width_length,width_length))
#        interfero = (deformation+wn)
#        print 'interfero before wrap: '+str(amax(interfero))
#        print 'Maximum WN: '+str(amax(wn))
#        wn_wrap = wrap(wn)
#        print 'WN after wrap: '+str(amax(wn))
#        print 'interfero before wrap: '+str(amax(interfero))
#        interfero = wrap(interfero)
#        print 'interfero after wrap: '+str(amax(interfero))
        wrapped = wrap(deformation+wn)
#        unwrapped = wn
        unwrapped = (deformation+wn)
#        unwrapped = unwrap_phase(wrapped)
#        unwrapped = reshape(unwrapped,(width_length,width_length))
    
        unw_filename = 'filt_'+str(k)+'-'+str(j)+'-sim_HDR_1rlks_c10.unw'      
    
        cor  = ones((width_length,width_length))
        ###Save to file###   
        unw_filename = 'filt_'+str(date[k])+'-'+str(date[j])+'-sim_HDR_1rlks_c10.unw'
        wrp_filename = 'filt_'+str(date[k])+'-'+str(date[j])+'-sim_HDR_1rlks_c10.int'
        cor_filename = 'filt_'+str(date[k])+'-'+str(date[j])+'-sim_HDR_1rlks.cor'    
    
        f = open(os.path.join(directory,unw_filename), 'wb')
        unwrapped.astype('float32').tofile(f)
        f.close()
        
        f = open(os.path.join(directory,wrp_filename), 'wb')
        wrapped.astype('float32').tofile(f)
        f.close()
        
        f = open(os.path.join(directory,cor_filename), 'wb')
        cor.astype('float32').tofile(f)
        f.close()
    
        ###PLOT###
        fig = plt.figure()
        ax = fig.add_axes([0.1,0.1,0.8,0.8])
        title=fig.suptitle(str(date[k])+'_'+str(date[j]))
#        im = ax.imshow(interfero,cmap='jet',vmin=-pi,vmax=pi)
        im = ax.imshow(unwrapped,cmap='jet',origin='lower')
        cb = fig.colorbar(im)
        cb.set_label('Phase[radians]')
        fig.savefig(directory+unw_filename+'.png')
        plt.close()
        print 'Unwrapped interferogram saved: ' + str(date[k]) +'_'+ str(date[j])
        ##################
        fig = plt.figure()
        ax = fig.add_axes([0.1,0.1,0.8,0.8])
        title=fig.suptitle(str(date[k])+'_'+str(date[j]))
        im = plt.imshow(wrapped,cmap='jet',vmin=-pi,vmax=pi,origin='lower')
    #    im = ax.imshow(unwrapped,cmap='jet')
        cb = fig.colorbar(im)
        cb.set_label('Phase[radians]')
        fig.savefig(directory+wrp_filename+'.png')
        plt.close()
        print 'Wrapped interferogram saved: ' + str(date[k]) +'_'+ str(date[j])
        ##################
        fig = plt.figure()
        ax = fig.add_axes([0.1,0.1,0.8,0.8])
        title=fig.suptitle('White Noise gamma = '+str(gamma)+'\n Max= '+ str(amax(wn)))
        im = ax.imshow(wn,cmap='jet',origin = 'lower')
        cb = fig.colorbar(im)
        cb.set_label('Phase [radian]')
        fig.savefig(directory+unw_filename+'_WN.png')
        plt.close()
    
    #####MAKE RSC#####
    
    files_list = glob.glob(os.path.join(directory,'*.unw'))
    count = -1 
    for interferogram in sorted(files_list):
        if count < len(files_list):
            count+=1
            spl = interferogram.split('-')
            day1 = spl[0][-6:]
            day2 = spl[1]
            date1 = datetime.datetime.strptime(day1, "%y%m%d")
            date2 = datetime.datetime.strptime(day2, "%y%m%d")
            time_span_year = ((date2-date1).days/365.00000000)
            date1=str(date1)
            y=date1.split('-')[0]
            m=date1.split('-')[1]
            d=date1.split()[0]
            d=d.split('-')[2]
            
            rsc_file = open(interferogram+'.rsc','w')
            print 'Writing: '+str(interferogram+'.rsc')
            
            rsc_file.write('DELTA_LINE_UTC                           0.024206984            ')     
            rsc_file.write('\nRANGE_PIXEL_SIZE                         62.431792             ')        
            rsc_file.write('\nAZIMUTH_PIXEL_SIZE                       182.839205969682      ')        
            rsc_file.write('\nBASELINE_SRC                             HDR                   ')        
            rsc_file.write('\nRLOOKS                                   1                     ')        
            rsc_file.write('\nALOOKS                                   5                    ')        
            rsc_file.write('\nDOPPLER_RANGE0                           -0.31                 ')        
            rsc_file.write('\nDOPPLER_RANGE1                           -1.69e-05             ')        
            rsc_file.write('\nDOPPLER_RANGE2                           0                     ')        
            rsc_file.write('\nDOPPLER_RANGE3                           0.                    ')        
            rsc_file.write('\nHEADING_DEG                              -13.4690              ')        
            rsc_file.write('\nRGE_REF1                                 822.2015              ')        
            rsc_file.write('\nLOOK_REF1                                15.6107               ')        
            rsc_file.write('\nLAT_REF1                                 33.4789               ')        
            rsc_file.write('\nLON_REF1                                 -107.4901             ')        
            rsc_file.write('\nRGE_REF2                                 870.8046              ')        
            rsc_file.write('\nLOOK_REF2                                23.6596               ')       
            rsc_file.write('\nLAT_REF2                                 33.7160               ')        
            rsc_file.write('\nLON_REF2                                 -106.1368             ')        
            rsc_file.write('\nRGE_REF3                                 822.2015              ')        
            rsc_file.write('\nLOOK_REF3                                15.5926               ')        
            rsc_file.write('\nLAT_REF3                                 33.9455               ')       
            rsc_file.write('\nLON_REF3                                 -107.6137             ')        
            rsc_file.write('\nRGE_REF4                                 870.8046              ')        
            rsc_file.write('\nLOOK_REF4                                23.6482               ')        
            rsc_file.write('\nLAT_REF4                                 34.1829               ')        
            rsc_file.write('\nLON_REF4                                 -106.2521             ')        
            rsc_file.write('\nSTARTING_RANGE1                          822201.4957           ')        
            rsc_file.write('\nSTARTING_RANGE2                          822638.2239           ')        
            rsc_file.write('\nFIRST_LINE_UTC                           17825.7747807006      ')        
            rsc_file.write('\nCENTER_LINE_UTC                          17829.7283862433      ')        
            rsc_file.write('\nLAST_LINE_UTC                            17833.6819917859      ')        
            rsc_file.write('\nSLC_RELATIVE_YMIN                        59311                 ')        
            rsc_file.write('\nAZIMUTH_PIXEL_GROUND                     20.30335              ')        
            rsc_file.write('\nORBIT_NUMBER                             0-0                   ')        
            rsc_file.write('\nI_BIAS                                   15.5                  ')        
            rsc_file.write('\nQ_BIAS                                   15.5                  ')        
            rsc_file.write('\nPLATFORM                                 ENVISAT               ')        
            rsc_file.write('\nSTARTING_RANGE                           822201.4957           ')        
            rsc_file.write('\nPRF                                      1652.415692           ')        
            rsc_file.write('\nAZIMUTH_BANDWIDTH                        1321.93255            ')        
            rsc_file.write('\nAZIMUTH_WEIGHTING                        Kaiser                ')        
            rsc_file.write('\nAZIMUTH_WEIGHTING_PARA                   2.12                  ')        
            rsc_file.write('\nANTENNA_SIDE                             -1                    ')        
            rsc_file.write('\nRANGE_SAMPLING_FREQUENCY                 1.9207680e+07         ')        
            rsc_file.write('\nPLANET_GM                                398600448073000       ')        
            rsc_file.write('\nPLANET_SPINRATE                          7.29211573052e-05     ')        
            rsc_file.write('\nWAVELENGTH                               0.0562356467937372    ')        
            rsc_file.write('\nPULSE_LENGTH                             2.717663e-05          ')        
            rsc_file.write('\nCHIRP_SLOPE                              -588741135306.328     ')        
            rsc_file.write('\nFILE_START                               1                     ')       
            rsc_file.write('\nDATA_TYPE                                CI2                   ')        
            rsc_file.write('\nHEIGHT                                   0.7876952002E+06      ')        
            rsc_file.write('\nHEIGHT_DT                                10.2453660371756      ')        
            rsc_file.write('\nVELOCITY                                 7553.15927579776      ')        
            rsc_file.write('\nLATITUDE                                 32.2110381            ')        
            rsc_file.write('\nLONGITUDE                                -109.5728122          ')        
            rsc_file.write('\nHEADING                                  -13.4061659           ')        
            rsc_file.write('\nEQUATORIAL_RADIUS                        6378137               ')        
            rsc_file.write('\nECCENTRICITY_SQUARED                     0.00669437999014132   ')        
            rsc_file.write('\nEARTH_EAST_RADIUS                        6383781.12727023      ')        
            rsc_file.write('\nEARTH_NORTH_RADIUS                       6353558.09533522      ')        
            rsc_file.write('\nEARTH_RADIUS                             6372091.00925855      ')        
            rsc_file.write('\nORBIT_DIRECTION                          ascending             ')        
            rsc_file.write('\nSQUINT                                   0                     ')        
            rsc_file.write('\nHEIGHT_DS                                0.1356434528E-02      ')        
            rsc_file.write('\nHEIGHT_DDS                               0.2394645375E-09      ')        
            rsc_file.write('\nCROSSTRACK_POS                           -0.1254735009E+03     ')        
            rsc_file.write('\nCROSSTRACK_POS_DS                        0.1700868504E-02      ')        
            rsc_file.write('\nCROSSTRACK_POS_DDS                       -0.5807410711E-08     ')        
            rsc_file.write('\nVELOCITY_S                               7553.1531893          ')        
            rsc_file.write('\nVELOCITY_C                               -0.0000000            ')        
            rsc_file.write('\nVELOCITY_H                               9.5887714             ')        
            rsc_file.write('\nACCELERATION_S                           0.0039038             ')        
            rsc_file.write('\nACCELERATION_C                           -0.5868938            ')        
            rsc_file.write('\nACCELERATION_H                           -7.9639789            ')        
            rsc_file.write('\nVERT_VELOCITY                            0.9122660801E+01      ')        
            rsc_file.write('\nVERT_VELOCITY_DS                         0.3166881638E-05      ')        
            rsc_file.write('\nCROSSTRACK_VELOCITY                      0.1271604283E+02      ')        
            rsc_file.write('\nCROSSTRACK_VELOCITY_DS                   -0.8701671160E-04     ')       
            rsc_file.write('\nALONGTRACK_VELOCITY                      0.7553474724E+04      ')        
            rsc_file.write('\nALONGTRACK_VELOCITY_DS                   -0.2232574231E-05     ')      
            rsc_file.write('\nPEG_UTC                                  17811.782             ')
            
            rsc_file.write('\nXMIN                                     0                    ')        
            rsc_file.write('\nXMAX'+'                                     '+str(int(width_length-1)))        
            rsc_file.write('\nYMIN                                     0                     ')        
            rsc_file.write('\nYMAX'+'                                     '+str(int(width_length-1)))        
            rsc_file.write('\nWIDTH'+'                                    '+str(int(width_length)))        
            rsc_file.write('\nFILE_LENGTH'+'                              '+str(int(width_length)))        
            rsc_file.write('\nTIME_SPAN_YEAR'+'                           '+str(time_span_year))
            rsc_file.write('\nDATE12'+'                                   '+str(day1)+'-'+str(day2))        
            rsc_file.write('\nDATE'+'                                     '+str(day1))       
            rsc_file.write('\nFIRST_LINE_YEAR'+'                          '+str(y))        
            rsc_file.write('\nFIRST_LINE_MONTH_OF_YEAR'+'                 '+str(m))        
            rsc_file.write('\nFIRST_LINE_DAY_OF_MONTH'+'                  '+str(d))        
          
            rsc_file.close()
            
#    files_list = glob.glob(os.path.join(directory,'*.int'))
#    count = -1 
#    for interferogram in sorted(files_list):
#        if count < len(files_list):
#            count+=1
#            spl = interferogram.split('-')
#            day1 = spl[0][-6:]
#            day2 = spl[1]
#            date1 = datetime.datetime.strptime(day1, "%y%m%d")
#            date2 = datetime.datetime.strptime(day2, "%y%m%d")
#            time_span_year = ((date2-date1).days/365.00000000)
#            date1=str(date1)
#            y=date1.split('-')[0]
#            m=date1.split('-')[1]
#            d=date1.split()[0]
#            d=d.split('-')[2]
#            
#            rsc_file = open(interferogram+'.rsc','w')
#            print 'Writing: '+str(interferogram+'.rsc')
#            
#            rsc_file.write('DELTA_LINE_UTC                           0.024206984            ')     
#            rsc_file.write('\nRANGE_PIXEL_SIZE                         62.431792             ')        
#            rsc_file.write('\nAZIMUTH_PIXEL_SIZE                       182.839205969682      ')        
#            rsc_file.write('\nBASELINE_SRC                             HDR                   ')        
#            rsc_file.write('\nRLOOKS                                   1                     ')        
#            rsc_file.write('\nALOOKS                                   5                    ')        
#            rsc_file.write('\nDOPPLER_RANGE0                           -0.31                 ')        
#            rsc_file.write('\nDOPPLER_RANGE1                           -1.69e-05             ')        
#            rsc_file.write('\nDOPPLER_RANGE2                           0                     ')        
#            rsc_file.write('\nDOPPLER_RANGE3                           0.                    ')        
#            rsc_file.write('\nHEADING_DEG                              -13.4690              ')        
#            rsc_file.write('\nRGE_REF1                                 822.2015              ')        
#            rsc_file.write('\nLOOK_REF1                                15.6107               ')        
#            rsc_file.write('\nLAT_REF1                                 33.4789               ')        
#            rsc_file.write('\nLON_REF1                                 -107.4901             ')        
#            rsc_file.write('\nRGE_REF2                                 870.8046              ')        
#            rsc_file.write('\nLOOK_REF2                                23.6596               ')       
#            rsc_file.write('\nLAT_REF2                                 33.7160               ')        
#            rsc_file.write('\nLON_REF2                                 -106.1368             ')        
#            rsc_file.write('\nRGE_REF3                                 822.2015              ')        
#            rsc_file.write('\nLOOK_REF3                                15.5926               ')        
#            rsc_file.write('\nLAT_REF3                                 33.9455               ')       
#            rsc_file.write('\nLON_REF3                                 -107.6137             ')        
#            rsc_file.write('\nRGE_REF4                                 870.8046              ')        
#            rsc_file.write('\nLOOK_REF4                                23.6482               ')        
#            rsc_file.write('\nLAT_REF4                                 34.1829               ')        
#            rsc_file.write('\nLON_REF4                                 -106.2521             ')        
#            rsc_file.write('\nSTARTING_RANGE1                          822201.4957           ')        
#            rsc_file.write('\nSTARTING_RANGE2                          822638.2239           ')        
#            rsc_file.write('\nFIRST_LINE_UTC                           17825.7747807006      ')        
#            rsc_file.write('\nCENTER_LINE_UTC                          17829.7283862433      ')        
#            rsc_file.write('\nLAST_LINE_UTC                            17833.6819917859      ')        
#            rsc_file.write('\nSLC_RELATIVE_YMIN                        59311                 ')        
#            rsc_file.write('\nAZIMUTH_PIXEL_GROUND                     20.30335              ')        
#            rsc_file.write('\nORBIT_NUMBER                             0-0                   ')        
#            rsc_file.write('\nI_BIAS                                   15.5                  ')        
#            rsc_file.write('\nQ_BIAS                                   15.5                  ')        
#            rsc_file.write('\nPLATFORM                                 ENVISAT               ')        
#            rsc_file.write('\nSTARTING_RANGE                           822201.4957           ')        
#            rsc_file.write('\nPRF                                      1652.415692           ')        
#            rsc_file.write('\nAZIMUTH_BANDWIDTH                        1321.93255            ')        
#            rsc_file.write('\nAZIMUTH_WEIGHTING                        Kaiser                ')        
#            rsc_file.write('\nAZIMUTH_WEIGHTING_PARA                   2.12                  ')        
#            rsc_file.write('\nANTENNA_SIDE                             -1                    ')        
#            rsc_file.write('\nRANGE_SAMPLING_FREQUENCY                 1.9207680e+07         ')        
#            rsc_file.write('\nPLANET_GM                                398600448073000       ')        
#            rsc_file.write('\nPLANET_SPINRATE                          7.29211573052e-05     ')        
#            rsc_file.write('\nWAVELENGTH                               0.0562356467937372    ')        
#            rsc_file.write('\nPULSE_LENGTH                             2.717663e-05          ')        
#            rsc_file.write('\nCHIRP_SLOPE                              -588741135306.328     ')        
#            rsc_file.write('\nFILE_START                               1                     ')       
#            rsc_file.write('\nDATA_TYPE                                CI2                   ')        
#            rsc_file.write('\nHEIGHT                                   0.7876952002E+06      ')        
#            rsc_file.write('\nHEIGHT_DT                                10.2453660371756      ')        
#            rsc_file.write('\nVELOCITY                                 7553.15927579776      ')        
#            rsc_file.write('\nLATITUDE                                 32.2110381            ')        
#            rsc_file.write('\nLONGITUDE                                -109.5728122          ')        
#            rsc_file.write('\nHEADING                                  -13.4061659           ')        
#            rsc_file.write('\nEQUATORIAL_RADIUS                        6378137               ')        
#            rsc_file.write('\nECCENTRICITY_SQUARED                     0.00669437999014132   ')        
#            rsc_file.write('\nEARTH_EAST_RADIUS                        6383781.12727023      ')        
#            rsc_file.write('\nEARTH_NORTH_RADIUS                       6353558.09533522      ')        
#            rsc_file.write('\nEARTH_RADIUS                             6372091.00925855      ')        
#            rsc_file.write('\nORBIT_DIRECTION                          ascending             ')        
#            rsc_file.write('\nSQUINT                                   0                     ')        
#            rsc_file.write('\nHEIGHT_DS                                0.1356434528E-02      ')        
#            rsc_file.write('\nHEIGHT_DDS                               0.2394645375E-09      ')        
#            rsc_file.write('\nCROSSTRACK_POS                           -0.1254735009E+03     ')        
#            rsc_file.write('\nCROSSTRACK_POS_DS                        0.1700868504E-02      ')        
#            rsc_file.write('\nCROSSTRACK_POS_DDS                       -0.5807410711E-08     ')        
#            rsc_file.write('\nVELOCITY_S                               7553.1531893          ')        
#            rsc_file.write('\nVELOCITY_C                               -0.0000000            ')        
#            rsc_file.write('\nVELOCITY_H                               9.5887714             ')        
#            rsc_file.write('\nACCELERATION_S                           0.0039038             ')        
#            rsc_file.write('\nACCELERATION_C                           -0.5868938            ')        
#            rsc_file.write('\nACCELERATION_H                           -7.9639789            ')        
#            rsc_file.write('\nVERT_VELOCITY                            0.9122660801E+01      ')        
#            rsc_file.write('\nVERT_VELOCITY_DS                         0.3166881638E-05      ')        
#            rsc_file.write('\nCROSSTRACK_VELOCITY                      0.1271604283E+02      ')        
#            rsc_file.write('\nCROSSTRACK_VELOCITY_DS                   -0.8701671160E-04     ')       
#            rsc_file.write('\nALONGTRACK_VELOCITY                      0.7553474724E+04      ')        
#            rsc_file.write('\nALONGTRACK_VELOCITY_DS                   -0.2232574231E-05     ')      
#            rsc_file.write('\nPEG_UTC                                  17811.782             ')
#            
#            rsc_file.write('\nXMIN                                     0                    ')        
#            rsc_file.write('\nXMAX'+'                                     '+str(int(width_length-1)))        
#            rsc_file.write('\nYMIN                                     0                     ')        
#            rsc_file.write('\nYMAX'+'                                     '+str(int(width_length-1)))        
#            rsc_file.write('\nWIDTH'+'                                    '+str(int(width_length)))        
#            rsc_file.write('\nFILE_LENGTH'+'                              '+str(int(width_length)))        
#            rsc_file.write('\nTIME_SPAN_YEAR'+'                           '+str(time_span_year))
#            rsc_file.write('\nDATE12'+'                                   '+str(day1)+'-'+str(day2))        
#            rsc_file.write('\nDATE'+'                                     '+str(day1))       
#            rsc_file.write('\nFIRST_LINE_YEAR'+'                          '+str(y))        
#            rsc_file.write('\nFIRST_LINE_MONTH_OF_YEAR'+'                 '+str(m))        
#            rsc_file.write('\nFIRST_LINE_DAY_OF_MONTH'+'                  '+str(d))        
#          
#            rsc_file.close()        
            
            rsc_file = open(interferogram[0:-8]+'.cor.rsc','w')
            print 'Writing'+' '+str(interferogram[0:-8]+'.cor.rsc')
            
            rsc_file.write('DELTA_LINE_UTC                           0.024206984            ')     
            
            rsc_file.write('\nRANGE_PIXEL_SIZE                         62.431792             ')        
            rsc_file.write('\nAZIMUTH_PIXEL_SIZE                       182.839205969682      ')        
            rsc_file.write('\nBASELINE_SRC                             HDR                   ')        
            rsc_file.write('\nRLOOKS                                   1                     ')        
            rsc_file.write('\nALOOKS                                   5                     ')        
            rsc_file.write('\nDOPPLER_RANGE0                           -0.31                 ')        
            rsc_file.write('\nDOPPLER_RANGE1                           -1.69e-05             ')        
            rsc_file.write('\nDOPPLER_RANGE2                           0                     ')        
            rsc_file.write('\nDOPPLER_RANGE3                           0.                    ')        
            rsc_file.write('\nHEADING_DEG                              -13.4690              ')        
            rsc_file.write('\nRGE_REF1                                 822.2015              ')        
            rsc_file.write('\nLOOK_REF1                                15.6107               ')        
            rsc_file.write('\nLAT_REF1                                 33.4789               ')        
            rsc_file.write('\nLON_REF1                                 -107.4901             ')        
            rsc_file.write('\nRGE_REF2                                 870.8046              ')        
            rsc_file.write('\nLOOK_REF2                                23.6596               ')       
            rsc_file.write('\nLAT_REF2                                 33.7160               ')        
            rsc_file.write('\nLON_REF2                                 -106.1368             ')        
            rsc_file.write('\nRGE_REF3                                 822.2015              ')        
            rsc_file.write('\nLOOK_REF3                                15.5926               ')        
            rsc_file.write('\nLAT_REF3                                 33.9455               ')       
            rsc_file.write('\nLON_REF3                                 -107.6137             ')        
            rsc_file.write('\nRGE_REF4                                 870.8046              ')        
            rsc_file.write('\nLOOK_REF4                                23.6482               ')        
            rsc_file.write('\nLAT_REF4                                 34.1829               ')        
            rsc_file.write('\nLON_REF4                                 -106.2521             ')        
            rsc_file.write('\nSTARTING_RANGE1                          822201.4957           ')        
            rsc_file.write('\nSTARTING_RANGE2                          822638.2239           ')        
            rsc_file.write('\nFIRST_LINE_UTC                           17825.7747807006      ')        
            rsc_file.write('\nCENTER_LINE_UTC                          17829.7283862433      ')        
            rsc_file.write('\nLAST_LINE_UTC                            17833.6819917859      ')        
            rsc_file.write('\nSLC_RELATIVE_YMIN                        59311                 ')        
            rsc_file.write('\nAZIMUTH_PIXEL_GROUND                     20.30335              ')        
            rsc_file.write('\nORBIT_NUMBER                             0-0                   ')        
            rsc_file.write('\nI_BIAS                                   15.5                  ')        
            rsc_file.write('\nQ_BIAS                                   15.5                  ')        
            rsc_file.write('\nPLATFORM                                 ENVISAT               ')        
            rsc_file.write('\nSTARTING_RANGE                           822201.4957           ')        
            rsc_file.write('\nPRF                                      1652.415692           ')        
            rsc_file.write('\nAZIMUTH_BANDWIDTH                        1321.93255            ')        
            rsc_file.write('\nAZIMUTH_WEIGHTING                        Kaiser                ')        
            rsc_file.write('\nAZIMUTH_WEIGHTING_PARA                   2.12                  ')        
            rsc_file.write('\nANTENNA_SIDE                             -1                    ')        
            rsc_file.write('\nRANGE_SAMPLING_FREQUENCY                 1.9207680e+07         ')        
            rsc_file.write('\nPLANET_GM                                398600448073000       ')        
            rsc_file.write('\nPLANET_SPINRATE                          7.29211573052e-05     ')        
            rsc_file.write('\nWAVELENGTH                               0.0562356467937372    ')        
            rsc_file.write('\nPULSE_LENGTH                             2.717663e-05          ')        
            rsc_file.write('\nCHIRP_SLOPE                              -588741135306.328     ')        
            rsc_file.write('\nFILE_START                               1                     ')       
            rsc_file.write('\nDATA_TYPE                                CI2                   ')        
            rsc_file.write('\nHEIGHT                                   0.7876952002E+06      ')        
            rsc_file.write('\nHEIGHT_DT                                10.2453660371756      ')        
            rsc_file.write('\nVELOCITY                                 7553.15927579776      ')        
            rsc_file.write('\nLATITUDE                                 32.2110381            ')        
            rsc_file.write('\nLONGITUDE                                -109.5728122          ')        
            rsc_file.write('\nHEADING                                  -13.4061659           ')        
            rsc_file.write('\nEQUATORIAL_RADIUS                        6378137               ')        
            rsc_file.write('\nECCENTRICITY_SQUARED                     0.00669437999014132   ')        
            rsc_file.write('\nEARTH_EAST_RADIUS                        6383781.12727023      ')        
            rsc_file.write('\nEARTH_NORTH_RADIUS                       6353558.09533522      ')        
            rsc_file.write('\nEARTH_RADIUS                             6372091.00925855      ')        
            rsc_file.write('\nORBIT_DIRECTION                          ascending             ')        
            rsc_file.write('\nSQUINT                                   0                     ')        
            rsc_file.write('\nHEIGHT_DS                                0.1356434528E-02      ')        
            rsc_file.write('\nHEIGHT_DDS                               0.2394645375E-09      ')        
            rsc_file.write('\nCROSSTRACK_POS                           -0.1254735009E+03     ')        
            rsc_file.write('\nCROSSTRACK_POS_DS                        0.1700868504E-02      ')        
            rsc_file.write('\nCROSSTRACK_POS_DDS                       -0.5807410711E-08     ')        
            rsc_file.write('\nVELOCITY_S                               7553.1531893          ')        
            rsc_file.write('\nVELOCITY_C                               -0.0000000            ')        
            rsc_file.write('\nVELOCITY_H                               9.5887714             ')        
            rsc_file.write('\nACCELERATION_S                           0.0039038             ')        
            rsc_file.write('\nACCELERATION_C                           -0.5868938            ')        
            rsc_file.write('\nACCELERATION_H                           -7.9639789            ')        
            rsc_file.write('\nVERT_VELOCITY                            0.9122660801E+01      ')        
            rsc_file.write('\nVERT_VELOCITY_DS                         0.3166881638E-05      ')        
            rsc_file.write('\nCROSSTRACK_VELOCITY                      0.1271604283E+02      ')        
            rsc_file.write('\nCROSSTRACK_VELOCITY_DS                   -0.8701671160E-04     ')       
            rsc_file.write('\nALONGTRACK_VELOCITY                      0.7553474724E+04      ')        
            rsc_file.write('\nALONGTRACK_VELOCITY_DS                   -0.2232574231E-05     ')      
            rsc_file.write('\nPEG_UTC                                  17811.782             ')
            
            rsc_file.write('\nXMIN                                     0                    ')        
            rsc_file.write('\nXMAX'+'                                     '+str(int(width_length-1)))        
            rsc_file.write('\nYMIN                                     0                     ')        
            rsc_file.write('\nYMAX'+'                                     '+str(int(width_length-1)))        
            rsc_file.write('\nWIDTH'+'                                    '+str(int(width_length)))        
            rsc_file.write('\nFILE_LENGTH'+'                              '+str(int(width_length)))        
            rsc_file.write('\nTIME_SPAN_YEAR'+'                           '+str(time_span_year))
            rsc_file.write('\nDATE12'+'                                   '+str(day1)+'-'+str(day2))        
            rsc_file.write('\nDATE'+'                                     '+str(day1))       
            rsc_file.write('\nFIRST_LINE_YEAR'+'                          '+str(y))        
            rsc_file.write('\nFIRST_LINE_MONTH_OF_YEAR'+'                 '+str(m))        
            rsc_file.write('\nFIRST_LINE_DAY_OF_MONTH'+'                  '+str(d))        
            rsc_file.write('\nTIME_SPAN_YEAR'+'                           '+str(time_span_year))
            rsc_file.write('\nDATE12'+'                                   '+str(day1)+'-'+str(day2))        
            rsc_file.write('\nDATE'+'                                     '+str(day1))       
            rsc_file.write('\nFIRST_LINE_YEAR'+'                          '+str(y))        
            rsc_file.write('\nFIRST_LINE_MONTH_OF_YEAR'+'                 '+str(m))        
            rsc_file.write('\nFIRST_LINE_DAY_OF_MONTH'+'                  '+str(d))        
          
            rsc_file.close()
    
            rsc_file = open(os.path.join(directory,day1+'_'+day2+'_baseline.rsc'),'w')
            print 'Writing'+' '+str(day1+'_'+day2+'_baseline.rsc')
            
            rsc_file.write('TIME_SPAN_YEAR'+'            '+str(time_span_year))
            rsc_file.write('\nH_BASELINE_TOP_HDR            0') 
            rsc_file.write('\nH_BASELINE_RATE_HDR           0')      
            rsc_file.write('\nH_BASELINE_ACC_HDR            0')
            rsc_file.write('\nV_BASELINE_TOP_HDR            0')
            rsc_file.write('\nV_BASELINE_RATE_HDR           0')
            rsc_file.write('\nV_BASELINE_ACC_HDR            0')
            rsc_file.write('\nP_BASELINE_TOP_HDR'+'      '+str(bperp[count]))
            rsc_file.write('\nP_BASELINE_BOTTOM_HDR'+'   '+str(bperp[count]))
            rsc_file.write('\nORB_SLC_AZ_OFFSET_HDR         0')                          
            rsc_file.write('\nORB_SLC_AZ_OFFSET_TOP_HDR     0')
            rsc_file.write('\nORB_SLC_AZ_OFFSET_BOTTOM_HDR  0')
            rsc_file.write('\nORB_SLC_R_OFFSET_HDR          0')                           
            rsc_file.write('\nORB_SLC_R_OFFSET_NEAR_HDR     0')
            rsc_file.write('\nORB_SLC_R_OFFSET_FAR_HDR      0')
            rsc_file.write('\nRANGE_OFFSET_HDR              0')
            rsc_file.write('\nPHASE_CONST_HDR               -99999')
                           
            rsc_file.close()
        else: break
#######################################
if __name__ == '__main__':
    main(sys.argv[:])        

