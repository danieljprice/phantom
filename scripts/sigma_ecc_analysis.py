import pysplashsph 
# you need to pip install pysplash
# this is needed to read particle properties
# you can change this to whatever particle 'loader'
# just update compute_hzEnMaeRphi; load_main_quant; load_dump accordingly

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import RegularGridInterpolator,griddata,RectBivariateSpline
from scipy.interpolate import interp1d

def interpolator(x,y,opt):
    option=opt #in case other variables need to be specified from the line
    #first filter with

    #sistema come option che per sigma serve window stretta perche funzione e' molto steep

    y_filtered=savgol_filter(y,window_length=int(len(x)/option),polyorder=2)
    #then interpolate
    interpolation_method = 'cubic'      
    interp = interp1d(x, y_filtered, kind=interpolation_method)
    #interp=np.polynomial.Chebyshev.fit(x,y,opt)
    return interp



def compute_hzEnMaeRphi(dump):
    nsink=int(dump.headers['nptmass'])
    unitv=dump['headers']['udist']/dump['headers']['utime']
    v=np.array([dump['vx'],dump['vy'],dump['vz']])/unitv
    x=np.array([dump['x'],dump['y'],dump['z']])
    h=x[0,:]*v[1,:]-x[1,:]*v[0,:]
    Mtot=np.sum(dump['particlemass'][-nsink:]/dump.headers['umass'])
    #msink required if one does not want to consider
    #mass of acompanion (e.g. external companion)
    rad=np.sqrt(x[0,:]**2+x[1,:]**2)
    En=1./2.*(v[0,:]**2+v[1,:]**2)-Mtot/rad
    a=-Mtot/(2.*En)
    e=np.sqrt(1.-h**2/(Mtot*a))
    phase=np.arctan2(\
                    (-h[:]/Mtot*v[0,:]-x[1,:]/rad),\
                    (h[:]/Mtot*v[1,:]-x[0,:]/rad))
    return h,En,Mtot,a,e,rad,phase,v,x


def binning_rad_faster(Rin,Rout,rad,smoothl,nbin=300):
    partinbin=[]
    binrad=np.linspace(Rin,Rout,nbin)
    dr=(Rout-Rin)/(nbin-1)

    for i in range(0,nbin): #accreted and killed particles (h<0 accreted =0 dea
        partinbin.append(np.nonzero((smoothl[:]>0)*\
                 (rad[:]<binrad[i]+dr)*(rad[:]>binrad[i]))[0])
    #accreted particles
    partinbin.append(np.nonzero(smoothl[:]<0)[0])
    #dead particles
    partinbin.append(np.nonzero(smoothl[:]==0)[0])

    partbinarr=np.array(partinbin[:])
    return partbinarr,binrad,dr


def compute_prof(var,partinbin,cumul=False,complexArr=False):
    nrad=partinbin.shape[0]-2 
    # -2 is because don't want also accreted and dead part
    if complexArr:
        arr=np.zeros(nrad,dtype=np.complex_)
    else:
        arr=np.zeros(nrad)
    if cumul:
        for i in range(0,nrad):
            if(len(partinbin[i]) != 0): #check if array is empty to avoid Error
                arr[i]=np.sum(var[partinbin[i]])
            else:
                arr[i]=0.
    else:
        for i in range(0,nrad):
            if(len(partinbin[i]) !=0): #check if array is empty to avoid Error
                arr[i]=np.mean(var[partinbin[i]])
            else:
                arr[i]=0.
    return arr

def compute_sigma(rad,partinbin,pmass,aprof=False):
    nrad=partinbin.shape[0]-2
    totmass=np.zeros(nrad)
    for i in range(0,nrad):
        totmass[i]=pmass*len(partinbin[i]) 
        #count how many part in bin*pmass
    area=2*np.pi*rad[:]*np.gradient(rad[:])
    if(not aprof): 
        sigma=totmass[:]/area[:]
    else:
    #if calculated on aprof returns fraction of totmass in each bin. It would be better using grad(M(a)) 
        sigma=totmass[:]/sum(totmass)
    return sigma

def load_main_quant(name):
    dump=load_dump(name)
    h,En,Mtot,a,e,rad,phase,v,x=compute_hzEnMaeRphi(dump)
    smoothl=dump['h']
    pmass=dump['headers']['massoftype1']
    npart=int(dump['header']['nparttot'][()])
    Lx=v[2]*x[1]-v[1]*x[2]
    Ly=v[0]*x[2]-v[2]*x[0]
    Lz=v[1]*x[0]-v[0]*x[1]
    l=np.array([Lx,Ly,Lz])
    return x,v,rad,a,e,phase,l,smoothl,dump

def load_dump(name):
    dump=pysplashsph.read.read.read_data_binary(name)
    return dump

if __name__=="__main__":
    try:
        name=sys.argv[1]
        dofileexist=os.path.exists(name)
        if (not dofileexist): 
            print('file passed does not exist')
            sys.exit()
    
    except IndexError as err:
        print('Usage: pythonanalysis.py namedump_000xx')
        raise SystemExit()
    
    
    #this is how you load the particle arrays x,v and smoothing length.
    #density calculations are done from particles within a certain volume.
    Rin=20 
    Rout=170.
    x,v,rad,a,e,phase,l,smoothl,dump=load_main_quant(name)
    pmassArr=dump['particlemass']/dump['headers']['umass']
    #for access to other quantities see dump.headers dictionary variable and dump.labels 
    # e.g. smoothl=dump['h']
 
    #Section where computing all the variables
    
    h,En,Mtot,a,e,r,phase,v,x=compute_hzEnMaeRphi(dump)
    phi=np.arctan2(x[1,:],x[0,:])  

    #binning with radius
    partinbin,radbin,dr=binning_rad_faster(Rin,Rout,r,smoothl)
    #binning with semimajor axis
    partinbinA,radbinA,da=binning_rad_faster(Rin,Rout,a,smoothl)

    #radial binning
    sigma=compute_sigma(radbin,partinbin,pmassArr[0])
    ev=e*np.exp(1j*phase)

    #plots semimajor axis binning
    hprof=compute_prof(h,partinbinA)
    eprof=compute_prof(e,partinbinA)
    evecA=compute_prof(ev,partinbinA,complexArr=True)

    #creating M_a profile
    Maprof=compute_prof(pmassArr,partinbinA,cumul=True)
    Mcumul=[]
    for i in range(len(Maprof)):
        Mcumul.append(np.sum(Maprof[:i]))

    McumulArr=np.array(Mcumul)
    Ma=np.gradient(McumulArr,radbin)
    S0=Ma/(2*np.pi*radbin)

    plt.figure(1)
    plt.plot(radbin,sigma,label='$\Sigma(R)$')     
    plt.plot(radbin,S0,label='$\Sigma_0(a)$')
    plt.xlabel('$R$')
    plt.ylabel('$\Sigma(R)$')
    plt.legend()

    plt.figure(2)
    plt.plot(radbin,eprof,label='$e(a)$') 
    plt.plot(radbin,np.abs(evecA),label='$|e_{vec}|(a)$') 
    plt.legend()
    plt.xlabel('$a$')
    plt.ylabel('$e(a)$')


    plt.show()
