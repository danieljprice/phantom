import matplotlib.pyplot as plt
import numpy as np
import random
import math
from scipy.special import i0,i1,iv,k0,k1,kn
#import asciitable
from scipy.integrate import quad as spiq
from scipy.integrate import trapz as trap
from scipy.integrate import simps as simp
from scipy.stats import maxwell as mbd
import sys

def readinoptions():
  print "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
  print " __  __ _  ______  _    _ _      _____ "
  print "|  \/  | |/ /  _ \| |  | | |    / ____|"
  print "| \  / | ' /| |_) | |  | | |   | |  __ "
  print "| |\/| |  < |  _ <| |  | | |   | | |_ |"
  print "| |  | | . \| |_) | |__| | |___| |__| |"
  print "|_|  |_|_|\_\____/ \____/|______\_____|"
  print "                                       "
  print "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"

                                       
                                       
  inputfile = (sys.argv)[1]
  ans1      = (sys.argv)[2]
  ansS      = (sys.argv)[3]
  answers=open(inputfile)
  ans=[]
  for line in answers:
    ans.append((line.split())[0])
  Bulge_ID  = ans[0]
  ansP      = ans[1]
  Npart_PB  = int(ans[2])
  Npart_VB  = int(ans[3])
  Npart_PD  = int(ans[4])
  Npart_VD  = int(ans[5])
  startogas = float(ans[6])
  softlen   = float(ans[7])   #In pc
  bulgeR    = float(ans[8])   #In kpc
  ansJ      = ans[9]
  stream    = float(ans[10])
  rBtrunc   = float(ans[11])
  rscatter  = float(ans[12])
  DiscRatio = (ans[13])[0]
  GassDisc  = float(ans[14])
  
  print '-=-=-=-=-=-=-=-=-=-=-=-=-=-'
  if ansS=='d':
    print 'SETTING THE DISC STARS/GAS...'
  elif ansS=='b':
    print 'SETTING THE BULGE STARS...'
    print '**WARNING**'
    print 'Make sure the disc has been set first.'
  else:
    print "Fine. Be like that. I'm not setting anything then."
  
  if Bulge_ID=='h':
    print "Henquist bulge selected,  Phi = MG / (r+a)"
  elif Bulge_ID=='p':
    print "Plummer bulge selected,  Phi = MG / sqrt(r^2+a^2)"
  elif Bulge_ID=='t':
    print "Truncated Plummer bulge selected,  Phi = MG / sqrt(r^2+a^2)"
  else:
    print "**ERROR**"
    print "If you don't want a bulge use the other code dumbass..."
    exit()
  
  
  if ans1=='o':
    print 'Plotting RC curve alone'
  if ans1=='p':
    print 'Drawing positions'
  if ans1=='v':
    print 'Drawing velocities (must be done after positions)'

  if ansP=='p':
    print 'Plotting as we go'
  else:
    print 'Not plotting (likely too many paritcles for matplotlib!)'    
  
  print 'Setting:',Npart_PB,'bulge particle positions'
  print 'Setting:',Npart_VB,'bulge particle velocities (Npart_P<=Npart_V)'
  print 'Setting:',Npart_PD,'disc particle positions'
  print 'Setting:',Npart_VD,'disc particle velocities (Npart_P<=Npart_V)'

  if (Npart_VD>Npart_VD) or (Npart_VB>Npart_VB):
    print '**ERROR**'
    print 'Npart_V must be <= Npart_P or positions will not be available!'
    exit()

  print 'With a star to gas ratio of:',startogas  
  if startogas>0.5:
    print '**WARNING**'
    print 'Are you sure you want more stars than gas? Useful mainly for testing.'
  
  print 'Smoftening length for mass binning of:',softlen,'pc'
  print 'Extent of the bulge:',bulgeR,'kpc'
  print 'Truncation of the bulge (double exp decay):',rBtrunc,'kpc'
  if ansJ=='p':
    print 'Integtating Jeans eq. through potential derivatives (analytic)'
  elif ansJ=='m':
    print 'Integtating Jeans eq. through binning mass components (numeric)'
  else:
    print 'Pick a relevant method for solving the Jeans eq. please...'
    exit()

  if stream>0.5:
    print 'The fraction of bulge that is rotating is',stream
  elif stream<0.5:
    print 'Pretty sure you dont want the bulge to counter-rotate...'
    print 're-try with stream>0.5 (stream=',stream,'currently)'
    exit()
  elif stream==0.5:
    print 'The bulge is NOT rotating (stream=0.5)'

  print 'Distance for scattering in distance arrays:',rscatter,' kpc'
  
  if DiscRatio=='l':
    print 'Using the lowest Md/Mh value'
  elif DiscRatio=='n':
    print 'Using the normal Md/Mh value (from Wada paper)'
  elif DiscRatio=='h':
    print 'Using the heavy Md/Mh value (my fiducial value)'
  elif DiscRatio=='v':
    print 'Using the very heavy Md/Mh value'
  else:
    print "Ivalid Md/Mh value found, assuming fiducial value ('heavy')"
    DiscRatio = 'h'
  
  print 'Using a gas disc mass of:',GassDisc,'x10^9 Mo'
  
  print '-=-=-=-=-=-=-=-=-=-=-=-=-=-'

  return ansS,Bulge_ID,ans1,ansP,Npart_PB,Npart_VB,Npart_PD,Npart_VD,startogas,softlen,bulgeR,ansJ,stream,rBtrunc,rscatter,DiscRatio,GassDisc
 
"<><><><><><><><><><><><><><><><><><><><><><><><><><><>"
"<><><><><><><><><><><><><><><><><><><><><><><><><><><>"
"In SI units:"
Mo  = 1.99e30
pc  = 3.0856e16
kpc = 1.0e3 * pc
G   = 6.67e-11
km  = 1.0e3
"<><><><><><><><><><><><><><><><><><><><><><><><><><><>"
"<><><><><><><><><><><><><><><><><><><><><><><><><><><>"
ansS,Bulge_ID,ans1,ansP,Npart_PB,Npart_VB,Npart_PD,Npart_VD,StarToGas,softlen,bulgeR,ansJ,stream,rBtrunc,rscatter,DiscRatio,GassDisc = readinoptions()

"Variables (set to Baba et al. 2009 rather than Hernquist 1993 values)"
r1  = 0.0 * kpc
r2  = 13. * kpc         #The full disk radius we will use
r2b = bulgeR * kpc      #The radius of the placement of bulge particles
rBtrunc = rBtrunc * kpc 

rscatter = rscatter*kpc  #radius to shuffle particles in Gaussian to take into account finite r array

dsoft = softlen*pc
disc_contractor = 0.35 #Percentage to contract disc by...

z1  = 0.  
zfr = 0.50
z2  = +zfr*r2
Mg  = GassDisc * 1.0e9 * Mo
"halo:"
r200= 122 * kpc   #Also called the virial radius.
Cnfw= 5.
rh  = r200/Cnfw

if Bulge_ID=='h':
  "disk:"
  r0  = 3.0 * kpc
  z0  = 0.3 * kpc
  Md  = 3.2e10 * Mo
  "bulge:"
  a0  = 0.30 * kpc
  Mb  = 1.5e10 * Mo


  #From Kafle2012
  r0  = 3.5 * kpc
  z0  = 0.3 * kpc
  Md  = 6.5e10 * Mo
  "bulge:"
  a0  = 0.50 * kpc
  Mb  = 1.8e10 * Mo


if (Bulge_ID=='p') or (Bulge_ID=='t'):
  "disk:"
  r0  = 3.0 * kpc
  z0  = 0.3 * kpc
  "bulge:"
  a0  = 0.35 * kpc
  Mb  = 1.05e10 * Mo

  if DiscRatio=='l':
    print "Assigning 'light' disc mass values"
    Md = 2.5e10     * Mo  #x0.7 Wada value
    Mh =10.1e11     * Mo  #x1.3*1.3 Wada
  if DiscRatio=='n':
    print "Assigning 'normal' disc mass values"
    Md = 3.2e10     * Mo  #Wada value
    Mh = 8.3e11     * Mo  #x1.3 heavy
  if DiscRatio=='h':
    print "Assigning 'heavy' disc mass values"
    Md = 4.1e10     * Mo  #x1.3 Wada value
    Mh = 6.3e11     * Mo  #Wada value
  if DiscRatio=='v':
    print "Assigning 'Vheavy' disc mass values"
    Md = 5.3e10     * Mo  #x1.3 heavy
    Mh = 4.4e11     * Mo  #x0.7 heavy
  '''
  "Wada"
  #Md  = 2.5e10 * Mo
  #Md  = 3.2e10 * Mo
  "Used to use"
  #Md  = 4.1e10 * Mo#, Can have x1.3 to slightly better match the RC (halo uses atm)
  #Md  = 5.3e10 * Mo#, Can have x1.3 x1.3 to slightly better match the RC (halo uses atm)
  #Clares settings:
  #a0 = 0.5*kpc 
  #Mb = 1.4e10*Mo
  "New Try..., a middle ground between exact and flat"
  #Md  = 3.8e10 * Mo  
  #a0 = 0.40*kpc 
  #Mb = 1.2e10*Mo
  '''


rho0= (Mh/(4*math.pi*r200**3))  *  Cnfw**3/(np.log(1.+Cnfw)+Cnfw/(1.+Cnfw))
vescFac = 0.95

"Orbital Params:"
Qref= 1.0
Rref= 1.5 * r0

"<><><><><><><><><><><><><><><><><><><><><><><><><><><>"
"<><><><><><><><><><><><><><><><><><><><><><><><><><><>"
def Sofue12():
  f=open('../Table_GrandRC.dat.txt')
  r =np.zeros(129)
  v =np.zeros(129)
  ev=np.zeros(129)
  i=0
  for line in f:
    r[i] =line.split()[0]
    v[i] =line.split()[1]
    ev[i]=line.split()[2]
    i+=1
  plt.errorbar(r,v,yerr=ev,xerr=None,fmt='bo',label='Sofue12 (GRC)',markersize=2)

def sech(x):
  s = 1./np.cosh(x)
  return s

def RhoDisk(r,z):
  C = Md/(4.*math.pi*r0*r0*z0)
  rho = C*np.exp(-r/r0)*np.power(sech(z/z0),2.)
  return rho
  
def RhoHalo(r):
  rho = rho0 / ( (r/rh) * (1.+r/rh)**2 ) # *np.exp(-(r/rHtrunc)**2)
  return rho

def SigDisk(r):
  C = Md/(2.*math.pi*r0*r0)
  sigma = C*np.exp(-r/r0)
  return sigma

def SigDiskz(z):
  C = Md/(4.*math.pi*r0*r0)   #Not exact constant but it will be integrated out in pdf.
  sigma = C*np.power(sech(z/z0),2.)
  return sigma

def pdfDisk(r):
  NormConst = Normaliser(1)
  N_r   = r*np.exp(-r/r0)
  pdf   = NormConst*N_r
  return pdf

def pdfDiskz(z):
  NormConst = Normaliser(2)
  N_z   = np.power(sech(z/z0),2.)
  pdf   = NormConst*N_z
  return pdf
  
def cdfDisk(r):
  NormConst = Normaliser(1)
  N_rdr   = -r0*np.exp(-r/r0)*(r0+r)
  r   = 0.
  N_rdr0  = -r0*np.exp(-r/r0)*(r0+r)
  cdf = NormConst*(N_rdr-N_rdr0)
  return cdf

def cdfDiskz(z):
  "NOTE: the z-dist doesnt need a *z as r does, because it is already pre one units dist (Sigma was 1/r^2)"
  NormConst = Normaliser(2)
  N_zdz   = z0*np.tanh(z/z0)  #z0*z*np.tanh(z/z0) - z0*z0*np.log(np.cosh(z/z0))
  z   = 0.
  N_zdz0  = z0*np.tanh(z/z0)  #z0*z*np.tanh(z/z0) - z0*z0*np.log(np.cosh(z/z0))
  cdf = NormConst*(N_zdz-N_zdz0)
  return cdf
  
def RhoBulge(r):
  if Bulge_ID=='h':
    rho = Mb*a0/(2.*math.pi*r*(r+a0)**3)
  if Bulge_ID=='p':
    rho = 0.75*Mb*a0*a0 /(math.pi* ( r*r + a0*a0 )**(5./2.) )
  if Bulge_ID=='t':
    rratio = r/rBtrunc
    rho = 0.75*Mb*a0*a0 /(math.pi* ( r*r + a0*a0 )**(5./2.) ) * np.exp(-(r/rBtrunc)**2)
    #rho = 0.75*Mb*a0*a0 / ( r*r + a0*a0 )**(5./2.)  * (1.-np.power(rratio,2)+np.power(rratio,4)/2. - np.power(rratio,6)/6. + np.power(rratio,8)/24. - np.power(rratio,10)/120. + np.power(rratio,12)/720. - np.power(rratio,14)/5040. + np.power(rratio,16)/40320. - np.power(rratio,18)/362880.+ np.power(rratio,20)/3628800. - np.power(rratio,22)/39916800. + np.power(rratio,24)/479001600. - np.power(rratio,26)/6227020800. ) 
  return rho


def pdfBulge(r):
  "Norm const. contains the G, M and 4pi term form dS."
  if Bulge_ID=='h':
    NormConst = Normaliser(3)
    N_r   = (1./(r*(r+a0)**3))  *  r**2
    pdf   = NormConst*N_r
  if Bulge_ID=='p':
    NormConst = Normaliser(4)
    N_r   = (1./ ( r*r + a0*a0 )**(5./2.)) *  r**2
    pdf   = NormConst*N_r
  if Bulge_ID=='t':
    NormConst = Normaliser(5)
    rratio = r/rBtrunc
    N_r   = (1./ ( r*r + a0*a0 )**(5./2.)) * np.exp(-(r/rBtrunc)**2) *  r**2 
    #N_r   = (1./ ( r*r + a0*a0 )**(5./2.)) * (1.-np.power(rratio,2)+np.power(rratio,4)/2. - np.power(rratio,6)/6. + np.power(rratio,8)/24. - np.power(rratio,10)/120. + np.power(rratio,12)/720. - np.power(rratio,14)/5040. + np.power(rratio,16)/40320. - np.power(rratio,18)/362880.+ np.power(rratio,20)/3628800. - np.power(rratio,22)/39916800. + np.power(rratio,24)/479001600. - np.power(rratio,26)/6227020800. ) *  r**2 
    #The expansion of the truncation:  1-x^2/a^2+x^4/(2 a^4)-x^6/(6 a^6)+x^8/(24 a^8)-x^10/(120 a^10)+x^12/(720 a^12)-x^14/(5040 a^14)+x^16/(40320 a^16)-x^18/(362880 a^18)+x^20/(3628800 a^20)-x^22/(39916800 a^22)+x^24/(479001600 a^24)-x^26/(6227020800 a^26)+x^28/(87178291200 a^28)-x^30/(1307674368000 a^30)+x^32/(20922789888000 a^32)+O(x^33)
    #Then +O(x^18)
    pdf   = NormConst*N_r/(kpc*0.3)
  return pdf

def cdfBulge(r):
  if Bulge_ID=='h':
    NormConst = Normaliser(3)
    x      = r/a0
    N_rdr  = (-1.0/a0)*( x/(2.*(1+x)**2) + 1/(2.*(1+x)) )
    r      = 0.
    x      = r/a0
    N_rdr0 = (-1.0/a0)*( x/(2.*(1+x)**2) + 1/(2.*(1+x)) )
    cdf = NormConst*(N_rdr-N_rdr0)
  if Bulge_ID=='p':
    NormConst = Normaliser(4)
    x      = r/a0
    N_rdr  = (1/a0)*(x*x*x / ( 3.*(1.+x*x)**(3./2.) ))
    r      = 0.
    x      = r/a0
    N_rdr0 = (1/a0)*(x*x*x / ( 3.*(1.+x*x)**(3./2.) ))
    cdf = NormConst*(N_rdr-N_rdr0)
  if Bulge_ID=='t':
    NormConst = Normaliser(5)
    #Need to int: 1/sqrt(a)  *  (x^2/(x^2+1)^(5./2.)) * ( 1-x^2/a^2+x^4/(2 a^4)-x^6/(6 a^6)+x^8/(24 a^8)-x^10/(120 a^10)+x^12/(720 a^12)-x^14/(5040 a^14)+x^16/(40320 a^16)-x^18/(362880 a^18)+x^20/(3628800 a^20)-x^22/(39916800 a^22)+x^24/(479001600 a^24)-x^26/(6227020800 a^26) )
    x      = r/a0
    #rratio = x
    N_rdr  = spiq( lambda r: ( (1./ ( r*r + a0*a0 )**(5./2.)) * np.exp(-(r/rBtrunc)**2) *  r**2) , 0., r)[0] 
    #r      = 0.
    #x      = r/a0
    #N_rdr0 = 0.
    cdf = NormConst*(N_rdr)
    #print cdf
  return cdf

def Normaliser(flag):
  if flag==1:
    "Radial Disk"
    r=r2
    c2 = -r0*np.exp(-r/r0)*(r0+r)
    r=r1
    c1 = -r0*np.exp(-r/r0)*(r0+r)
    return 1./(c2-c1)
  if flag==2:
    "Vertical Disk"
    z=z2
    c2 = z0*np.tanh(z/z0)  #z0*z*np.tanh(z/z0) - z0*z0*np.log(np.cosh(z/z0))
    z=z1
    c1 = z0*np.tanh(z/z0)  #z0*z*np.tanh(z/z0) - z0*z0*np.log(np.cosh(z/z0))
    return 1./(c2-c1)
  if flag==3:
    "Hernquist Bulge normalisation"
    r=r2b
    x=r/a0
    c2 =  (-1.0/a0)*( x/(2.*(1+x)**2) + 1/(2.*(1+x)) )
    r=r1
    x=r/a0
    c1 =  (-1.0/a0)*( x/(2.*(1+x)**2) + 1/(2.*(1+x)) )
    return 1./(c2-c1)
  if flag==4:
    "Plummer Bulge normalisation"
    r=r2b
    x=r/a0
    c2 =  (1/a0)*(x*x*x / ( 3.*(1.+x*x)**(3./2.) ))
    r=r1
    x=r/a0
    c1 =  (1/a0)*(x*x*x / ( 3.*(1.+x*x)**(3./2.) ))
    return 1./(c2-c1)
  if flag==5:
    "Truncated Plummer Bulge normalisation"
    r=r2b
    x1=r/a0
    #c2 =  (1/a0)*(x*x*x / ( 3.*(1.+x*x)**(3./2.) ))
    r=r1
    x2=r/a0
    #c1 =  (1/a0)*(x*x*x / ( 3.*(1.+x*x)**(3./2.) ))
    ansSPIQ = 1./spiq( lambda r: ( (1./ ( r*r + a0*a0 )**(5./2.)) * np.exp(-(r/rBtrunc)**2) *  r**2) , r1, r2b)[0] 
    return ansSPIQ

def randomdist(r,cd,N):
  Rpred=np.zeros(N)
  "Pull N random particles to place"
  j=0
  while j<N:
	z=np.random.rand()
	k=0
	while  (k<len(r)) and (z>=cd[k]):
	  k+=1
	if k==len(r):
	  Rpred[j] = r[-1]
	else:
	  Rpred[j] = r[k]
	j+=1
  return Rpred

def gausfunc(x,mu,sigma):
  y=np.exp(-(x/sigma)**2)/np.sqrt(2.*math.pi*sigma**2)
  return y

def readindisk(diskfile,Ndisc):
  #Assign masses randomly (shouldn't matter):
  diskarray = np.zeros((Ndisc,2))
  i=0
  #while i<np.shape(diskarray)[0]:
  #  x_di= discput[i][0]
  #  y_di= discput[i][1]
  #  z_di= discput[i][2]
  #  m_di= discput[i][3]
  discput = open('asciifile_D')
  print discput
  for line in discput:
    x_di= float((line.split())[0])
    y_di= float((line.split())[1])
    z_di= float((line.split())[2])
    m_di= float((line.split())[3])
    #Convert cylindrical polars to spherical polar co-ordinates:
    r_sph= np.sqrt(x_di*x_di+y_di*y_di+z_di*z_di)
    diskarray[i,0] = r_sph * (1. - disc_contractor)
    diskarray[i,1] = m_di
    i+=1
  discput2 = open('asciifile_G')
  print discput2
  for line in discput2:
    x_di= float((line.split())[0])
    y_di= float((line.split())[1])
    z_di= float((line.split())[2])
    m_di= float((line.split())[3])
    #Convert cylindrical polars to spherical polar co-ordinates:
    r_sph= np.sqrt(x_di*x_di+y_di*y_di+z_di*z_di)
    diskarray[i,0] = r_sph * (1. - disc_contractor)
    diskarray[i,1] = m_di
    i+=1
  return diskarray

def COUNTtheDISK(r):
  #Want to count the mass within a certain radius, r.
  "THIS ARRAY IS RADIUS BY MASS"
  InnerArray = diskarray[( diskarray[:,0] <= (r+dsoft) )]
  countedMass= InnerArray[:,1].sum()
  #print countedMass/(1.e10*Mo),r/kpc
  return countedMass

def MassDisk(r):
  #MassinD = COUNTtheDISK(r)  #Would need to count EACH radius call, slow.
  
  #iii  =0
  #while r_int[iii]<r:  #Move up until break radius of that mass bin.
  #  iii+=1
  #MassinD = MassinDisc_Ar[iii-1]
  
  MassinD = MassinDisc_Ar[int(r/drint)] 
  
  return MassinD

def MassBulge(r):
  if Bulge_ID=='h':
    MassinB=MassHERNBulge(r)
  if Bulge_ID=='p':
    MassinB=MassPLUMBulge(r)
  if Bulge_ID=='t':
    MassinB=MassPLUMBulge(r)
  return MassinB

def MassHERNBulge(r):
  "Hernquist Bulge"
  MassinB = Mb*(r**2)/(r+a0)**2
  return MassinB

def MassPLUMBulge(r):
  "Plummer Bulge"
  MassinB = Mb*(r**3)/(r*r+a0*a0)**(3./2.)
  return MassinB

def MassHalo(r):
  MassinH = 4.*math.pi*rho0*rh**3 *(np.log((rh+r)/rh) - r/(rh+r))
  "Above taken from: Zavala et al. 2006."
  return MassinH

def PlotMasses(r):
  #Test massinobjects
  diskarray = readindisk(diskfile,Ndisc)   #Returns an array of [r,m] in spherical polars
  m1=np.zeros(len(r))
  ij=0
  for line in m1:
    m1[ij] = MassDisk(r[ij])
    ij+=1
  m2 = MassBulge(r)
  m3 = MassHalo(r)
  #plt.figure(20)
  #plt.plot(r/kpc,m1/(1.e10*Mo),'ko',label='Disc')
  #plt.plot(r/kpc,m2/(1.e10*Mo),'go',label='Bulge')
  #plt.plot(r/kpc,m3/(1.e10*Mo),'ro',label='Halo')
  #plt.legend()
  #plt.savefig('MassCDF.png')
  
def PotDisk(n,r):
  "n: order of differentiation"
  rratio  =0.5*r/r0
  if n==1:
    dPhidr_disk  =(G*Md*r/(2.*r0*r0*r0)) * (i0(rratio)*k0(rratio)-i1(rratio)*k1(rratio))
    return dPhidr_disk
  elif n==2:
    d2Phidr2_disk=(G*Md  /(4.*r0*r0*r0)) * (-rratio*iv(2,rratio)*k1(rratio) + i0(rratio)*(2*k0(rratio)-3.*rratio*k1(rratio)) + i1(rratio)*(3.*rratio*k0(rratio)-2.*k1(rratio)+rratio*kn(2,rratio)) )
    return d2Phidr2_disk
  elif n=='v':
    dPhidr_disk  =(G*Md*r/(2.*r0*r0*r0)) * (i0(rratio)*k0(rratio)-i1(rratio)*k1(rratio))
    Vc  =  np.sqrt(r*dPhidr_disk)
    return Vc
  else:
    Phi_disk     =(-G*Md*r/(2.*r0*r0)) * (i0(rratio)*k1(rratio)-i1(rratio)*k0(rratio))
    return Phi_disk
  
def PotHalo(n,r):
  "n: order of differentiation"
  if n==1:
    dPhidr_halo  =( -np.log(1.+r/rh)/(r*r) + 1./((rh*r)*(1.+r/rh)) ) * -G*Mh/(np.log(Cnfw+1.)-Cnfw/(Cnfw+1.))
    return dPhidr_halo
  elif n==2:
    d2Phidr2_halo=( +2*np.log(1.+r/rh)/(r*r*r) - 1./((rh*rh*r)*(1.+r/rh)*(1.+r/rh)) - 2./((rh*r*r)*(1.+r/rh)) ) * -G*Mh/(np.log(Cnfw+1.)-Cnfw/(Cnfw+1.))
    return d2Phidr2_halo
  elif n=='v':
    dPhidr_halo  =( -np.log(1.+r/rh)/(r*r) + 1./((rh*r)*(1.+r/rh)) ) * -G*Mh/(np.log(Cnfw+1.)-Cnfw/(Cnfw+1.))
    Vc  =  np.sqrt(r*dPhidr_halo)
    return Vc
  else:
    Phi_halo     =(  np.log(1.+r/rh)/r ) * -G*Mh/(np.log(Cnfw+1.)-Cnfw/(Cnfw+1.))
    return Phi_halo

def PotBulge(n,r):
  if Bulge_ID=='h':
    potout=PotHERNBulge(n,r)
  if Bulge_ID=='p':
    potout=PotPLUMBulge(n,r)
  if Bulge_ID=='t':
    potout=PotPLUMBulge(n,r)
  return potout

def PotHERNBulge(n,r):
  "n: order of differentiation"
  "Hernquist Bulge"
  print Mb/Mo,a0/kpc
  if n==1:
    dPhidr_bulge  =+G*Mb/(r+a0)**2
    return dPhidr_bulge
  elif n==2:
    d2Phidr2_bulge=-2.*G*Mb/(r+a0)**3
    return d2Phidr2_bulge
  elif n=='v':
    dPhidr_bulge  =+G*Mb/(r+a0)**2
    Vc  =  np.sqrt(r*dPhidr_bulge)
    return Vc
  else:
    Phi_bulge     =-G*Mb/(r+a0)
    return Phi_bulge

def PotPLUMBulge(n,r):
  "n: order of differentiation"
  "Plummer Bulge"
  r2term = np.sqrt(r*r + a0*a0)
  if n==1:
    dPhidr_bulge  =+G*Mb*r/(r2term*r2term*r2term)
    return dPhidr_bulge
  elif n==2:
    d2Phidr2_bulge=-G*Mb*(2.*r*r-a0*a0)/(r2term*r2term*r2term*r2term*r2term)
    return d2Phidr2_bulge
  elif n=='v':
    dPhidr_bulge  =+G*Mb*r/(r2term*r2term*r2term)
    Vc  =  np.sqrt(r*dPhidr_bulge)
    return Vc
  else:
    Phi_bulge     =-G*Mb/r2term
    return Phi_bulge

def kappa_freq(r):
  "Epicycle frequency"
  dPhidr_disk   =PotDisk(1,r)
  d2Phidr2_disk =PotDisk(2,r)
  dPhidr_halo   =PotHalo(1,r)
  d2Phidr2_halo =PotHalo(2,r)
  dPhidr_bulge  =PotBulge(1,r)
  d2Phidr2_bulge=PotBulge(2,r)
  kappa = np.sqrt( (dPhidr_disk+dPhidr_halo+dPhidr_bulge)*3./r  +  (d2Phidr2_disk+d2Phidr2_halo+d2Phidr2_bulge)  )
  return kappa

def kappa_freqDisc(r):
  "Epicycle frequency"
  dPhidr_disk   =PotDisk(1,r)
  d2Phidr2_disk =PotDisk(2,r)
  dPhidr_halo   =PotHalo(1,r)
  d2Phidr2_halo =PotHalo(2,r)
  dPhidr_bulge  =0.#PotBulge(1,r)
  d2Phidr2_bulge=0.#PotBulge(2,r)
  kappa = np.sqrt( (dPhidr_disk+dPhidr_halo+dPhidr_bulge)*3./r  +  (d2Phidr2_disk+d2Phidr2_halo+d2Phidr2_bulge)  )
  return kappa

def omega_freq(r):
  "Circular frequency"
  Vcirc = Vc(r)
  omega = Vcirc/r
  return omega

def Vc(r):
  "Circular velocity"
  Vcirc = np.sqrt(r*(PotDisk(1,r) + PotHalo(1,r) + PotBulge(1,r)))
  return Vcirc
    
def ToomreConstFn(r):
  kappa_Rref = kappa_freq(Rref)
  Vref  = Qref * 3.36 * G * SigDisk(Rref)/kappa_Rref
  Const = (Vref*Vref) * np.exp(Rref/r0)
  return Const

def BulgeJeansEq_P(x):
  fnout = RhoBulge(x) * (PotDisk(1,x) + PotHalo(1,x) + PotBulge(1,x))
  return fnout

def BulgeJeansEq_M(x):
  fnout = RhoBulge(x) * G * ( MassHalo(x) + MassBulge(x) + MassDisk(x) )/(x*x)
  return fnout

print 'The bar stability parameter;'
print 'eta_b = Vmax / SQRT( G*Md/Rd ), is' 
Vmax = np.max( Vc(np.arange(0.1,13.0,0.1)*kpc) )
print 'Vmax=',Vmax/1000.,'km/s'
etab = Vmax/np.sqrt( G*Md/r0 )
print etab
if (etab<1.1):
  print 'Unstable to bar formation (etabar<1.1)'
else:
  print 'Stable to bar formation (etabar>1.1)'

"<><><><><><><><><><><><><><><><><><><><><><><><><><><>"
"RUN TIME"
"<><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  
if ans1=='o':
  print '------------------------------------'
  print 'Plotting Rc curve, no initialisation'
  print '------------------------------------'
  r=np.arange(0.001,20.0,0.001)*kpc
  p1 = PotDisk('v',r)
  p2 = PotHalo('v',r)
  p3 = PotBulge('v',r)
  pt = np.sqrt(np.power(p1,2)+np.power(p2,2)+np.power(p3,2))
  plt.figure(1)
  plt.plot(r/kpc,p1/km,'b--')
  plt.plot(r/kpc,p2/km,'r--')
  plt.plot(r/kpc,p3/km,'g--')
  plt.plot(r/kpc,pt/km,'k-')
  #Sofue12()
  plt.xlabel('R [kpc]')
  plt.ylabel('Vc [km/s]')


  e1 = PotDisk('p',r)
  e2 = PotHalo('p',r)
  e3 = PotBulge('p',r)
  et = np.sqrt( -2.*(e1+e2+e3) )
  plt.plot(r/kpc,vescFac*np.sqrt(-2.*e1)/km,'b-')
  plt.plot(r/kpc,vescFac*np.sqrt(-2.*e2)/km,'r-')
  plt.plot(r/kpc,vescFac*np.sqrt(-2.*e3)/km,'g-')
  plt.plot(r/kpc,vescFac*et/(km),'c-')

  plt.xlim(0,14)
  plt.ylim(0,300)
  plt.savefig('rc_curveB.png')

  plt.figure(12)
  plt.plot(r/kpc,e1,label='disc')
  plt.plot(r/kpc,e2,label='halo')
  plt.plot(r/kpc,e3,label='bulge')
  plt.legend()
  plt.savefig('pots.png')
  
  plt.figure(13)
  r1 = RhoDisk(r,0.)
  r3 = RhoBulge(r)
  plt.semilogy(r/kpc,(r1)/km,'b-',label='disc')
  plt.semilogy(r/kpc,(r3)/km,'g-',label='bulge')
  plt.ylim(1e-28,)
  plt.legend()
  plt.savefig('rhos.png')
  
  print 'Plotting swing amplitude as a fn of R'
  plt.figure(7)
  swingM = r * kappa_freq(r)**2 / (4. * 3.14159 * G * SigDisk(r))
  meanM  = np.mean(swingM)
  solarM = np.mean( swingM[(r>6.9*kpc) & (r<7.1*kpc)]  )
  plt.plot(r/kpc,swingM,'ko',markersize=1.0,label='$ m(R) $')
  plt.plot(((0.,15.)),meanM*np.ones(2),'r--',label='$<m> = $'+str(np.round(meanM,3)))
  plt.plot(((0.,15.)),solarM*np.ones(2),'g--',label='$m(R=7\\rm kpc) = $'+str(np.round(solarM,3)))
  "Limitted to rs"
  swingLim = 0.0*kpc
  rs       = r[r>swingLim]
  swingMrs = swingM[r>swingLim]
  meanMrs  = np.mean(swingMrs)
  if swingLim<0.:
    plt.plot(((0.,15.)),meanMrs*np.ones(2),'r-',label='$<m> = $'+str(np.round(meanMrs,3)))
  plt.plot(swingLim*np.ones(2)/kpc,((0,20)),'b--')
  "Disc only Kappa"
  swingM = r * kappa_freqDisc(r)**2 / (4. * 3.14159 * G * SigDisk(r))
  meanMD = np.mean(swingM)
  plt.plot(r/kpc,swingM,'go',markersize=1.0,label='$ m_d(R) $')
  plt.plot(((0.,15.)),meanMD*np.ones(2),'b--',label='$<m_d> = $'+str(np.round(meanMD,3)))

  plt.xlabel('R[kpc]')
  plt.ylabel('m')
  plt.ylim(0,10)
  plt.xlim(0,13)
  plt.legend()
  plt.savefig('swingAmp.png')
  print 'Mean Swing =',meanM
  print 'Mean Swing =',meanMrs,'(r>',swingLim/kpc,'kpc )'
  print 'Mean Swing =',meanMD,', kappaD only'
  print 'Solar Swing=',solarM

elif ans1=='p' and ansS=='d':
  print '------------------------------'
  print 'Setting disk star positions...'
  Npart =Npart_PD

  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  #DRAW THE RADII:
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  r  = np.arange(r1,r2,0.001*(r2-r1))
  
  cd1=np.zeros(len(r))
  pd1=np.zeros(len(r))
  s1 =np.zeros(len(r))
  i  =0
  for line in cd1:
      cd1[i] = cdfDisk(r[i])
      pd1[i] = pdfDisk(r[i])
      s1[i]  = SigDisk(r[i])
      i+=1
  
 
  if ansP=='p':
    Am =plt.figure(1,figsize=(14,6))
    A1=Am.add_subplot(1,2,1)
    A1.plot(r/kpc,cd1,label='cdf')
    A1.plot(r/kpc,pd1*kpc,label='pdf')
    A1.plot(r/kpc,s1, label='Sigma')
    plt.xlabel('R [kpc]')
    plt.xlim(r1/kpc,r2/kpc)
    plt.ylim(0,1.5)
    plt.legend(loc=2)
  
  
  rpred=randomdist(r,cd1,Npart)
  
  if rscatter>0:
    print 'Add some random offset to account for finite r array of',rscatter/kpc,' kpc'
    i=0
    while i<len(rpred):
      sig_posr = random.gauss(0.,rscatter)
      rpred[i] +=sig_posr
      #Ensure that this doesn't make stuff go r<0kpc:
      rpred[i] = abs(rpred[i])
      i+=1
  
  if ansP=='p':
    A1.hist(rpred/kpc,20,normed=True)
  
  phipred=np.ones(len(rpred))
  i=0
  for line in phipred:
    phipred[i]=2.*math.pi*random.random()
    i+=1
  
  dr   =0.05*(r2-r1)
  rlow =r1
  rhigh=r2
  rb1  =np.arange(rlow,rhigh,dr)
  rb2  =np.arange(rlow+dr,rhigh+dr,dr)
  s    =np.zeros(len(rb1))
  
  i=0
  for thing in s:
    s[i]= len(rpred[(rpred>rb1[i]) & (rpred<rb2[i])]) / (3.14159 * ((rb2[i]/kpc)**2-(rb1[i]/kpc)**2))
    i+=1
  scaler=60.0

  if ansP=='p':
    A1.errorbar(0.5*(rb1+rb2)/kpc,s/scaler,yerr=np.sqrt(s)/scaler,xerr=0.,fmt='k-o')
    A2=Am.add_subplot(1,2,2)
    A2.plot(rpred*np.cos(phipred)/kpc,rpred*np.sin(phipred)/kpc,'ko',markersize=2)
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    plt.savefig('DistR_D.png')
    
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  #DRAW THE VERTICAL DISTANCE:
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  z  = np.arange(z1,z2,0.001*(z2-z1))

  cd2=np.zeros(len(z))
  pd2=np.zeros(len(z))
  s2 =np.zeros(len(z))
  i  =0
  for line in cd2:
      cd2[i] = cdfDiskz(z[i])
      pd2[i] = pdfDiskz(z[i])
      s2[i]  = SigDiskz(z[i])
      i+=1
  if ansP=='p':
    Bm =plt.figure(2,figsize=(14,6))
    B1=Bm.add_subplot(1,2,1)
    B1.plot(z/kpc,cd2,label='cdf')
    B1.plot(z/kpc,0.5*pd2*kpc,label='pdf')   #Must be half to count for only using half dist.
    B1.plot(z/kpc,s2, label='Sigma')
    plt.xlabel('z [kpc]')
    plt.xlim(-z2/kpc,z2/kpc)
    plt.ylim(0,2.2)
    plt.legend(loc=2)

  
  "The z-dist is symmetric, so now randomly assign above/below z=0."  
  zpred=randomdist(z,cd2,Npart)
  i=0
  while i<len(zpred):
    if (random.random()>=0.5):
      zpred[i]=-zpred[i]
    i+=1
  
  if rscatter>0.:
    print 'Adding scatter to the z-dist of',rscatter/kpc,' kpc'
    i=0
    while i<len(zpred):
      sig_posz = random.gauss(0.,rscatter)
      zpred[i] +=sig_posz
      i+=1
  
  if ansP=='p':
    B1.hist(zpred/kpc,20,normed=True)
    B1.plot(-z/kpc,0.5*pd2*kpc,c='green')

  dz   =0.05*(2*z2)
  zlow =-z2
  zhigh=+z2
  zb1  =np.arange(zlow,zhigh,dz)
  zb2  =np.arange(zlow+dz,zhigh+dz,dz)
  s    =np.zeros(len(zb1))
  
  i=0
  for thing in s:
    s[i]= len(zpred[(zpred>zb1[i]) & (zpred<zb2[i])]) / (3.14159 * ((zb2[i]/kpc)-(zb1[i]/kpc)))
    i+=1
  scaler=3800.0

  if ansP=='p':
    B1.errorbar(0.5*(zb1+zb2)/kpc,s/scaler,yerr=np.sqrt(s)/scaler,xerr=0.,fmt='k-o')
    B2=Bm.add_subplot(1,2,2)
    B2.plot(rpred*np.cos(phipred)/kpc,zpred/kpc,'ko',markersize=2)
    plt.xlabel('x [kpc]')
    plt.ylabel('z [kpc]')
    plt.savefig('DistZ_D.png')

  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  #Save Outputs
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  np.savez('positionsD_'+str(Npart),r=rpred,phi=phipred,z=zpred)
  exit()

elif ans1=='p' and ansS=='b':
  print '-------------------------------'
  print 'Setting bulge star positions...'
  Npart =Npart_PB
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  #DRAW THE RADII:
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  r  = np.arange(r1,r2b,0.001*(r2b-r1))
  'r2 is really only for full disk stars, but prob will be so low beyond 5kpc none should be drawn'
  cd1=np.zeros(len(r))
  pd1=np.zeros(len(r))
  s1 =np.zeros(len(r))
  i  =0
  for line in cd1:
      cd1[i] = cdfBulge(r[i])
      pd1[i] = pdfBulge(r[i])
      #s1[i]  = RhoBulge(r[i]) * r[i]*2.0*3.14159 / 10.    #A rho not a sigma
      s1[i]  = RhoBulge(r[i])*kpc
      i+=1
  
  if ansP=='p':
    Am =plt.figure(1,figsize=(14,6))
    A1=Am.add_subplot(1,2,1)
    A1.plot(r/kpc,cd1,label='cdf')
    rhoscaler = kpc*kpc/2.5
    A1.plot(r/kpc,pd1*rhoscaler,label='pdf')
    A1.semilogy(r/kpc,s1, label='Rho')
    plt.xlabel('R [kpc]')
    plt.xlim(r1/kpc,r2b/kpc)
    plt.ylim(0,1.5)
    plt.legend(loc=1)
  
  rpred=randomdist(r,cd1,Npart)
  
  if rscatter>0:
    print 'Add some random offset to account for finite r array'
    print 'Scattering by',rscatter/kpc,'[kpc]'
    i=0
    while i<len(rpred):
      sig_posr = random.gauss(0.,rscatter)
      rpred[i] +=sig_posr
      #Ensure that this doesn't make stuff go r<0kpc:
      rpred[i] = abs(rpred[i])
      i+=1
  
  if ansP=='p':
    A1.hist(rpred/kpc,500,normed=True)
  
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  #DRAW THE ANGLES (phi & theta) FROM RAN.U.DIST.:
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  #if ansP=='p':
  #  B1.hist(zpred/kpc,20,normed=True)
  #  B1.plot(-z/kpc,0.5*pd2*kpc,c='green')
 
  phipred=np.zeros(len(rpred))
  i=0
  for line in phipred:
    phipred[i]=2.*math.pi*random.random()
    i+=1
 
  thetapred=np.zeros(len(rpred))
  i=0
  for line in thetapred:
    #The below can't just be random in phi, there needs to be some weighting:
    #see: http://mathworld.wolfram.com/SpherePointPicking.html
    thetapred[i]=np.arccos(2.*random.random()-1.)   #math.pi*random.random()
    i+=1
  
  dr   =0.01*(r2b-r1)
  rlow =r1
  rhigh=r2b
  rb1  =np.arange(rlow,rhigh,dr)
  rb2  =np.arange(rlow+dr,rhigh+dr,dr)
  s    =np.zeros(len(rb1))
  
  i=0
  for thing in s:
    s[i]= len(rpred[(rpred>rb1[i]) & (rpred<rb2[i])]) / (4./3. * 3.14159 * ((rb2[i]/kpc)**3-(rb1[i]/kpc)**3))
    i+=1
  scaler=150.0

  if ansP=='p':
    x=rpred*np.cos(phipred)*np.sin(thetapred)
    y=rpred*np.sin(phipred)*np.sin(thetapred)
    z=rpred                *np.cos(thetapred)
    
    #A1.errorbar(0.5*(rb1+rb2)/kpc,s/scaler,yerr=np.sqrt(s)/scaler,xerr=0.,fmt='k-o')
    A1.plot(0.5*(rb1+rb2)/kpc,s/scaler,'k-o')
    #plt.ylim(5e-6,1e4)
    plt.ylim(1e-3,1e4)
    plt.xlim(0.,10.)
    
    ig=0
    for line in s:
      plt.annotate(len(rpred[(rpred>rb1[ig]) & (rpred<rb2[ig])]),xy=((0.5*(rb1+rb2)/kpc)[ig],(s/scaler)[ig]))
      ig+=1
    A2=Am.add_subplot(1,2,2)
    A2.plot(x/kpc,y/kpc,'ko',markersize=2)
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    plt.xlim(-15,15)
    plt.ylim(-15,15)
    plt.savefig('DistR_B.png')

   
    B=plt.figure(3,figsize=(14,6))
    B1=B.add_subplot(1,2,1)
    B1.plot(x/kpc,y/kpc,'ko',markersize=2)
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    plt.xlim(-15,15)
    plt.ylim(-15,15)
    B2=B.add_subplot(1,2,2)
    B2.plot(x/kpc,z/kpc,'ko',markersize=2)
    plt.xlabel('x [kpc]')
    plt.ylabel('z [kpc]')
    plt.xlim(-15,15)
    plt.ylim(-15,15)
    plt.savefig('xyzPlotB.png')
    
    plt.figure(4)
    plt.hist(x/kpc,30,normed=True,alpha=0.5,label='x')
    plt.hist(y/kpc,30,normed=True,alpha=0.5,label='y')
    plt.hist(z/kpc,30,normed=True,alpha=0.5,label='z')
    plt.legend()
    plt.savefig('xyzHistB.png')
  
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  #Save Outputs
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  np.savez('positionsB_'+str(Npart),r=rpred,phi=phipred,theta=thetapred)
  exit()

elif ans1=='v' and ansS=='d':
  print '------------------------------------------'
  print 'Setting disk star velocities and masses...'
  
  Npart=Npart_VD    #Work with just a set for velocity setting (Npart_V<=Npart_P).
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  #DRAW THE RADII:
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  inputz=np.load('positionsD_'+str(Npart_PD)+'.npz')
  r   = inputz['r'][0:Npart]
  phi = inputz['phi'][0:Npart]
  z   = inputz['z'][0:Npart]
  x   = r*np.cos(phi)
  y   = r*np.sin(phi)
  if ansP=='p':
    plt.figure(1)
    plt.plot(x/kpc,y/kpc,'ko',markersize=2)
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    plt.savefig('xyPlotD.png')
  
  '''
  First we need the radom components of the radial, azimuthal, and verical velocities.
  These a taken from Gaussians centred on 0, 0 and v_phi respectively.
  '''
  sig_z  =np.zeros(len(r))
  sig_r  =np.zeros(len(r))
  sig_phi=np.zeros(len(r))
  v_phi  =np.zeros(len(r))

  "Begin with assigning sig_z"
  sig_z2 = math.pi*G*z0*SigDisk(r)
  i=0
  while i<len(sig_z):
    sig_z[i] = random.gauss(0.,np.sqrt(sig_z2[i]))
    i+=1
  if ansP=='p':
    plt.figure(2,figsize=(14,6))
    plt.subplot(1,2,1)
    plt.scatter(x/kpc,y/kpc,c=sig_z/km)
    cbar=plt.colorbar()
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    plt.subplot(1,2,2)
    plt.plot(r/kpc,sig_z/km,'ko',alpha=0.2)
    plt.xlabel('r [kpc]')
    plt.ylabel('$V_z \\rm [km/s]$',fontsize=16)
    plt.savefig('sigZ_D.png')
  
  "Then assign sig_r"
  SigVConst =ToomreConstFn(r)
  print 'The Sig_v constant is:',np.sqrt(SigVConst)/km,'km/s'
  sig_r2 = SigVConst * np.exp(-r/r0)
  i=0
  while i<len(sig_r):
    sig_r[i] = random.gauss(0.,np.sqrt(sig_r2[i]))
    i+=1  
  if ansP=='p':
    plt.figure(3,figsize=(14,6))
    plt.subplot(1,2,1)
    plt.scatter(x/kpc,y/kpc,c=sig_r/km)
    cbar=plt.colorbar()
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    plt.subplot(1,2,2)
    plt.plot(r/kpc,sig_r/km,'ko',alpha=0.2)
    plt.xlabel('r [kpc]')
    plt.ylabel('$V_r \\rm [km/s]$',fontsize=16)
    plt.savefig('sigR_D.png')
  
  "Then sig_phi, and finally the azimuthal streaming velocity, v_phi"
  sig_phi2 = sig_r2 * kappa_freq(r)**2 / (4.*omega_freq(r)**2)
  v_phi2   = sig_r2 - sig_phi2 - sig_r2*r/r0 + (Vc(r))**2
  i=0
  while i<len(sig_phi):
    sig_phi[i] = random.gauss(np.sqrt(v_phi2[i]),np.sqrt(sig_phi2[i]))
    i+=1  
  if ansP=='p':
    plt.figure(4,figsize=(14,6))
    plt.subplot(1,2,1)
    plt.scatter(x/kpc,y/kpc,c=sig_phi/km)
    cbar=plt.colorbar()
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    plt.subplot(1,2,2)
    plt.plot(r/kpc,sig_phi/km,'ko',alpha=0.2)
    plt.xlabel('r [kpc]')
    plt.ylabel('$\sigma_\phi \\rm [km/s]$',fontsize=16)
    plt.savefig('sigR_D.png')
  
  sig_phi = -sig_phi
  vx = sig_r*np.cos(phi) - sig_phi*np.sin(phi)
  vy = sig_r*np.sin(phi) + sig_phi*np.cos(phi)
  vz = sig_z

  MWgas = Mg #1.0E9*Mo  From Wada paper stuff.
  MWstr = Md
  Nstr  = float(StarToGas)*float(Npart)
  Ngas  = float(Npart)-float(Nstr)
  print 'Ngas :',Nstr
  print 'Nstar:',Ngas

  if Ngas>0:
    mi_g  = MWgas/Ngas
  else:
    mi_g  = 0.
  mi_s  = MWstr/Nstr
  p_g   =0
  p_s   =10
  
  "Set masses, in kg:"
  mass  =np.zeros(len(r))
  phase =np.zeros(len(mass))          #0 for gas, 10 for stars
  i=0
  while i<len(mass):  #Auto fill as stars
     mass[i] = mi_s
     phase[i]= p_s
     if i>=Nstr:
       mass[i] =mi_g
       phase[i]=p_g
     i+=1
  
  "Disc only Kappa"
  swingM = r * kappa_freq(r)**2 / (4. * 3.14159 * G * SigDisk(r))
  meanM  = np.mean(swingM)
  swingMD= r * kappa_freqDisc(r)**2 / (4. * 3.14159 * G * SigDisk(r))
  meanMD = np.mean(swingMD)
  print 'Mean Swing =',meanM
  print 'Mean Swing =',meanMD,', kappaD only'

  if ansP=='p':
    plt.figure(6)
    vr = np.sqrt(vx*vx+vy*vy+vz*vz)
    plt.plot(r/kpc,vr/1000.,'ko',markersize=4.0)
    vtester =np.arange(0.,np.max(r)/kpc,0.02)*kpc
    vcirc   =Vc(vtester)
    plt.plot(vtester/kpc,vcirc/km,linewidth=3)
    plt.xlabel('R[kpc]')
    plt.ylabel('Vr[km/s]')
    plt.savefig('rc_setD.png')
    
    print 'Plotting swing amplitude as a fn of R'
    plt.figure(7)
    plt.plot(r/kpc,swingM,'ko',markersize=1.0,label='$ m(R) $')
    plt.plot(((0.,15.)),meanM*np.ones(2),'r--',label='$ m_{avg,i} = $'+str(np.round(meanM,3)))
    plt.plot(r/kpc,swingMD,'go',markersize=1.0,label='$ m_d(R) $')
    plt.plot(((0.,15.)),meanMD*np.ones(2),'b--',label='$<m_d> = $'+str(np.round(meanMD,3)))
    plt.xlabel('R[kpc]')
    plt.ylabel('m')
    plt.ylim(0,10)
    plt.xlim(0,13)
    plt.legend()
    plt.savefig('swingAmpPost.png')
 
  outd=open('asciifile_D','w')
  outg=open('asciifile_G','w')
  iout=0
  while iout<Nstr:
    outd.write(str(x[iout])+' '+str(y[iout])+' '+str(z[iout])+' '+str(mass[iout])+' '+str(vx[iout])+' '+str(vy[iout])+' '+str(vz[iout])+' '+str(phase[iout])+'\n')
    iout+=1
  outd.close()
  while iout<Nstr+Ngas:
    outg.write(str(x[iout])+' '+str(y[iout])+' '+str(z[iout])+' '+str(mass[iout])+' '+str(vx[iout])+' '+str(vy[iout])+' '+str(vz[iout])+' '+str(phase[iout])+'\n')
    iout+=1
  outg.close()
  exit()

elif ans1=='v' and ansS=='b':
  print '-------------------------------------------'
  print 'Setting bulge star velocities and masses...'
  Npart=Npart_VB

  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  #DRAW THE RADII:
  "<><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
  inputz   = np.load('positionsB_'+str(Npart_PB)+'.npz')
  Ndisc    = Npart_VD
  diskfile = 'asciifile_D'
  print 'reading from:',inputz,' & ',diskfile,'[',Ndisc,']'
  r    = inputz['r'][0:Npart]
  phi  = inputz['phi'][0:Npart]
  theta= inputz['theta'][0:Npart]
  x   = r*np.cos(phi)*np.sin(theta)
  y   = r*np.sin(phi)*np.sin(theta)
  z   = r*            np.cos(theta)
  if ansP=='p':
    plt.figure(1)
    plt.subplot(221)
    plt.plot(x/kpc,y/kpc,'ko',markersize=2)
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    plt.xlim(-13,13)
    plt.ylim(-13,13)
    plt.subplot(222)
    plt.plot(x/kpc,z/kpc,'ko',markersize=2)
    plt.xlabel('x [kpc]')
    plt.ylabel('z [kpc]')
    plt.xlim(-13,13)
    plt.ylim(-13,13)
    plt.subplot(223)
    plt.plot(y/kpc,z/kpc,'ko',markersize=2)
    plt.xlabel('y [kpc]')
    plt.ylabel('z [kpc]')
    plt.xlim(-13,13)
    plt.ylim(-13,13)
    plt.savefig('xyzPlotB_in.png')
  '''
  First we need the radom components of the radial, azimuthal, and verical velocities.
  These a taken from Gaussians centred on 0, 0 and v_phi respectively.
  '''
  sig_r =np.zeros(len(r))
  vpred =np.zeros(len(r))


  "Begin with assigning r velocity dispersion"
  "We assume sperhical isotropy so v_r = v_phi = v_theta"
  rlim = 100.*kpc
  i    = 0
  if ansJ=='m':
    print ':IMPORTANT::IMPORTANT::IMPORTANT::IMPORTANT::IMPORTANT::IMPORTANT:'
    print '--------PLEASE REMEMBER TO MAKE THE DISC BEFORE THE BULGE---------'
    print '-------PLEASE CHECK DISC STRARtoGAS RATIO BEFORE THE BULGE--------'
    print ':IMPORTANT::IMPORTANT::IMPORTANT::IMPORTANT::IMPORTANT::IMPORTANT:'
    print Ndisc
    diskarray = readindisk(diskfile,Ndisc)   #Returns an array of [r,m] in spherical polars
    
    #Should call the disk counter FIRST here:
    drint   = rlim/10000.#was 100000. for good res.
    r_int   = np.arange(0.,rlim,drint)
    MassinDisc_Ar  =  np.zeros(len(r_int))
    BJ_eq          =  np.zeros(len(r_int))
    i_mass=0
    while i_mass<len(r_int):
      MassinDisc_Ar[i_mass] = COUNTtheDISK(r_int[i_mass])
      BJ_eq[i_mass]  =  BulgeJeansEq_M(r_int[i_mass])
      i_mass+=1
    PlotMasses(r)

    while i<len(sig_r):
      working_r_Ar  =  r_int[ r_int>r[i] ]
      working_BJ_Ar =  BJ_eq[ r_int>r[i] ]
      inout = np.sqrt(simp(working_BJ_Ar,x=working_r_Ar)/RhoBulge(r[i]))
      sig_r[i] =inout
      i+=1


  if ansJ=='p':
    while i<len(sig_r):
      dxint   =(rlim-r[i])/10000.
      x_int   =np.arange(r[i],rlim,dxint)
      intg    =BulgeJeansEq_P(x_int)

      inout  = np.sqrt(  dxint*intg.sum()/RhoBulge(r[i]) )
      #inout = np.sqrt(trap(intg,x=x_int)/RhoBulge(r[i]))
      #inout = np.sqrt(spiq(BulgeJeansEq_P,r[i],rlim)/RhoBulge(r[i]))[0]
      
      sig_r[i]=inout
      i+=1
    plt.figure(2)
    plt.plot(r/kpc,sig_r/km,'ko')
    plt.xlabel('r [kpc]')
    plt.ylabel('sigma_r [km/s]')
    plt.ylim(0,160.)
    plt.savefig('sigR_B.png')

  
  '''
  How to use scipy.stats.maxwell, a is the mean velocity factor (a^2 = kB*T/m)
  x0 is the centring of the function. 
  rv = maxwell(x0,a)
  rv.cdf(x)  or  rv.pdf(x)
  '''
  "Make a dummy array to loop through to find the velocties form random variables."
  var = np.arange(0.0,700.0*km,10.)
  i    = 0
  while i<len(sig_r):
    #Check to see if bulge velocity breaks the systems escape velocity
    e1 = PotDisk('p',r[i])
    e2 = PotHalo('p',r[i])
    e3 = PotBulge('p',r[i])
    VescTemp = vescFac * np.sqrt( -2.*(e1+e2+e3) )

    MB_cdf   = mbd(0,sig_r[i]).cdf(var)
    thisVel  = 1e30   #Arbitrary large velocity.
    while thisVel > VescTemp:
      #Draw until new vel is lower than escape vel.
      thisVel = randomdist(var,MB_cdf,1)
    vpred[i] = thisVel
    i+=1
  if ansP=='p':
    plt.figure(3)
    plt.hist(vpred/km,20,normed=True)
    vmean= np.mean(sig_r)/km
    vark = var/km
    MaxB=4*math.pi*(1.0/(2.0*math.pi*vmean*vmean))**1.5 *vark*vark*np.exp(-vark*vark/(2.0*vmean*vmean))
    plt.plot(var/km,MaxB,'k--',linewidth=4)
    vmean= np.min(sig_r)/km
    vark = var/km
    MaxB=4*math.pi*(1.0/(2.0*math.pi*vmean*vmean))**1.5 *vark*vark*np.exp(-vark*vark/(2.0*vmean*vmean))
    plt.plot(var/km,MaxB,'r--',linewidth=4)
    vmean= np.max(sig_r)/km
    vark = var/km
    MaxB=4*math.pi*(1.0/(2.0*math.pi*vmean*vmean))**1.5 *vark*vark*np.exp(-vark*vark/(2.0*vmean*vmean))
    plt.plot(var/km,MaxB,'g--',linewidth=4)
    plt.xlabel('|V| [km/s]')
    plt.savefig('MBdists.png')

    plt.figure(4,figsize=(14,6))
    plt.subplot(1,2,1)
    plt.scatter(x/kpc,y/kpc,c=vpred/km)
    cbar=plt.colorbar()
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    plt.subplot(1,2,2)
    plt.plot(r/kpc,vpred/km,'ko',alpha=0.2)
    plt.xlabel('r [kpc]')
    plt.ylabel('$V_r \\rm [km/s]$',fontsize=16)   
    plt.savefig('sigTotB.png')
    
  
  #Now (assuming isotropy) pluck a random theta and phi this time to define the velocity vector...
  ij=0
  vx=np.zeros(len(vpred))
  vy=np.zeros(len(vpred))
  vz=np.zeros(len(vpred))
  print 'Assigning velocities and streaming fraction'
  while ij<len(vpred):
    phi_vel  =2.*math.pi*random.random()
    theta_vel=np.arccos(2.*random.random()-1.)
    vx[ij] = vpred[ij]*np.sin(theta_vel)*np.cos(phi_vel)
    vy[ij] = vpred[ij]*np.sin(theta_vel)*np.sin(phi_vel)
    vz[ij] = vpred[ij]*np.cos(theta_vel)
    
    #Streaming (copied verbatim from MKKD95 in NEMO)
    #This is done in cylindrical polars, whereas above is in spherical.
    #Pluck a random number, if it beats the streaming fraction, then assign those with positive rotation.
    r_xy = np.sqrt(x[ij]*x[ij] + y[ij]*y[ij])
    vrtemp = (+vx[ij]*x[ij] + vy[ij]*y[ij])/r_xy
    vptemp = (-vx[ij]*y[ij] + vy[ij]*x[ij])/r_xy
    rantemp= random.random()
    if rantemp<stream:
      vptemp = -abs(vptemp)
    else:
      vptemp = +abs(vptemp)
    vx[ij] = (vrtemp * x[ij] - vptemp * y[ij])/r_xy
    vy[ij] = (vrtemp * y[ij] + vptemp * x[ij])/r_xy
    
    ij+=1
  

  "The stellar mass in the bulge"
  print 'Assigning masses'
  MWgas = Mg 
  MWstr = Mb 
  Nstr  = Npart        #In the bulge we only use star particles
  Ngas  = Npart-Nstr

  if Ngas>0:
    mi_g  = MWgas/Ngas
  else:
    mi_g  = 0.

  if Nstr>0:
    mi_s  = MWstr/Nstr
  else:
    mi_s  = 0.
  "Phase flags"
  p_g   =0
  p_s   =10
  
  "Set masses, in kg:"
  mass  =np.zeros(len(r))
  phase =np.zeros(len(mass))          #0 for gas, 10 for stars
  i=0
  while i<len(mass):  #Auto fill as stars
     mass[i] = mi_s
     phase[i]= p_s
     if i>Nstr:
       mass[i] =mi_g
       phase[i]=p_g
     i+=1
  
  if ansP=='p':
    plt.figure(6)
    vr = np.sqrt(vx*vx+vy*vy+vz*vz)
    plt.plot(r/kpc,vr/km,'ko',markersize=4.0)
    plt.xlabel('R[kpc]')
    plt.ylabel('Vr[km/s]')
    vtester =np.arange(0.,np.max(r)/kpc,0.1)*kpc
    vcirc   =Vc(vtester)
    plt.plot(vtester/kpc,vcirc/km)
    plt.savefig('rc_setB.png')

  outp=open('asciifile_B','w')
  iout=0
  while iout<Npart:
    outp.write(str(x[iout])+' '+str(y[iout])+' '+str(z[iout])+' '+str(mass[iout])+' '+str(vx[iout])+' '+str(vy[iout])+' '+str(vz[iout])+' '+str(phase[iout])+'\n')
    iout+=1
  outp.close()
  print "Done setting masses and vels..."
  print '-------------------------------------------'
  #print "Now use:"
  #print "cat asciifile_B asciifile_D > asciifile_T"
  #print "to concatenate the disk and bulge together."
  exit()

elif ans1=='f':
  print 'Last step, so write out the file needed for phantom to setup.'
  Nstr  = int(float(StarToGas)*float(Npart_PD))
  Ngas  = int(float(Npart_PD)-float(Nstr))
  Nbulge= int(Npart_PB)
  Nhalo = 0
  outnom = 'galsetic.txt'
  outo=open(outnom,'w')
  iout=0
  outo.write('Ngas'+' '+str(Ngas)+'\n')
  outo.write('Nstar'+' '+str(Nstr)+'\n')
  outo.write('Nbulge'+' '+str(Nbulge)+'\n')
  outo.write('Ndark'+' '+str(Nhalo)+'\n')
  outo.close()
  print '<<<<<<<<<<<<<<DONE>>>>>>>>>>>>>>>'
  print 'Now place',outnom,'and the asciifiles in the directory of phantomsetup.'


else:
  print "Nothing good can come of this..."
  print "you must be reading stuff in..."
  npart=10000
  x =np.zeros(npart)
  y =np.zeros(npart)
  z =np.zeros(npart)
  vx=np.zeros(npart)
  vy=np.zeros(npart)
  vz=np.zeros(npart)
  Readin = open('asciifile_b')
  i=0
  for line in Readin:
	  x[i]=float(line.split()[0])/kpc
	  y[i]=float(line.split()[1])/kpc
	  z[i]=float(line.split()[2])/kpc
	  vx[i]=float(line.split()[4])*1e3
	  vy[i]=float(line.split()[5])*1e3
	  vz[i]=float(line.split()[6])*1e3
	  i+=1
  
  
  r=np.sqrt(x*x+y*y+z*z)
  theta=np.arctan(y/x)
  phi  =np.arccos(z/r)
  vr  =vx*x/r + vy*y/r + vz*z/r
  vphi=vx*(-y/np.sqrt(x*x+y*y)) + vy*(x/np.sqrt(x*x+y*y))
  vtheta=vx*(x*z/np.sqrt(x*x+y*y))/r + vy*(y*z/np.sqrt(x*x+y*y))/r - vy*((x*x+y*y)/np.sqrt(x*x+y*y))/r

  plt.figure(1)
  plt.plot(theta,phi,'ko',markersize=1.0)
  plt.figure(2)
  plt.plot(y,vy,'go',markersize=3)
  plt.plot(z,vz,'bo',markersize=3)
  plt.plot(x,vx,'ro',markersize=3)

  plt.show() 

  exit()

  

  
